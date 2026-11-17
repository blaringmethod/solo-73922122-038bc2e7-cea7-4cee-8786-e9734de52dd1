"""FastAPI application: versioned VCF upload + indexed interval queries,
plus a versioned gene catalog (BED) and indexed gene annotation of variants.

Layout of the shared data directory:

    <data_dir>/meta.sqlite3          metadata, current-version pointers, ref pins
    <data_dir>/versions/<uuid>/      published VCF file sets (variants.vcf.gz + .tbi)
    <data_dir>/catalogs/<uuid>/      published catalog file sets (genes.bed.gz + .tbi)
    <data_dir>/staging/<instance>/<uuid>/   in-progress builds

Publishing is atomic: the file set is completed in staging, renamed into
versions/ (or catalogs/), and only then the current-version pointer is
flipped inside one SQLite transaction. Readers pin a version before opening
its files, so an in-flight query keeps its (version, data, index) trio
consistent while a replacement is published; annotation queries pin the
(VCF version, catalog) pair in a single transaction, so they always see one
moment's snapshot of both. The reaper frees old files once no pin remains.
"""
from __future__ import annotations

import asyncio
import os
import shutil
import time
import uuid
from contextlib import asynccontextmanager
from dataclasses import dataclass

import pysam
from fastapi import FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import RedirectResponse

from . import bedbuild
from .db import Conflict, Store
from .vcfbuild import GZ_NAME, VCFError, validate_and_build

DEFAULT_MAX_UPLOAD = 10 * 1024 * 1024  # 10 MiB


@dataclass
class Config:
    data_dir: str
    instance_id: str
    max_upload_bytes: int = DEFAULT_MAX_UPLOAD
    reaper_interval: float = 2.0      # seconds between reaper runs (0 disables)
    heartbeat_interval: float = 2.0   # seconds between heartbeats (0 disables)
    dead_after: float = 15.0          # instance considered dead after this silence
    orphan_grace: float = 60.0        # min age before orphan dirs are swept


def _config_from_env() -> Config:
    e = os.environ
    return Config(
        data_dir=e.get("VCF_DATA_DIR", "/data"),
        instance_id=e.get("VCF_INSTANCE_ID") or uuid.uuid4().hex,
        max_upload_bytes=int(e.get("VCF_MAX_UPLOAD_BYTES", DEFAULT_MAX_UPLOAD)),
        reaper_interval=float(e.get("VCF_REAPER_INTERVAL", "2")),
        heartbeat_interval=float(e.get("VCF_HEARTBEAT_INTERVAL", "2")),
        dead_after=float(e.get("VCF_DEAD_AFTER", "15")),
        orphan_grace=float(e.get("VCF_ORPHAN_GRACE", "60")),
    )


def create_app(config: Config | None = None) -> FastAPI:
    cfg = config or _config_from_env()
    store = Store(os.path.join(cfg.data_dir, "meta.sqlite3"))
    versions_dir = os.path.join(cfg.data_dir, "versions")
    catalogs_dir = os.path.join(cfg.data_dir, "catalogs")
    staging_base = os.path.join(cfg.data_dir, "staging")
    staging_root = os.path.join(staging_base, cfg.instance_id)
    os.makedirs(versions_dir, exist_ok=True)
    os.makedirs(catalogs_dir, exist_ok=True)
    # Leftover builds of a previous incarnation of this instance can never be
    # published any more; drop them before serving.
    shutil.rmtree(staging_root, ignore_errors=True)
    os.makedirs(staging_root, exist_ok=True)

    def reap_once() -> None:
        now = time.time()
        # 1) DB: drop dead instances' pins, delete reclaimable version/catalog rows.
        for rel in store.reap_metadata(now, cfg.dead_after):
            shutil.rmtree(os.path.join(cfg.data_dir, rel), ignore_errors=True)
        # 2) FS: sweep version/catalog dirs with no DB row (crashed publish/reaper).
        cutoff = now - cfg.orphan_grace
        for base_dir, known in (
            (versions_dir, store.all_version_paths()),
            (catalogs_dir, store.all_catalog_paths()),
        ):
            for name in os.listdir(base_dir):
                path = os.path.join(base_dir, name)
                rel = os.path.relpath(path, cfg.data_dir)
                try:
                    old = os.path.getmtime(path) < cutoff
                except OSError:
                    continue
                if rel not in known and old:
                    shutil.rmtree(path, ignore_errors=True)
        # 3) FS: sweep staging dirs of dead instances.
        live = store.live_instances(now, cfg.dead_after)
        for inst in os.listdir(staging_base):
            path = os.path.join(staging_base, inst)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if inst not in live and old:
                shutil.rmtree(path, ignore_errors=True)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        store.heartbeat(cfg.instance_id, time.time())
        stopped = asyncio.Event()

        async def heartbeat_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.heartbeat_interval)
                try:
                    store.heartbeat(cfg.instance_id, time.time())
                except Exception:
                    pass

        async def reaper_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.reaper_interval)
                try:
                    await asyncio.to_thread(reap_once)
                except Exception:
                    pass

        tasks = []
        if cfg.heartbeat_interval > 0:
            tasks.append(asyncio.create_task(heartbeat_loop()))
        if cfg.reaper_interval > 0:
            tasks.append(asyncio.create_task(reaper_loop()))
        try:
            yield
        finally:
            stopped.set()
            for t in tasks:
                t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            try:
                store.unregister(cfg.instance_id)
            except Exception:
                pass

    app = FastAPI(
        title="Versioned VCF variant store",
        version="1.1.0",
        lifespan=lifespan,
        description=(
            "Upload single-chromosome biallelic VCFs (max 10 MiB, no sample "
            "columns, REF/ALT in ACGT); records are sorted, BGZF-compressed and "
            "Tabix-indexed. Query 1-based closed intervals; every record whose "
            "reference span [POS, POS+len(REF)-1] intersects the interval is "
            "returned via the index. Uploads are optimistic-concurrency "
            "controlled through `base_version`.\n\n"
            "A versioned gene catalog complements the variants: upload a "
            "single-chromosome 4-column BED (0-based half-open intervals, at "
            "most 100 entries, unique gene names) and annotate the published "
            "variants with every gene whose region intersects a record's REF "
            "span. Annotation queries pin the (VCF version, catalog) pair at "
            "one moment and report both versions; while no usable catalog is "
            "published they report `pending_update` instead of stale gene "
            "names."
        ),
    )
    app.state.cfg = cfg
    app.state.store = store
    app.state.reap_once = reap_once

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse("/docs")

    @app.post("/vcf", status_code=201)
    async def upload_vcf(
        file: UploadFile = File(...),
        base_version: int = Query(..., ge=0, description="Version the client last saw (0 when empty)"),
    ):
        staging = os.path.join(staging_root, uuid.uuid4().hex)
        os.makedirs(staging)
        try:
            upload_path = os.path.join(staging, "upload.vcf")
            size = 0
            with open(upload_path, "wb") as out:
                while True:
                    chunk = await file.read(1 << 20)
                    if not chunk:
                        break
                    size += len(chunk)
                    if size > cfg.max_upload_bytes:
                        raise HTTPException(
                            413, f"upload exceeds the {cfg.max_upload_bytes}-byte limit"
                        )
                    out.write(chunk)
            try:
                chrom, n_records = await asyncio.to_thread(
                    validate_and_build, upload_path, staging
                )
            except VCFError as exc:
                raise HTTPException(400, f"invalid VCF: {exc}")
            os.remove(upload_path)

            # File set is complete; move it into place, then flip the pointer.
            final_rel = os.path.join("versions", uuid.uuid4().hex)
            os.rename(staging, os.path.join(cfg.data_dir, final_rel))
            try:
                new_version = await asyncio.to_thread(
                    store.publish, base_version, chrom, final_rel, n_records, time.time()
                )
            except Conflict as c:
                shutil.rmtree(os.path.join(cfg.data_dir, final_rel), ignore_errors=True)
                raise HTTPException(
                    409,
                    f"base_version {base_version} does not match current version "
                    f"{c.current}; re-read and retry",
                )
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)
        return {"version": new_version, "chrom": chrom, "n_records": n_records}

    @app.get("/variants")
    def get_variants(
        chrom: str = Query(..., description="Chromosome name as in the VCF"),
        start: int = Query(..., ge=1, description="1-based interval start (inclusive)"),
        end: int = Query(..., ge=1, description="1-based interval end (inclusive)"),
    ):
        if end < start:
            raise HTTPException(400, "end must be >= start")
        pinned = store.pin(cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        version, vchrom, rel, ref_id = pinned
        try:
            if chrom != vchrom:
                raise HTTPException(
                    400, f"the published dataset contains chromosome '{vchrom}' only"
                )
            gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
            try:
                tb = pysam.TabixFile(gz_path)
                try:
                    # 1-based closed [start, end] -> 0-based half-open fetch.
                    lines = tb.fetch(vchrom, start - 1, end)
                    records = []
                    for ln in lines:
                        f = ln.split("\t", 4)
                        pos = int(f[1])
                        # Exact span check on top of the index iterator.
                        if pos <= end and pos + len(f[3]) - 1 >= start:
                            records.append(ln)
                finally:
                    tb.close()
            except HTTPException:
                raise
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
            return {
                "version": version,
                "chrom": vchrom,
                "start": start,
                "end": end,
                "count": len(records),
                "records": records,
            }
        finally:
            store.unpin(ref_id)

    @app.post("/catalog", status_code=201)
    async def upload_catalog(
        file: UploadFile = File(...),
        base_version: int = Query(..., ge=0, description="Catalog version the client last saw (0 when empty)"),
    ):
        staging = os.path.join(staging_root, uuid.uuid4().hex)
        os.makedirs(staging)
        try:
            upload_path = os.path.join(staging, "upload.bed")
            size = 0
            with open(upload_path, "wb") as out:
                while True:
                    chunk = await file.read(1 << 20)
                    if not chunk:
                        break
                    size += len(chunk)
                    if size > cfg.max_upload_bytes:
                        raise HTTPException(
                            413, f"upload exceeds the {cfg.max_upload_bytes}-byte limit"
                        )
                    out.write(chunk)
            try:
                chrom, n_genes = await asyncio.to_thread(
                    bedbuild.validate_and_build, upload_path, staging
                )
            except bedbuild.BEDError as exc:
                raise HTTPException(400, f"invalid BED: {exc}")
            os.remove(upload_path)

            # File set is complete; move it into place, then flip the pointer.
            final_rel = os.path.join("catalogs", uuid.uuid4().hex)
            os.rename(staging, os.path.join(cfg.data_dir, final_rel))
            try:
                new_catalog = await asyncio.to_thread(
                    store.publish_catalog, base_version, chrom, final_rel, n_genes, time.time()
                )
            except Conflict as c:
                shutil.rmtree(os.path.join(cfg.data_dir, final_rel), ignore_errors=True)
                raise HTTPException(
                    409,
                    f"base_version {base_version} does not match current catalog "
                    f"version {c.current}; re-read and retry",
                )
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)
        return {"catalog_version": new_catalog, "chrom": chrom, "n_genes": n_genes}

    @app.get("/annotations")
    def get_annotations(
        chrom: str = Query(..., description="Chromosome name as in the VCF"),
        start: int = Query(..., ge=1, description="1-based interval start (inclusive)"),
        end: int = Query(..., ge=1, description="1-based interval end (inclusive)"),
    ):
        if end < start:
            raise HTTPException(400, "end must be >= start")
        # Pin the (VCF, catalog) pair in one transaction: the whole query then
        # reads one moment's files, index and catalog, however many VCF
        # corrections or catalog replacements are published meanwhile.
        pinned = store.pin_pair(cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        version, vchrom, vrel, vref_id, cat, cchrom, crel, cref_id = pinned
        try:
            if chrom != vchrom:
                raise HTTPException(
                    400, f"the published dataset contains chromosome '{vchrom}' only"
                )
            if cat is None:
                return {
                    "status": "pending_update",
                    "reason": "no_catalog",
                    "detail": "no gene catalog has been published yet; "
                              "upload a BED catalog to enable annotation",
                    "vcf_version": version,
                    "catalog_version": None,
                }
            if cchrom != vchrom:
                # Never annotate with a catalog built for another chromosome.
                return {
                    "status": "pending_update",
                    "reason": "chromosome_mismatch",
                    "detail": f"the published catalog covers chromosome '{cchrom}' "
                              f"but the current VCF covers '{vchrom}'; upload a "
                              "matching catalog to enable annotation",
                    "vcf_version": version,
                    "catalog_version": cat,
                }
            gz_path = os.path.join(cfg.data_dir, vrel, GZ_NAME)
            bed_path = os.path.join(cfg.data_dir, crel, bedbuild.GZ_NAME)
            try:
                tb = pysam.TabixFile(gz_path)
                bed = pysam.TabixFile(bed_path)
                try:
                    # 1-based closed [start, end] -> 0-based half-open fetch.
                    lines = tb.fetch(vchrom, start - 1, end)
                    records = []
                    for ln in lines:
                        f = ln.split("\t", 4)
                        pos = int(f[1])
                        ref = f[3]
                        span_end = pos + len(ref) - 1
                        # Exact span check on top of the index iterator.
                        if not (pos <= end and span_end >= start):
                            continue
                        genes = set()
                        # REF span [pos, span_end] 1-based closed -> 0-based
                        # half-open [pos-1, span_end) gene fetch via the index.
                        for gln in bed.fetch(vchrom, pos - 1, span_end):
                            g = gln.split("\t")
                            g_start, g_end, g_name = int(g[1]), int(g[2]), g[3]
                            # Exact overlap check on top of the index iterator.
                            if g_start < span_end and g_end > pos - 1:
                                genes.add(g_name)
                        records.append({
                            "pos": pos,
                            "ref": ref,
                            "span": [pos, span_end],
                            "genes": sorted(genes),
                        })
                finally:
                    tb.close()
                    bed.close()
            except HTTPException:
                raise
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
            return {
                "status": "ok",
                "vcf_version": version,
                "catalog_version": cat,
                "chrom": vchrom,
                "start": start,
                "end": end,
                "count": len(records),
                "records": records,
            }
        finally:
            store.unpin_pair(vref_id, cref_id)

    @app.get("/catalog/version")
    def get_catalog_version():
        info = store.catalog_info()
        if info is None:
            return {"published": False, "current_version": 0}
        return {
            "published": True,
            "current_version": info["catalog"],
            "chrom": info["chrom"],
            "n_genes": info["n_genes"],
            "created_at": info["created_at"],
        }

    @app.get("/version")
    def get_version():
        info = store.version_info()
        if info is None:
            return {"published": False, "current_version": 0}
        return {
            "published": True,
            "current_version": info["version"],
            "chrom": info["chrom"],
            "n_records": info["n_records"],
            "created_at": info["created_at"],
        }

    return app


# Run with: uvicorn --factory app.main:create_app
