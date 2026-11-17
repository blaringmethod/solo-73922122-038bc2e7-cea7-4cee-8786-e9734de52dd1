"""Catalog versioning: pair pinning, reclamation, cross-instance interleavings."""
import os
import threading
import time

import pysam
from conftest import CATALOG, SAMPLE, bed, make_config, upload, upload_catalog, vcf

from app.bedbuild import GZ_NAME as BED_GZ
from app.db import Conflict, Store
from app.main import create_app
from app.vcfbuild import GZ_NAME as VCF_GZ


def _data(app, rel):
    return os.path.join(app.state.cfg.data_dir, rel)


def test_concurrent_catalog_publish_only_one_wins(tmp_path):
    """Two writers racing with the same base_version: exactly one succeeds."""
    store = Store(os.path.join(str(tmp_path), "meta.sqlite3"))
    store.publish_catalog(0, "chr1", "catalogs/c1", 3)

    results = {"ok": 0, "conflict": 0}
    lock = threading.Lock()
    barrier = threading.Barrier(8)

    def race(i):
        barrier.wait()
        try:
            store.publish_catalog(1, "chr1", f"catalogs/racer{i}", 1)
            with lock:
                results["ok"] += 1
        except Conflict:
            with lock:
                results["conflict"] += 1

    threads = [threading.Thread(target=race, args=(i,)) for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert results == {"ok": 1, "conflict": 7}
    assert store.catalog_info()["catalog"] == 2


def test_pin_pair_snapshots_one_moment(tmp_path):
    store = Store(os.path.join(str(tmp_path), "meta.sqlite3"))
    assert store.pin_pair("inst") is None, "no VCF -> nothing to pin"
    store.publish(0, "chr1", "versions/v1", 4)
    pinned = store.pin_pair("inst")
    assert pinned[0] == 1 and pinned[4] is None, "VCF only -> catalog fields are None"
    store.unpin_pair(pinned[3], pinned[7])
    store.publish_catalog(0, "chr1", "catalogs/c1", 2)
    pinned = store.pin_pair("inst")
    assert (pinned[0], pinned[4]) == (1, 1)
    store.unpin_pair(pinned[3], pinned[7])


def test_in_flight_annotation_keeps_vcf_catalog_pair(client, app):
    """A slow query pins (VCF, catalog); updates proceed but must not touch it."""
    upload(client, SAMPLE, base_version=0)
    upload_catalog(client, CATALOG, base_version=0)
    store = app.state.store
    v1_path = store.version_info()["path"]
    c1_path = store.catalog_info()["path"]

    # simulate an in-flight annotation query on the (v1, c1) pair
    pinned = store.pin_pair("test-inst")
    vnum, vchrom, vrel, vref, cnum, cchrom, crel, cref = pinned
    assert (vnum, cnum) == (1, 1)

    # corrections are published while the old query is still open: neither the
    # VCF replacement nor the catalog swap is blocked by the slow query
    r = upload(client, vcf([(500, "A", "G")]), base_version=1)
    assert r.status_code == 201 and r.json()["version"] == 2
    r = upload_catalog(client, bed([(499, 501, "GENE_NEW")]), base_version=1)
    assert r.status_code == 201 and r.json()["catalog_version"] == 2

    # the reaper must not touch either pinned file set
    app.state.reap_once()
    assert os.path.isdir(_data(app, v1_path))
    assert os.path.isdir(_data(app, c1_path))

    # the in-flight query still reads v1 variants and c1 genes via their indexes
    tb = pysam.TabixFile(os.path.join(_data(app, v1_path), VCF_GZ))
    try:
        poss = [ln.split("\t")[1] for ln in tb.fetch("chr1", 99, 106)]
    finally:
        tb.close()
    assert poss == ["100", "105"]
    bed_tb = pysam.TabixFile(os.path.join(_data(app, c1_path), BED_GZ))
    try:
        names = [ln.split("\t")[3] for ln in bed_tb.fetch("chr1", 99, 104)]
    finally:
        bed_tb.close()
    assert names == ["GENE_A"]

    # new queries see the new pair
    body = client.get("/annotations", params={"chrom": "chr1", "start": 1, "end": 1000}).json()
    assert (body["vcf_version"], body["catalog_version"]) == (2, 2)
    assert body["records"][0]["genes"] == ["GENE_NEW"]

    # old query finishes -> reaper removes both old file sets and rows
    store.unpin_pair(vref, cref)
    app.state.reap_once()
    assert not os.path.exists(_data(app, v1_path))
    assert not os.path.exists(_data(app, c1_path))
    assert store.current_version() == 2
    assert store.catalog_info()["catalog"] == 2


def test_in_flight_query_unaffected_by_chromosome_changing_correction(client, app):
    """VCF corrected to another chromosome mid-query: the in-flight query keeps
    the old (VCF, catalog) moment; new queries report pending_update."""
    upload(client, SAMPLE, base_version=0)
    upload_catalog(client, CATALOG, base_version=0)
    store = app.state.store

    pinned = store.pin_pair("test-inst")  # in-flight query on (v1 chr1, c1 chr1)
    vnum, vchrom, vrel, vref, cnum, cchrom, crel, cref = pinned

    # VCF corrected to chr2; the chr1 catalog is still current
    r = upload(client, vcf(["chr2\t50\t.\tA\tG\t.\t.\t."]), base_version=1)
    assert r.status_code == 201

    # new queries must report pending_update, not stale chr1 gene names
    body = client.get("/annotations", params={"chrom": "chr2", "start": 1, "end": 100}).json()
    assert body["status"] == "pending_update"
    assert body["reason"] == "chromosome_mismatch"
    assert "records" not in body
    assert "GENE" not in str(body)

    # the in-flight query still reads the old pair consistently
    tb = pysam.TabixFile(os.path.join(app.state.cfg.data_dir, vrel, VCF_GZ))
    bed_tb = pysam.TabixFile(os.path.join(app.state.cfg.data_dir, crel, BED_GZ))
    try:
        assert [ln.split("\t")[1] for ln in tb.fetch(vchrom, 99, 106)] == ["100", "105"]
        assert [ln.split("\t")[3] for ln in bed_tb.fetch(cchrom, 99, 104)] == ["GENE_A"]
    finally:
        tb.close()
        bed_tb.close()
    store.unpin_pair(vref, cref)


def test_dead_instances_catalog_pins_do_not_block_reclaim(client, app):
    upload(client, SAMPLE, base_version=0)
    upload_catalog(client, CATALOG, base_version=0)
    store = app.state.store
    c1_path = store.catalog_info()["path"]

    now = time.time()
    store.heartbeat("ghost", now - 1000)  # long past dead_after
    with store._txn() as con:
        con.execute(
            "INSERT INTO catalog_refs(catalog, instance, created_at) VALUES (1, 'ghost', ?)",
            (now - 1000,),
        )
    upload_catalog(client, bed([(0, 10, "GENE_NEW")]), base_version=1)
    app.state.reap_once()
    assert not os.path.exists(_data(app, c1_path)), "dead ref must be ignored"


def test_orphan_catalog_dirs_are_swept_and_published_catalog_survives(client, app):
    upload(client, SAMPLE, base_version=0)
    upload_catalog(client, CATALOG, base_version=0)
    cfg = app.state.cfg

    # half-published junk: a catalogs/ dir with data but no DB row
    junk = os.path.join(cfg.data_dir, "catalogs", "junk-half-set")
    os.makedirs(junk)
    open(os.path.join(junk, BED_GZ), "wb").write(b"not really bgzf")
    old = time.time() - 10 * cfg.orphan_grace
    os.utime(junk, (old, old))

    app.state.reap_once()
    assert not os.path.exists(junk)

    # the published catalog is untouched and still annotates
    q = client.get("/annotations", params={"chrom": "chr1", "start": 100, "end": 103})
    assert q.json()["records"][0]["genes"] == ["GENE_A"]


def test_two_instances_share_catalog(tmp_path):
    """Two app instances (as in docker compose) over one shared data dir."""
    from fastapi.testclient import TestClient

    app1 = create_app(make_config(tmp_path, instance_id="inst-1"))
    app2 = create_app(make_config(tmp_path, instance_id="inst-2"))
    with TestClient(app1) as c1, TestClient(app2) as c2:
        assert upload(c1, SAMPLE, base_version=0).status_code == 201
        assert upload_catalog(c1, CATALOG, base_version=0).status_code == 201
        # instance 2 annotates through instance 1's catalog
        q = c2.get("/annotations", params={"chrom": "chr1", "start": 100, "end": 103})
        assert q.json()["records"][0]["genes"] == ["GENE_A"]
        # instance 2 replaces the catalog; instance 1 sees it
        r = upload_catalog(c2, bed([(0, 10, "GENE_NEW")]), base_version=1)
        assert r.status_code == 201 and r.json()["catalog_version"] == 2
        body = c1.get("/annotations", params={"chrom": "chr1", "start": 1, "end": 1000}).json()
        assert body["catalog_version"] == 2
        # a stale replacement attempt from instance 1 fails
        assert upload_catalog(c1, CATALOG, base_version=1).status_code == 409


def test_restart_with_same_data_dir_keeps_catalog_queryable(tmp_path):
    """Simulated process restart: new app object over the same data dir."""
    from fastapi.testclient import TestClient

    app1 = create_app(make_config(tmp_path, instance_id="inst-a"))
    with TestClient(app1) as c1:
        assert upload(c1, SAMPLE, base_version=0).status_code == 201
        assert upload_catalog(c1, CATALOG, base_version=0).status_code == 201

    app2 = create_app(make_config(tmp_path, instance_id="inst-a"))
    with TestClient(app2) as c2:
        body = c2.get("/annotations", params={"chrom": "chr1", "start": 100, "end": 103}).json()
        assert body["catalog_version"] == 1
        assert body["records"][0]["genes"] == ["GENE_A"]
        assert c2.get("/catalog/version").json()["current_version"] == 1
