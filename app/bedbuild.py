"""Gene-catalog BED validation and indexed build (sort -> BGZF -> Tabix).

The accepted format is a strict 4-column BED:

    chrom<TAB>start<TAB>end<TAB>name

with 0-based half-open intervals [start, end), one chromosome per file, at
most 100 entries and unique gene names. Blank lines and '#' comment lines
are ignored. The published file set is genes.bed.gz + genes.bed.gz.tbi, so
gene lookups go through the Tabix index exactly like variant lookups.
"""
from __future__ import annotations

import os
import re

import pysam

CHROM_RE = re.compile(r"\S+")
NAME_RE = re.compile(r"\S+")
# TBI binning supports positions up to 2^29-1 (covers all human chromosomes).
MAX_POS = 2**29 - 1
MAX_ENTRIES = 100

GZ_NAME = "genes.bed.gz"
TBI_NAME = GZ_NAME + ".tbi"


class BEDError(ValueError):
    """The uploaded BED violates the acceptance rules."""


def _fail(lineno: int, msg: str):
    raise BEDError(f"line {lineno}: {msg}")


def _parse(upload_path: str):
    """Parse and validate the upload. Returns (entries, chrom).

    entries is a list of (start, end, name) with 0-based half-open
    coordinates, in file order.
    """
    entries = []
    names = set()
    chrom = None

    try:
        fh = open(upload_path, "r", encoding="utf-8", errors="strict", newline=None)
    except OSError as exc:  # pragma: no cover - defensive
        raise BEDError(f"cannot read upload: {exc}")
    with fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue  # blank lines and comments carry no entries
            fields = line.split("\t")
            if len(fields) != 4:
                _fail(
                    lineno,
                    "expected 4 tab-separated fields (chrom, start, end, name), "
                    f"got {len(fields)}",
                )
            c, start_s, end_s, name = fields
            if not CHROM_RE.fullmatch(c):
                _fail(lineno, "invalid chrom value")
            if chrom is None:
                chrom = c
            elif c != chrom:
                _fail(lineno, f"multiple chromosomes: expected '{chrom}', found '{c}'")
            if not start_s.isdigit():
                _fail(lineno, "start must be a non-negative integer (0-based)")
            if not end_s.isdigit():
                _fail(lineno, "end must be a non-negative integer")
            start = int(start_s)
            end = int(end_s)
            if start >= end:
                _fail(
                    lineno,
                    f"0-based half-open interval requires start < end, "
                    f"got [{start}, {end})",
                )
            if end > MAX_POS:
                _fail(lineno, f"end out of supported range 1..{MAX_POS}")
            if not NAME_RE.fullmatch(name):
                _fail(lineno, "gene name must be non-empty and contain no whitespace")
            if name in names:
                _fail(lineno, f"duplicate gene name '{name}'")
            names.add(name)
            entries.append((start, end, name))
            if len(entries) > MAX_ENTRIES:
                raise BEDError(
                    f"too many entries: at most {MAX_ENTRIES} gene regions allowed"
                )

    if not entries:
        raise BEDError("BED contains no gene entries")
    if chrom is None:  # pragma: no cover - entries imply a chromosome
        raise BEDError("could not determine chromosome")
    return entries, chrom


def validate_and_build(upload_path: str, staging_dir: str):
    """Validate the uploaded BED and build sorted BGZF + Tabix index.

    Writes genes.bed.gz and genes.bed.gz.tbi into staging_dir.
    Returns (chrom, n_entries). Raises BEDError on invalid input.
    """
    try:
        entries, chrom = _parse(upload_path)
    except UnicodeDecodeError as exc:
        raise BEDError(f"upload is not valid UTF-8 text: {exc}")

    # Sort by (start, end); single chromosome, ties keep upload order.
    entries.sort(key=lambda e: (e[0], e[1]))

    sorted_path = os.path.join(staging_dir, "sorted.bed")
    with open(sorted_path, "w", encoding="utf-8", newline="\n") as out:
        for start, end, name in entries:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    gz_path = os.path.join(staging_dir, GZ_NAME)
    pysam.tabix_compress(sorted_path, gz_path, force=True)
    os.remove(sorted_path)
    pysam.tabix_index(gz_path, preset="bed", force=True)

    # Best-effort durability of the finished file set.
    for p in (gz_path, gz_path + ".tbi"):
        try:
            fd = os.open(p, os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
        except OSError:
            pass
    return chrom, len(entries)
