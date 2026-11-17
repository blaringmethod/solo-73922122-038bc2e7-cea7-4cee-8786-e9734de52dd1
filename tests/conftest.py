"""Shared fixtures/helpers for the test suite."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from app.main import Config, create_app

META = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
CHROM_LINE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
HEADER = META + CHROM_LINE


def vcf(records, header=HEADER):
    """Build a VCF body from (pos, ref, alt) tuples or raw strings."""
    lines = []
    for r in records:
        if isinstance(r, str):
            lines.append(r)
        else:
            pos, ref, alt = r
            lines.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")
    return header + "\n".join(lines) + "\n"


# chr1:100 deletion spanning 100-103, SNP at 105, insertion at 110, del 200-204
SAMPLE = vcf([
    (100, "ACGT", "A"),
    (105, "G", "T"),
    (110, "C", "CT"),
    (200, "TAAAA", "T"),
])


def bed(entries, chrom="chr1"):
    """Build a BED body from (start, end, name) tuples (0-based half-open)."""
    return "".join(f"{chrom}\t{s}\t{e}\t{n}\n" for s, e, n in entries)


# chr1 genes (0-based half-open):
#   GENE_A [89,103)  -> 1-based 90-103, intersects the deletion at 100-103
#   GENE_B [104,115) -> 1-based 105-115, intersects 105 and 110
#   GENE_C [199,204) -> 1-based 200-204, intersects the deletion at 200-204
CATALOG = bed([
    (89, 103, "GENE_A"),
    (104, 115, "GENE_B"),
    (199, 204, "GENE_C"),
])


def make_config(tmp_path, **kw):
    kw.setdefault("instance_id", "test-inst")
    kw.setdefault("reaper_interval", 0)    # deterministic: tests call reap_once()
    kw.setdefault("heartbeat_interval", 0)
    return Config(data_dir=str(tmp_path / "data"), **kw)


@pytest.fixture()
def app(tmp_path):
    return create_app(make_config(tmp_path))


@pytest.fixture()
def client(app):
    from fastapi.testclient import TestClient

    with TestClient(app) as c:
        yield c


def upload(client, body, base_version):
    return client.post(
        "/vcf",
        params={"base_version": base_version},
        files={"file": ("variants.vcf", body if isinstance(body, bytes) else body.encode())},
    )


def upload_catalog(client, body, base_version):
    return client.post(
        "/catalog",
        params={"base_version": base_version},
        files={"file": ("genes.bed", body if isinstance(body, bytes) else body.encode())},
    )
