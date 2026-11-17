"""BED catalog upload validation rules."""
import os

import pytest
from conftest import CATALOG, SAMPLE, bed, upload, upload_catalog, vcf

from app.bedbuild import GZ_NAME, MAX_ENTRIES, TBI_NAME


def test_valid_bed_upload(client, app):
    r = upload_catalog(client, CATALOG, base_version=0)
    assert r.status_code == 201, r.text
    assert r.json() == {"catalog_version": 1, "chrom": "chr1", "n_genes": 3}
    info = app.state.store.catalog_info()
    d = os.path.join(app.state.cfg.data_dir, info["path"])
    assert os.path.isfile(os.path.join(d, GZ_NAME))
    assert os.path.isfile(os.path.join(d, TBI_NAME))


def test_exactly_max_entries_accepted(client):
    body = bed([(i * 10, i * 10 + 5, f"G{i}") for i in range(MAX_ENTRIES)])
    r = upload_catalog(client, body, base_version=0)
    assert r.status_code == 201, r.text
    assert r.json()["n_genes"] == MAX_ENTRIES


def test_too_many_entries_rejected(client, app):
    body = bed([(i * 10, i * 10 + 5, f"G{i}") for i in range(MAX_ENTRIES + 1)])
    r = upload_catalog(client, body, base_version=0)
    assert r.status_code == 400 and "invalid BED" in r.json()["detail"]
    assert app.state.store.catalog_info() is None


BAD_BEDS = {
    "two_chromosomes": "chr1\t0\t10\tG1\nchr2\t0\t10\tG2\n",
    "start_equals_end": "chr1\t10\t10\tG1\n",
    "start_after_end": "chr1\t20\t10\tG1\n",
    "negative_start": "chr1\t-1\t10\tG1\n",
    "start_not_a_number": "chr1\tx\t10\tG1\n",
    "end_not_a_number": "chr1\t0\ty\tG1\n",
    "duplicate_names": "chr1\t0\t10\tG1\nchr1\t20\t30\tG1\n",
    "three_columns": "chr1\t0\t10\n",
    "five_columns": "chr1\t0\t10\tG1\textra\n",
    "name_with_space": "chr1\t0\t10\tGE NE\n",
    "empty_name": "chr1\t0\t10\t\n",
    "empty_file": "",
    "comments_only": "# just a comment\n#track nothing\n",
    "end_beyond_tabix_limit": "chr1\t0\t536870912\tG1\n",
}


@pytest.mark.parametrize("name", sorted(BAD_BEDS))
def test_invalid_beds_are_rejected(client, app, name):
    r = upload_catalog(client, BAD_BEDS[name], base_version=0)
    assert r.status_code == 400, f"{name}: {r.status_code} {r.text}"
    assert "invalid BED" in r.json()["detail"]
    # nothing published, nothing left behind
    assert app.state.store.catalog_info() is None
    assert os.listdir(os.path.join(app.state.cfg.data_dir, "catalogs")) == []


def test_zero_based_start_allowed(client):
    r = upload_catalog(client, bed([(0, 50, "G1")]), base_version=0)
    assert r.status_code == 201, r.text


def test_comments_and_blank_lines_are_ignored(client):
    body = "# genome browser comment\n\nchr1\t0\t10\tG1\n# another comment\nchr1\t20\t30\tG2\n"
    r = upload_catalog(client, body, base_version=0)
    assert r.status_code == 201, r.text
    assert r.json()["n_genes"] == 2


def test_unsorted_bed_is_sorted_and_queryable(client):
    upload(client, SAMPLE, base_version=0)
    body = bed([(199, 204, "GENE_C"), (89, 103, "GENE_A"), (104, 115, "GENE_B")])
    assert upload_catalog(client, body, base_version=0).status_code == 201
    q = client.get("/annotations", params={"chrom": "chr1", "start": 1, "end": 1000})
    by_pos = {rec["pos"]: rec["genes"] for rec in q.json()["records"]}
    assert by_pos[100] == ["GENE_A"]
    assert by_pos[200] == ["GENE_C"]


def test_binary_upload_rejected_as_invalid_bed(client):
    import gzip

    r = upload_catalog(client, gzip.compress(CATALOG.encode()), base_version=0)
    assert r.status_code == 400
    assert "invalid BED" in r.json()["detail"]


def test_first_catalog_upload_requires_base_version_zero(client):
    r = upload_catalog(client, CATALOG, base_version=3)
    assert r.status_code == 409
    assert "current catalog version 0" in r.json()["detail"]


def test_replace_catalog_requires_fresh_base_version(client):
    assert upload_catalog(client, CATALOG, base_version=0).status_code == 201
    r = upload_catalog(client, CATALOG, base_version=0)
    assert r.status_code == 409
    assert "current catalog version 1" in r.json()["detail"]
    r = upload_catalog(client, CATALOG, base_version=1)
    assert r.status_code == 201 and r.json()["catalog_version"] == 2


def test_invalid_bed_keeps_original_catalog(client):
    """A failed replacement must leave the published catalog fully usable."""
    upload(client, SAMPLE, base_version=0)
    assert upload_catalog(client, CATALOG, base_version=0).status_code == 201
    r = upload_catalog(client, BAD_BEDS["duplicate_names"], base_version=1)
    assert r.status_code == 400
    assert client.get("/catalog/version").json()["current_version"] == 1
    q = client.get("/annotations", params={"chrom": "chr1", "start": 100, "end": 103})
    assert q.json()["status"] == "ok"
    assert q.json()["records"][0]["genes"] == ["GENE_A"]


def test_failed_upload_leaves_no_staging(client, app):
    upload_catalog(client, BAD_BEDS["two_chromosomes"], base_version=0)
    staging = os.path.join(app.state.cfg.data_dir, "staging", app.state.cfg.instance_id)
    assert os.listdir(staging) == []
