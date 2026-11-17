"""Annotation queries: genes per variant REF span, versions, pending states."""
import json

from conftest import CATALOG, SAMPLE, bed, upload, upload_catalog, vcf


def _annotate(client, chrom="chr1", start=1, end=1000):
    return client.get("/annotations", params={"chrom": chrom, "start": start, "end": end})


def _setup(client):
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    assert upload_catalog(client, CATALOG, base_version=0).status_code == 201


def test_annotate_flow_lists_genes_per_ref_span(client):
    _setup(client)
    r = _annotate(client)
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["status"] == "ok"
    assert body["vcf_version"] == 1
    assert body["catalog_version"] == 1
    assert body["chrom"] == "chr1"
    assert body["count"] == 4
    by_pos = {rec["pos"]: rec for rec in body["records"]}
    # deletion at 100 (REF len 4) spans 100-103 -> GENE_A (1-based 90-103)
    assert by_pos[100]["span"] == [100, 103]
    assert by_pos[100]["genes"] == ["GENE_A"]
    # SNP at 105 and insertion at 110 both fall in GENE_B (1-based 105-115)
    assert by_pos[105]["span"] == [105, 105]
    assert by_pos[105]["genes"] == ["GENE_B"]
    assert by_pos[110]["genes"] == ["GENE_B"]
    # deletion at 200 (REF len 5) spans 200-204 -> GENE_C (1-based 200-204)
    assert by_pos[200]["span"] == [200, 204]
    assert by_pos[200]["genes"] == ["GENE_C"]


def test_zero_based_half_open_convention(client):
    """BED [start, end) is 0-based half-open: it covers 1-based start+1..end."""
    boundary = bed([
        (102, 103, "TAIL"),       # 1-based 103: tail of the 100-103 deletion
        (103, 104, "AFTER_DEL"),  # 1-based 104: just past the deletion
        (104, 105, "AT_SNP"),     # 1-based 105: exactly the SNP
        (105, 106, "AFTER_SNP"),  # 1-based 106: just past the SNP
    ])
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    assert upload_catalog(client, boundary, base_version=0).status_code == 201
    body = _annotate(client).json()
    by_pos = {rec["pos"]: rec["genes"] for rec in body["records"]}
    assert by_pos[100] == ["TAIL"], "deletion span ends at 103; AFTER_DEL must not appear"
    assert by_pos[105] == ["AT_SNP"], "SNP at 105; AFTER_SNP (1-based 106) must not appear"
    assert by_pos[110] == []
    assert by_pos[200] == []


def test_multiple_overlapping_genes_all_listed(client):
    catalog = bed([
        (99, 102, "G1"),   # 1-based 100-102
        (100, 103, "G2"),  # 1-based 101-103
        (90, 110, "G3"),   # 1-based 91-110, also covers 105 and 110
    ])
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    assert upload_catalog(client, catalog, base_version=0).status_code == 201
    body = _annotate(client).json()
    by_pos = {rec["pos"]: rec["genes"] for rec in body["records"]}
    assert by_pos[100] == ["G1", "G2", "G3"]
    assert by_pos[105] == ["G3"]
    assert by_pos[110] == ["G3"]
    assert by_pos[200] == []


def test_annotation_is_scoped_to_query_interval(client):
    _setup(client)
    body = _annotate(client, start=105, end=110).json()
    assert [rec["pos"] for rec in body["records"]] == [105, 110]
    body = _annotate(client, start=106, end=109).json()
    assert body["records"] == []


def test_pending_update_when_no_catalog(client):
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    r = _annotate(client)
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["status"] == "pending_update"
    assert body["reason"] == "no_catalog"
    assert body["vcf_version"] == 1
    assert body["catalog_version"] is None
    assert "records" not in body


def test_pending_update_on_chromosome_mismatch_never_returns_old_genes(client):
    """VCF corrected to another chromosome: old catalog genes must not leak."""
    _setup(client)  # VCF v1 chr1, catalog c1 chr1
    assert _annotate(client).json()["status"] == "ok"

    # VCF corrected: now chr2, catalog still describes chr1.
    r = upload(client, vcf(["chr2\t50\t.\tA\tG\t.\t.\t."]), base_version=1)
    assert r.status_code == 201 and r.json()["version"] == 2

    body = _annotate(client, chrom="chr2", start=1, end=100).json()
    assert body["status"] == "pending_update"
    assert body["reason"] == "chromosome_mismatch"
    assert body["vcf_version"] == 2
    assert body["catalog_version"] == 1
    assert "records" not in body
    assert "GENE" not in json.dumps(body), "stale gene names must not be reported"

    # A matching catalog flips annotation back to ok.
    r = upload_catalog(client, bed([(40, 60, "GENE_Z")], chrom="chr2"), base_version=1)
    assert r.status_code == 201 and r.json()["catalog_version"] == 2
    body = _annotate(client, chrom="chr2", start=1, end=100).json()
    assert body["status"] == "ok"
    assert body["vcf_version"] == 2
    assert body["catalog_version"] == 2
    assert body["records"][0]["genes"] == ["GENE_Z"]


def test_pending_update_when_catalog_on_other_chromosome(client):
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    r = upload_catalog(client, bed([(0, 50, "GENE_X")], chrom="chr2"), base_version=0)
    assert r.status_code == 201
    body = _annotate(client).json()
    assert body["status"] == "pending_update"
    assert body["reason"] == "chromosome_mismatch"
    assert "GENE_X" not in str(body)


def test_annotation_without_vcf_is_404(client):
    assert _annotate(client).status_code == 404


def test_annotation_query_validation(client):
    _setup(client)
    assert _annotate(client, start=0).status_code == 422
    assert _annotate(client, start=10, end=5).status_code == 400
    r = _annotate(client, chrom="chr2")
    assert r.status_code == 400 and "chr1" in r.json()["detail"]


def test_response_tracks_both_version_sequences(client):
    _setup(client)
    body = _annotate(client).json()
    assert (body["vcf_version"], body["catalog_version"]) == (1, 1)
    # replace only the VCF
    assert upload(client, vcf([(500, "A", "G")]), base_version=1).status_code == 201
    body = _annotate(client).json()
    assert (body["vcf_version"], body["catalog_version"]) == (2, 1)
    # replace only the catalog
    assert upload_catalog(client, bed([(499, 501, "GENE_Q")]), base_version=1).status_code == 201
    body = _annotate(client).json()
    assert (body["vcf_version"], body["catalog_version"]) == (2, 2)
    assert body["records"][0]["genes"] == ["GENE_Q"]


def test_catalog_version_endpoint(client):
    assert client.get("/catalog/version").json() == {"published": False, "current_version": 0}
    _setup(client)
    body = client.get("/catalog/version").json()
    assert body["published"] is True
    assert body["current_version"] == 1
    assert body["chrom"] == "chr1"
    assert body["n_genes"] == 3
