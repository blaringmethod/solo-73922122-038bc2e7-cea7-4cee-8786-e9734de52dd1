# Versioned VCF Variant Store

FastAPI + pysam + SQLite service for uploading single-chromosome VCF files and
querying them by 1-based closed chromosome intervals, with atomic versioned
replacement and safe concurrent access from two app instances. A versioned
gene catalog (BED) annotates the published variants with the genes each
record's reference span overlaps.

## Run

```bash
docker compose up --build
```

Two instances come up sharing one data volume:

- instance 1: http://localhost:8001/docs
- instance 2: http://localhost:8002/docs

Both serve the same published data; uploads may go to either one.

## API

### `POST /vcf?base_version=<n>` (multipart file field `file`)

Accepts a VCF that is:

- at most 10 MiB,
- single-chromosome (all records share one CHROM),
- exactly 8 columns (`#CHROM POS ID REF ALT QUAL FILTER INFO`, no FORMAT/sample
  columns),
- biallelic with `REF`/`ALT` containing only `A/C/G/T` (uppercase; no `N`, no
  symbolic alleles, no `*`, no multi-allelic `A,G`).

Records are sorted by position, BGZF-compressed and Tabix-indexed. The upload
carries the version the client last saw (`base_version`, `0` when the store is
empty). If two replacements race — even across the two instances — exactly one
wins; the loser gets `409` and must re-read `/version` and retry.
Returns `201` with `{"version", "chrom", "n_records"}`.

### `GET /variants?chrom=<c>&start=<s>&end=<e>`

1-based closed interval `[start, end]`. Returns every raw VCF record whose
reference span `[POS, POS+len(REF)-1]` intersects the interval — including
deletions that start before the interval. Fetching goes through the Tabix
index (`TabixFile.fetch`), never a full scan. The response carries the
`version` that was actually read; that version's data and index are pinned for
the duration of the query, so version, data and index always match.

### `POST /catalog?base_version=<n>` (multipart file field `file`)

Publishes a new gene catalog: a 4-column BED (`chrom start end name`) that is

- at most 100 entries (blank lines and `#` comments are ignored),
- single-chromosome,
- 0-based half-open intervals `[start, end)` with `start < end`,
- unique gene names (no whitespace).

Entries are sorted, BGZF-compressed and Tabix-indexed. The catalog has its own
version sequence with the same optimistic-concurrency rule as the VCF
(`base_version` is the catalog version the client last saw, `0` when none).
An invalid BED or a failed write leaves the previously published catalog
untouched. Returns `201` with `{"catalog_version", "chrom", "n_genes"}`.

### `GET /annotations?chrom=<c>&start=<s>&end=<e>`

Annotates the published variants with genes. Variants are fetched through the
VCF Tabix index; for each record, every catalog gene whose 0-based half-open
region intersects the record's reference span `[POS, POS+len(REF)-1]` is
looked up through the catalog's own Tabix index. The query pins the
(VCF version, catalog version) pair in a single transaction, so the whole
response reflects one moment's files, index and catalog — even if a VCF
correction or a catalog swap is published concurrently. The response carries
both `vcf_version` and `catalog_version`.

When no catalog has been published yet, or the published catalog's chromosome
does not match the current VCF's, the response reports
`{"status": "pending_update", "reason": "no_catalog" | "chromosome_mismatch"}`
instead of annotating with stale gene names.

### `GET /version` / `GET /catalog/version`

Current published VCF version / catalog version, chromosome and record/gene
counts.

## Consistency design

- **Atomic publish**: a new version (or catalog) is fully built (sorted,
  compressed, indexed) in a per-instance staging directory, moved into
  `versions/` (or `catalogs/`), and only then the current-version pointer in
  SQLite is flipped inside a single `BEGIN IMMEDIATE` transaction. Readers
  never see a half-built file set.
- **Optimistic concurrency**: the pointer flip checks `base_version` inside
  the same transaction, so concurrent replacements from both instances are
  serialized and only one succeeds.
- **Old queries finish**: a query pins its version in a `refs` table before
  opening files; annotation queries pin the (VCF, catalog) pair in one
  transaction. A background reaper deletes a superseded file set only when no
  pin references it, so slow queries never block updates and updates never
  disturb slow queries. Instances heartbeat; pins of crashed instances expire
  and never block reclamation.
- **Crash recovery**: if the process dies mid-build or mid-publish, the
  staging/orphan directories left behind are swept (with a grace period) and
  the previously published version keeps serving. The pointer flip is the last
  step, so a restart never exposes a partial file set.

## Local development

```bash
pip install -r requirements.txt
python -m pytest tests/          # 81 tests
uvicorn --factory app.main:create_app   # single dev instance (uses /data or $VCF_DATA_DIR)
```

Configuration via environment: `VCF_DATA_DIR`, `VCF_INSTANCE_ID`,
`VCF_MAX_UPLOAD_BYTES`, `VCF_REAPER_INTERVAL`, `VCF_HEARTBEAT_INTERVAL`,
`VCF_DEAD_AFTER`, `VCF_ORPHAN_GRACE`.
