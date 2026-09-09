# Protein-Content Accumulation Curve Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the accumulation-curve analysis that quantifies how much new BFD protein content (internal, vs. BFD's own accumulated set; external, vs. public references) is revealed per genome added, with per-clade attribution, a pan-proteome power-law fit, and the seven required validation checks, all computed by resampling a genome × cluster incidence matrix built once from an existing clustering run.

**Architecture:** A small pixi-managed Python+R toolchain of focused scripts, each independently unit-testable against synthetic fixtures shaped like the real pilot data: a version-tag/provenance layer (Tasks 2-4), an incidence-matrix builder that joins clustering output + taxonomy + QC metrics from `BFD.duckdb` (Task 5), a pure-function resampling/attribution/fit library (Tasks 6-8) wrapped by a thin CLI (Task 9), two validation scripts consuming the CLI's raw per-permutation output (Tasks 10-11), plotting scripts (Tasks 12-13), and finally two new Nextflow processes wiring it all into `protein_novelty.nf` (Task 14).

**Tech Stack:** Python 3.11+ (numpy, pyarrow, duckdb — no scipy: the power-law fit uses log-log OLS via numpy, not `scipy.optimize.curve_fit`, to keep the dependency footprint minimal, consistent with this repo's existing stdlib-only `bin/` scripts), R (ggplot2, arrow, jsonlite) for publication plots, Plotly (Python) for the three interactive figures. Managed via a new `pixi.toml` (`[workspace]` table per project convention).

**Spec:** `docs/superpowers/specs/2026-09-09-protein-novelty-accumulation-curve-design.md`

## Global Constraints

- Seed the permutation/genome universe from ASMIDs present in `bfd_proteins_provenance.tsv`, **never** from `samples.csv` directly (repeats the L-8 false-singleton bug at this layer otherwise).
- Duplicate ASMIDs are a data-integrity error: assert uniqueness, fail loudly, never silently dedupe.
- A genome missing its taxonomy value at the density-adaptive rank rolls up to the next coarser rank; if unresolvable to phylum, bucket as `unclassified`.
- Version-tag hash = first 7 hex chars of SHA-256 over the manifest **with the `version_tag` field itself excluded** (hash-then-embed, never circular).
- "External novelty" always means no-hit against **both** NCBI nr (`nr_cluster_seq.dmnd`, datestamped `20260128`) and the in-house UniProt-fungi set — record both identities in every manifest.
- BUSCO/N50 come from `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/db/BFD.duckdb` (`busco_genome.parquet.complete_pct`, `asm_stats.parquet.N50_bp`, both ASMID-keyed) — not from `samples.csv`. Verify ASMID coverage against the genome set in use each run; this DuckDB (built 2026-08-05) will need rebuilding before scale-up.
- CI band convention everywhere in this design: **95% CI = 2.5/97.5 percentiles** of the per-permutation distribution. Default `n_permutations = 1000`, explicitly provisional pending the convergence check (Task 9's `--n-permutations` sweep).
- Pan-proteome power-law fit is `n_new(N) ≈ κ·N^(−α)` fit to the **per-permutation marginal (rate) series**, never to the cumulative curve (cumulative R² is uninformative — see spec). Decision rule: `α ≤ 1` open, `α > 1` closed/closing.
- Large tabular outputs (incidence matrix, raw per-permutation marginals) are Parquet + zstd; the three small summary TSVs stay plain uncompressed TSV.
- Plotting: ggplot2 (R) is the primary implementation covering all 10 figures; Plotly (Python) covers only figures 1-3.
- Every output file/figure caption carries the run's version tag.

---

## File Structure

New files, all in the `DeltaGain` repo:

- `pixi.toml` — workspace environment (Task 1)
- `bin/version_tag.py` — short hash-based version-tag builder, shared by UniProt and BFD manifests (Task 2)
- `bin/parse_uniprot_relnotes.py` — parses UniProt's `relnotes.txt` for release version/date (Task 3)
- `bin/build_uniprot_manifest.py` — writes `uniprot_fungi_nr.manifest.json` (Task 3)
- `modules/refdb/FETCH_UNIPROT_FUNGI_FTP/main.nf` — **modify**: also fetch `relnotes.txt` (Task 3)
- `bin/build_cluster_manifest.py` — writes `bfd_cluster_manifest.json` (Task 4)
- `bin/build_cluster_incidence_matrix.py` — CLI: parses cluster/provenance/novelty/taxonomy files + joins `BFD.duckdb`, applies genome-universe edge-case rules, writes `incidence_matrix.parquet`, `cluster_external_status.parquet`, `genome_metadata.parquet` (Task 5)
- `bin/lib_accumulation.py` — pure-function library: permutation walk, clade Shapley attribution, power-law fit (Tasks 6-8)
- `bin/accumulation_curve.py` — CLI wrapping `lib_accumulation.py`: writes `accumulation_curve_summary.tsv`, `pangenome_powerlaw_fit.tsv`, `clade_contribution.tsv`, `permutation_marginals.parquet` (Task 9)
- `bin/validate_positive_control.py` — validation #4, conditioned check (Task 10)
- `bin/check_qc_confounds.py` — validation #3, confound correlations (Task 11)
- `bin/plot_accumulation_ggplot.R` — figures 1-10, static PDF (Task 12)
- `bin/plot_accumulation_plotly.py` — figures 1-3, interactive HTML (Task 13)
- `modules/binning/BUILD_CLUSTER_INCIDENCE_MATRIX/main.nf`, `modules/binning/ACCUMULATION_CURVE/main.nf` — new Nextflow processes; `protein_novelty.nf` — **modify**, wire them in after `MMSEQS_LINCLUST_BFD` + `BIN_NOVELTY_HITS` (Task 14)

Test fixtures, mirroring real column formats confirmed from the actual pilot outputs (`DeltaGain_Fungi/results/protein_novelty/`):

- `tests/data/accumulation/bfd_proteins_clusters.tsv` — 2-column, no header, MMseqs2 `createtsv` format (`repr_id \t member_id`)
- `tests/data/accumulation/bfd_proteins_provenance.tsv` — header `new_id\tasmid\toriginal_id`
- `tests/data/accumulation/bfd_novelty_bins.tsv` — header `query_id\t...\toverall_status`
- `tests/data/accumulation/samples.csv` — same columns as the real file (`ASMID,...,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,...`)
- `tests/data/accumulation/genome_classification.tsv` — header `asmid\tcomponent_id\tcomponent_size\t...\tcluster_class`

---

## Task 1: Pixi workspace environment

**Files:**
- Create: `pixi.toml`
- Test: manual environment activation (no test framework needed for this task)

**Interfaces:**
- Produces: a `pixi run python`, `pixi run Rscript` environment with `numpy`, `pyarrow`, `duckdb` (Python) and `r-base`, `r-ggplot2`, `r-arrow`, `r-jsonlite` (R) available, for every later task.

- [ ] **Step 1: Write `pixi.toml`**

```toml
[workspace]
name = "deltagain-accumulation-curve"
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64"]

[dependencies]
python = ">=3.11"
numpy = ">=1.26"
pyarrow = ">=15"
duckdb = ">=1.0"
matplotlib = ">=3.8"
plotly = ">=5.20"
r-base = ">=4.3"
r-ggplot2 = ">=3.5"
r-arrow = ">=15"
r-jsonlite = ">=1.8"
pytest = ">=8"

[tasks]
test = "pytest tests/ -v"
```

- [ ] **Step 2: Verify the environment resolves and activates**

Run: `pixi install && pixi run python -c "import numpy, pyarrow, duckdb; print('ok')"`
Expected: `ok` printed, no resolution errors.

Run: `pixi run Rscript -e 'library(ggplot2); library(arrow); library(jsonlite); cat("ok\n")'`
Expected: `ok` printed.

- [ ] **Step 3: Commit**

```bash
git add pixi.toml
git commit -m "Add pixi workspace environment for accumulation-curve analysis

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 2: Version-tag utility

**Files:**
- Create: `bin/version_tag.py`
- Test: `tests/test_version_tag.py`

**Interfaces:**
- Produces: `build_version_tag(fields: dict, prefix: str) -> str` — `fields` is an ordered dict of human-readable identifying values (e.g. `{"date": "20260909", "n": "11024", "params": "mmseqs95c90"}`); returns e.g. `"bfd-v20260909-n11024-mmseqs95c90-a1b2c3d"`. Also `manifest_hash(manifest: dict) -> str` — SHA-256 of the JSON-serialized manifest **with any `version_tag` key removed first**, truncated to 7 hex chars.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_version_tag.py
import json
from bin.version_tag import build_version_tag, manifest_hash

def test_manifest_hash_excludes_version_tag_field():
    manifest_without_tag = {"date": "20260909", "n": 100}
    manifest_with_tag = dict(manifest_without_tag, version_tag="whatever")
    assert manifest_hash(manifest_without_tag) == manifest_hash(manifest_with_tag)

def test_manifest_hash_is_7_hex_chars():
    h = manifest_hash({"a": 1})
    assert len(h) == 7
    assert all(c in "0123456789abcdef" for c in h)

def test_manifest_hash_changes_with_content():
    assert manifest_hash({"a": 1}) != manifest_hash({"a": 2})

def test_build_version_tag_format():
    fields = {"date": "20260909", "n": "11024", "params": "mmseqs95c90"}
    tag = build_version_tag(fields, prefix="bfd", manifest={"date": "20260909", "n": 11024})
    parts = tag.split("-")
    assert parts[0] == "bfd"
    assert parts[1] == "v20260909"
    assert parts[2] == "n11024"
    assert parts[3] == "mmseqs95c90"
    assert len(parts[4]) == 7
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_version_tag.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'bin.version_tag'`

- [ ] **Step 3: Implement `bin/version_tag.py`**

```python
#!/usr/bin/env python3
"""Short hash-based version-tag builder shared by UniProt reference and BFD
dataset provenance manifests (design doc: Data provenance & versioning).

The hash is always computed over the manifest with any existing
`version_tag` key removed first -- hashing a manifest that already embeds
its own tag would make the tag circular (its value would depend on itself).
"""
import hashlib
import json


def manifest_hash(manifest: dict) -> str:
    clean = {k: v for k, v in manifest.items() if k != "version_tag"}
    serialized = json.dumps(clean, sort_keys=True, default=str)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:7]


def build_version_tag(fields: dict, prefix: str, manifest: dict) -> str:
    """fields: ordered {label: value} pairs rendered as `label+value`
    segments (e.g. {"date": "20260909"} -> "v20260909"). The first field is
    rendered with a leading 'v', the rest bare -- matching the spec's
    example `bfd-v20260909-n11024-mmseqs95c90-a1b2c3d`."""
    segments = [prefix]
    for i, (label, value) in enumerate(fields.items()):
        segments.append(f"v{value}" if i == 0 else f"{label}{value}")
    segments.append(manifest_hash(manifest))
    return "-".join(segments)
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/test_version_tag.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add bin/version_tag.py tests/test_version_tag.py
git commit -m "Add shared version-tag builder for provenance manifests

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 3: UniProt release provenance

**Files:**
- Create: `bin/parse_uniprot_relnotes.py`
- Create: `bin/build_uniprot_manifest.py`
- Modify: `modules/refdb/FETCH_UNIPROT_FUNGI_FTP/main.nf`
- Test: `tests/test_parse_uniprot_relnotes.py`, `tests/test_build_uniprot_manifest.py`
- Test fixture: `tests/data/accumulation/relnotes_excerpt.txt`

**Interfaces:**
- Consumes: `manifest_hash`, `build_version_tag` from Task 2.
- Produces: `parse_relnotes(text: str) -> dict` returning `{"release": "2026_03", "release_date": "2026-09-02"}`; `build_uniprot_manifest(release_info: dict, mmseqs_version: str, mmseqs_params: dict, seq_count_in: int, seq_count_out: int) -> dict` returning the full manifest dict (including `version_tag`), used by later Nextflow wiring (Task 14).

- [ ] **Step 1: Write the fixture**

```
# tests/data/accumulation/relnotes_excerpt.txt
UniProt Release 2026_03
==========================

Header of the UniProt Knowledgebase Release 2026_03 (02-Sept-2026)
========================================================================

This release contains the following:
...
```

- [ ] **Step 2: Write the failing tests**

```python
# tests/test_parse_uniprot_relnotes.py
from bin.parse_uniprot_relnotes import parse_relnotes

def test_parses_real_relnotes_format():
    text = open("tests/data/accumulation/relnotes_excerpt.txt").read()
    info = parse_relnotes(text)
    assert info["release"] == "2026_03"
    assert info["release_date"] == "2026-09-02"

def test_raises_on_unparseable_text():
    import pytest
    with pytest.raises(ValueError):
        parse_relnotes("nothing recognizable here")
```

```python
# tests/test_build_uniprot_manifest.py
from bin.build_uniprot_manifest import build_uniprot_manifest

def test_manifest_has_required_fields_and_tag():
    release_info = {"release": "2026_03", "release_date": "2026-09-02"}
    manifest = build_uniprot_manifest(
        release_info=release_info,
        mmseqs_version="15.6f452",
        mmseqs_params={"min_seq_id": 0.95, "c": 0.9, "cov_mode": 0},
        seq_count_in=17_500_000,
        seq_count_out=15_089_757,
    )
    assert manifest["release"] == "2026_03"
    assert manifest["seq_count_out"] == 15_089_757
    assert "version_tag" in manifest
    assert manifest["version_tag"].startswith("uniprot-fungi-v2026_03-")
```

- [ ] **Step 3: Run to verify failures**

Run: `pixi run pytest tests/test_parse_uniprot_relnotes.py tests/test_build_uniprot_manifest.py -v`
Expected: FAIL, both modules missing.

- [ ] **Step 4: Implement `bin/parse_uniprot_relnotes.py`**

```python
#!/usr/bin/env python3
"""Parse UniProt's relnotes.txt for the current release version and date.

Real format confirmed 2026-09-09 against
https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt:
    Header of the UniProt Knowledgebase Release 2026_03 (02-Sept-2026)
"""
import re
from datetime import datetime

_RELEASE_RE = re.compile(
    r"UniProt Knowledgebase Release (\d{4}_\d{2}) \((\d{1,2}-[A-Za-z]+-\d{4})\)"
)
_MONTHS = {m: i for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
     "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"], start=1)}


def parse_relnotes(text: str) -> dict:
    match = _RELEASE_RE.search(text)
    if not match:
        raise ValueError("could not find a release line in relnotes.txt")
    release, date_str = match.groups()
    day, month_str, year = date_str.split("-")
    release_date = datetime(int(year), _MONTHS[month_str], int(day)).strftime("%Y-%m-%d")
    return {"release": release, "release_date": release_date}
```

- [ ] **Step 5: Implement `bin/build_uniprot_manifest.py`**

```python
#!/usr/bin/env python3
"""Build the uniprot_fungi_nr.manifest.json provenance sidecar (design doc:
UniProt reference versioning)."""
from version_tag import build_version_tag


def build_uniprot_manifest(release_info: dict, mmseqs_version: str,
                            mmseqs_params: dict, seq_count_in: int,
                            seq_count_out: int) -> dict:
    manifest = {
        "release": release_info["release"],
        "release_date": release_info["release_date"],
        "mmseqs_version": mmseqs_version,
        "mmseqs_params": mmseqs_params,
        "seq_count_in": seq_count_in,
        "seq_count_out": seq_count_out,
    }
    manifest["version_tag"] = build_version_tag(
        {"release": release_info["release"], "n": seq_count_out},
        prefix="uniprot-fungi",
        manifest=manifest,
    )
    return manifest


if __name__ == "__main__":
    import argparse
    import json
    from parse_uniprot_relnotes import parse_relnotes

    ap = argparse.ArgumentParser()
    ap.add_argument("--relnotes", required=True)
    ap.add_argument("--mmseqs-version", required=True)
    ap.add_argument("--min-seq-id", type=float, required=True)
    ap.add_argument("--cov", type=float, required=True)
    ap.add_argument("--cov-mode", type=int, required=True)
    ap.add_argument("--seq-count-in", type=int, required=True)
    ap.add_argument("--seq-count-out", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    release_info = parse_relnotes(open(args.relnotes).read())
    manifest = build_uniprot_manifest(
        release_info, args.mmseqs_version,
        {"min_seq_id": args.min_seq_id, "c": args.cov, "cov_mode": args.cov_mode},
        args.seq_count_in, args.seq_count_out,
    )
    json.dump(manifest, open(args.out, "w"), indent=2)
```

Note: the `build_version_tag` example format in the spec uses `v<date>` as the
first segment; here the first identifying field is the release string, so the
tag renders as `uniprot-fungi-v2026_03-n15089757-<hash>` — consistent with
`build_version_tag`'s "first field gets the `v` prefix" convention from Task 2.

- [ ] **Step 6: Run to verify tests pass**

Run: `pixi run pytest tests/test_parse_uniprot_relnotes.py tests/test_build_uniprot_manifest.py -v`
Expected: PASS (3 tests)

- [ ] **Step 7: Modify `FETCH_UNIPROT_FUNGI_FTP` to also fetch `relnotes.txt`**

Read the current module first (`modules/refdb/FETCH_UNIPROT_FUNGI_FTP/main.nf`) to preserve its existing `-C -` resume behavior and stub block exactly. Add a second `curl` for relnotes.txt and emit it as a new output channel:

```groovy
// (inside the existing process, after the existing curl for ${remote_name}):
    curl -fsSL --retry 5 --retry-delay 30 \\
        "https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt" \\
        -o relnotes.txt
```

Add `path("relnotes.txt"), emit: relnotes` to the `output:` block (only needs
to be emitted once, not once per division — if this process is called per
division via a channel of tuples as `protein_novelty.nf` currently does,
move the relnotes fetch to run once, e.g. gated on `division == 'sprot'`, to
avoid a redundant duplicate emit), and extend the `stub:` block with
`printf "Header of the UniProt Knowledgebase Release 2026_03 (02-Sept-2026)\n" > relnotes.txt`.

- [ ] **Step 8: Commit**

```bash
git add bin/parse_uniprot_relnotes.py bin/build_uniprot_manifest.py \
    tests/test_parse_uniprot_relnotes.py tests/test_build_uniprot_manifest.py \
    tests/data/accumulation/relnotes_excerpt.txt \
    modules/refdb/FETCH_UNIPROT_FUNGI_FTP/main.nf
git commit -m "Capture UniProt relnotes.txt release version for provenance

Closes the gap where FETCH_UNIPROT_FUNGI_FTP hit current_release/ blind,
with no release version or date recorded anywhere in its output.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 4: BFD dataset manifest

**Files:**
- Create: `bin/build_cluster_manifest.py`
- Test: `tests/test_build_cluster_manifest.py`

**Interfaces:**
- Consumes: `build_version_tag` from Task 2.
- Produces: `build_cluster_manifest(run_date, n_genomes_total, n_genomes_annotated, mmseqs_version, mmseqs_params, qc_thresholds, screened_out, zero_protein_asmids) -> dict`, used by Task 14's Nextflow wiring.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_build_cluster_manifest.py
from bin.build_cluster_manifest import build_cluster_manifest

def test_manifest_records_screen_out_reasons_separately():
    manifest = build_cluster_manifest(
        run_date="2026-09-09",
        n_genomes_total=23683,
        n_genomes_annotated=10922,
        mmseqs_version="15.6f452",
        mmseqs_params={"min_seq_id": 0.95, "c": 0.9, "cov_mode": 0},
        qc_thresholds={"busco_complete_pct_min": 70.0, "n50_bp_min": 10000},
        screened_out={"busco_complete_pct_min": ["ASM1", "ASM2"], "n50_bp_min": ["ASM3"]},
        zero_protein_asmids=["ASM4"],
    )
    assert manifest["n_genomes_annotated"] == 10922
    assert manifest["screened_out"]["busco_complete_pct_min"] == ["ASM1", "ASM2"]
    assert manifest["zero_protein_asmids"] == ["ASM4"]
    assert "version_tag" in manifest
    assert manifest["version_tag"].startswith("bfd-v20260909-")
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_build_cluster_manifest.py -v`
Expected: FAIL, module missing.

- [ ] **Step 3: Implement `bin/build_cluster_manifest.py`**

```python
#!/usr/bin/env python3
"""Build the bfd_cluster_manifest.json provenance sidecar for one BFD
clustering run (design doc: BFD dataset versioning). Records genome count,
QC screen-outs (by reason, not just the final kept set), and clustering
params, so multiple clustering iterations with different parameter/QC
choices stay distinguishable."""
from version_tag import build_version_tag


def build_cluster_manifest(run_date: str, n_genomes_total: int,
                            n_genomes_annotated: int, mmseqs_version: str,
                            mmseqs_params: dict, qc_thresholds: dict,
                            screened_out: dict, zero_protein_asmids: list) -> dict:
    date_compact = run_date.replace("-", "")
    manifest = {
        "run_date": run_date,
        "n_genomes_total": n_genomes_total,
        "n_genomes_annotated": n_genomes_annotated,
        "mmseqs_version": mmseqs_version,
        "mmseqs_params": mmseqs_params,
        "qc_thresholds": qc_thresholds,
        "screened_out": screened_out,
        "zero_protein_asmids": zero_protein_asmids,
    }
    params_short = f"mmseqs{int(mmseqs_params['min_seq_id']*100)}c{int(mmseqs_params['c']*100)}"
    manifest["version_tag"] = build_version_tag(
        {"date": date_compact, "n": n_genomes_annotated, "params": params_short},
        prefix="bfd",
        manifest=manifest,
    )
    return manifest
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/test_build_cluster_manifest.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add bin/build_cluster_manifest.py tests/test_build_cluster_manifest.py
git commit -m "Add BFD clustering-run provenance manifest builder

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 5: Cluster incidence matrix builder

**Files:**
- Create: `bin/build_cluster_incidence_matrix.py`
- Test: `tests/test_build_cluster_incidence_matrix.py`
- Test fixtures: `tests/data/accumulation/bfd_proteins_clusters.tsv`, `bfd_proteins_provenance.tsv`, `bfd_novelty_bins.tsv`, `samples.csv`, `busco_genome.parquet`, `asm_stats.parquet`, `genome_classification.tsv` — all in the exact real column formats confirmed against the actual pilot outputs

**Interfaces:**
- Produces: `build_incidence_matrix(clusters_tsv, provenance_tsv, novelty_bins_tsv, samples_csv, busco_parquet, asm_stats_parquet, genome_classification_tsv, min_clade_n) -> (incidence_df, external_status_df, genome_metadata_df)` — three pandas/pyarrow-backed DataFrames, written to Parquet by the CLI. `genome_metadata_df` columns: `asmid, n_proteins, complete_pct, n50_bp, component_id, cluster_class, clade_rank, clade_label, phylum..species`.

- [ ] **Step 1: Write the fixtures** (small, hand-constructed, matching real formats exactly)

```
# tests/data/accumulation/bfd_proteins_clusters.tsv (no header, repr<TAB>member)
GENOME_A__p1	GENOME_A__p1
GENOME_A__p1	GENOME_B__p1
GENOME_A__p2	GENOME_A__p2
GENOME_B__p2	GENOME_B__p2
GENOME_C__p1	GENOME_C__p1
```

```
# tests/data/accumulation/bfd_proteins_provenance.tsv
new_id	asmid	original_id
GENOME_A__p1	GENOME_A	p1
GENOME_A__p2	GENOME_A	p2
GENOME_B__p1	GENOME_B	p1
GENOME_B__p2	GENOME_B	p2
GENOME_C__p1	GENOME_C	p1
```

```
# tests/data/accumulation/bfd_novelty_bins.tsv (only overall_status matters here)
query_id	overall_status
GENOME_A__p1	strong_hit
GENOME_A__p2	no_hit
GENOME_B__p2	weak_hit
GENOME_C__p1	no_hit
```

```
# tests/data/accumulation/samples.csv (header matches real samples.csv exactly)
ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,TRANSL_TABLE,LOCUSTAG
GENOME_A,Foo bar,S1,PRJ1,1,dikarya,Ascomycota,,Dothideomycetes,,Pleosporales,,Foo,Foo bar,1,X1
GENOME_B,Foo baz,S2,PRJ1,2,dikarya,Ascomycota,,Dothideomycetes,,Pleosporales,,Foo,Foo baz,1,X2
GENOME_C,Qux quux,S3,PRJ1,3,dikarya,Basidiomycota,,Agaricomycetes,,Agaricales,,Qux,Qux quux,1,X3
```

```
# tests/data/accumulation/genome_classification.tsv
asmid	component_id	component_size	existing_species	majority_species_in_component	best_within_component_ani	best_ani_partner	cluster_class
GENOME_A	comp_1	2	Foo bar	Foo bar	99.5	GENOME_B	redundant_strain
GENOME_B	comp_1	2	Foo baz	Foo bar	99.5	GENOME_A	redundant_strain
GENOME_C	comp_2	1	Qux quux	Qux quux			singleton_isolated
```

Build `busco_genome.parquet` and `asm_stats.parquet` fixtures with a tiny Python snippet run once in the test itself (pyarrow `Table.from_pydict`), rather than committing binary files:

```python
# inside tests/test_build_cluster_incidence_matrix.py, a fixture helper
import pyarrow as pa
import pyarrow.parquet as pq

def _write_qc_fixtures(tmp_path):
    busco = pa.table({"ASMID": ["GENOME_A", "GENOME_B", "GENOME_C"],
                       "complete_pct": [98.5, 97.0, 60.0]})
    asm = pa.table({"ASMID": ["GENOME_A", "GENOME_B", "GENOME_C"],
                     "N50_bp": [500000, 480000, 20000]})
    pq.write_table(busco, tmp_path / "busco_genome.parquet")
    pq.write_table(asm, tmp_path / "asm_stats.parquet")
    return tmp_path / "busco_genome.parquet", tmp_path / "asm_stats.parquet"
```

- [ ] **Step 2: Write the failing tests**

```python
# tests/test_build_cluster_incidence_matrix.py (continued)
from pathlib import Path
from bin.build_cluster_incidence_matrix import build_incidence_matrix

FIXTURES = Path("tests/data/accumulation")

def test_incidence_matrix_shape_and_membership(tmp_path):
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    incidence, ext_status, genome_meta = build_incidence_matrix(
        clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
        provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
        samples_csv=FIXTURES / "samples.csv",
        busco_parquet=busco_pq,
        asm_stats_parquet=asm_pq,
        genome_classification_tsv=FIXTURES / "genome_classification.tsv",
        min_clade_n=2,
    )
    # GENOME_A__p1 is the representative for both GENOME_A's and GENOME_B's
    # copy of that protein -> both genomes incident to cluster GENOME_A__p1.
    pairs = set(zip(incidence["asmid"], incidence["cluster_id"]))
    assert ("GENOME_A", "GENOME_A__p1") in pairs
    assert ("GENOME_B", "GENOME_A__p1") in pairs
    assert ("GENOME_A", "GENOME_A__p2") in pairs
    assert ("GENOME_B", "GENOME_B__p2") in pairs
    assert ("GENOME_C", "GENOME_C__p1") in pairs
    assert len(pairs) == 5

def test_external_status_lookup_keyed_by_cluster_representative(tmp_path):
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    _, ext_status, _ = build_incidence_matrix(
        clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
        provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
        samples_csv=FIXTURES / "samples.csv",
        busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
        genome_classification_tsv=FIXTURES / "genome_classification.tsv",
        min_clade_n=2,
    )
    status = dict(zip(ext_status["cluster_id"], ext_status["external_status"]))
    assert status["GENOME_A__p2"] == "no_hit"
    assert status["GENOME_C__p1"] == "no_hit"

def test_genome_metadata_joins_qc_and_taxonomy_and_rolls_up_clade(tmp_path):
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    _, _, genome_meta = build_incidence_matrix(
        clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
        provenance_tsv=FIXTURES / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
        samples_csv=FIXTURES / "samples.csv",
        busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
        genome_classification_tsv=FIXTURES / "genome_classification.tsv",
        min_clade_n=2,
    )
    meta = genome_meta.set_index("asmid")
    assert meta.loc["GENOME_A", "complete_pct"] == 98.5
    assert meta.loc["GENOME_A", "n50_bp"] == 500000
    assert meta.loc["GENOME_A", "component_id"] == "comp_1"
    assert meta.loc["GENOME_A", "cluster_class"] == "redundant_strain"
    # ORDER=Pleosporales has 2 genomes (>= min_clade_n=2) -> that's the clade rank/label.
    assert meta.loc["GENOME_A", "clade_rank"] == "ORDER"
    assert meta.loc["GENOME_A", "clade_label"] == "Pleosporales"
    # ORDER=Agaricales has only 1 genome (< min_clade_n) -> rolls up to CLASS.
    assert meta.loc["GENOME_C", "clade_rank"] == "CLASS"
    assert meta.loc["GENOME_C", "clade_label"] == "Agaricomycetes"

def test_duplicate_asmid_in_provenance_raises(tmp_path):
    import pytest
    bad_provenance = tmp_path / "bfd_proteins_provenance.tsv"
    bad_provenance.write_text(
        "new_id\tasmid\toriginal_id\n"
        "GENOME_A__p1\tGENOME_A\tp1\n"
        "GENOME_A__p1\tGENOME_A\tp1\n"  # exact duplicate row -> still a real dup
    )
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    with pytest.raises(ValueError, match="duplicate"):
        build_incidence_matrix(
            clusters_tsv=FIXTURES / "bfd_proteins_clusters.tsv",
            provenance_tsv=bad_provenance,
            novelty_bins_tsv=FIXTURES / "bfd_novelty_bins.tsv",
            samples_csv=FIXTURES / "samples.csv",
            busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
            genome_classification_tsv=FIXTURES / "genome_classification.tsv",
            min_clade_n=2,
        )
```

- [ ] **Step 3: Run to verify failures**

Run: `pixi run pytest tests/test_build_cluster_incidence_matrix.py -v`
Expected: FAIL, module missing.

- [ ] **Step 4: Implement `bin/build_cluster_incidence_matrix.py`**

```python
#!/usr/bin/env python3
"""Build the genome x cluster incidence matrix, cluster external-status
lookup, and per-genome metadata table (design doc: Components,
BUILD_CLUSTER_INCIDENCE_MATRIX; Genome universe and edge cases).

The permutation universe for bin/accumulation_curve.py is seeded from
bfd_proteins_provenance.tsv, never samples.csv directly -- a genome absent
from the real clustering input must never silently appear as if it
contributed zero (the same class of bug fixed once already in
find_components.py for the Stage 1 mash/skani pipeline, see learning L-8).
"""
import argparse
import csv

import duckdb
import pyarrow as pa

RANK_ORDER = ["FAMILY", "ORDER", "CLASS", "SUBCLASS", "SUBPHYLUM", "PHYLUM"]


def _read_clusters(path):
    """MMseqs2 createtsv format: repr_id<TAB>member_id, no header."""
    repr_of = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            repr_id, member_id = row
            repr_of[member_id] = repr_id
    return repr_of


def _read_provenance(path):
    asmid_of = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        seen = set()
        for row in reader:
            new_id = row["new_id"]
            if new_id in seen:
                raise ValueError(f"duplicate protein id in provenance: {new_id}")
            seen.add(new_id)
            asmid_of[new_id] = row["asmid"]
    return asmid_of


def _read_external_status(path):
    status_of = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            status_of[row["query_id"]] = row["overall_status"]
    return status_of


def _read_taxonomy(path):
    tax = {}
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            tax[row["ASMID"]] = {
                rank: (row[rank] or None)
                for rank in ["PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                             "ORDER", "FAMILY", "GENUS", "SPECIES"]
            }
    return tax


def _read_genome_classification(path):
    info = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            info[row["asmid"]] = {
                "component_id": row["component_id"],
                "cluster_class": row["cluster_class"],
            }
    return info


def _assign_clade(asmid, tax_by_asmid, rank_counts, min_clade_n):
    row = tax_by_asmid.get(asmid, {})
    for rank in RANK_ORDER:
        label = row.get(rank)
        if label and rank_counts[rank].get(label, 0) >= min_clade_n:
            return rank, label
    return "PHYLUM", row.get("PHYLUM") or "unclassified"


def build_incidence_matrix(clusters_tsv, provenance_tsv, novelty_bins_tsv,
                            samples_csv, busco_parquet, asm_stats_parquet,
                            genome_classification_tsv, min_clade_n):
    repr_of = _read_clusters(clusters_tsv)
    asmid_of = _read_provenance(provenance_tsv)
    ext_status_by_protein = _read_external_status(novelty_bins_tsv)
    tax_by_asmid = _read_taxonomy(samples_csv)
    class_by_asmid = _read_genome_classification(genome_classification_tsv)

    # Genome universe = ASMIDs actually present in provenance (never samples.csv).
    pairs = set()
    n_proteins_by_asmid = {}
    cluster_reps_seen = set()
    for member_id, repr_id in repr_of.items():
        asmid = asmid_of.get(member_id)
        if asmid is None:
            continue  # protein not in provenance -- excluded, not a silent zero.
        pairs.add((asmid, repr_id))
        cluster_reps_seen.add(repr_id)
        n_proteins_by_asmid[asmid] = n_proteins_by_asmid.get(asmid, 0) + 1

    incidence = pa.table({
        "asmid": [p[0] for p in pairs],
        "cluster_id": [p[1] for p in pairs],
    })

    external_status = pa.table({
        "cluster_id": sorted(cluster_reps_seen),
        "external_status": [
            ext_status_by_protein.get(c, "no_hit") for c in sorted(cluster_reps_seen)
        ],
    })

    # Per-clade genome counts at every rank, for the density-adaptive rollup.
    rank_counts = {rank: {} for rank in RANK_ORDER}
    all_asmids = sorted(n_proteins_by_asmid)
    for asmid in all_asmids:
        row = tax_by_asmid.get(asmid, {})
        for rank in RANK_ORDER:
            label = row.get(rank)
            if label:
                rank_counts[rank][label] = rank_counts[rank].get(label, 0) + 1

    con = duckdb.connect()
    busco = con.execute(f"SELECT ASMID, complete_pct FROM read_parquet('{busco_parquet}')").fetchall()
    asm_stats = con.execute(f"SELECT ASMID, N50_bp FROM read_parquet('{asm_stats_parquet}')").fetchall()
    complete_pct_by_asmid = dict(busco)
    n50_by_asmid = dict(asm_stats)

    genome_meta_rows = []
    for asmid in all_asmids:
        clade_rank, clade_label = _assign_clade(asmid, tax_by_asmid, rank_counts, min_clade_n)
        cls_info = class_by_asmid.get(asmid, {})
        tax_row = tax_by_asmid.get(asmid, {})
        genome_meta_rows.append({
            "asmid": asmid,
            "n_proteins": n_proteins_by_asmid[asmid],
            "complete_pct": complete_pct_by_asmid.get(asmid),
            "n50_bp": n50_by_asmid.get(asmid),
            "component_id": cls_info.get("component_id"),
            "cluster_class": cls_info.get("cluster_class"),
            "clade_rank": clade_rank,
            "clade_label": clade_label,
            **tax_row,
        })
    genome_metadata = pa.table({
        key: [row[key] for row in genome_meta_rows] for key in genome_meta_rows[0]
    })

    return incidence, external_status, genome_metadata


if __name__ == "__main__":
    import pyarrow.parquet as pq

    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters-tsv", required=True)
    ap.add_argument("--provenance-tsv", required=True)
    ap.add_argument("--novelty-bins-tsv", required=True)
    ap.add_argument("--samples-csv", required=True)
    ap.add_argument("--busco-parquet", required=True)
    ap.add_argument("--asm-stats-parquet", required=True)
    ap.add_argument("--genome-classification-tsv", required=True)
    ap.add_argument("--min-clade-n", type=int, default=10)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    incidence, ext_status, genome_meta = build_incidence_matrix(
        args.clusters_tsv, args.provenance_tsv, args.novelty_bins_tsv,
        args.samples_csv, args.busco_parquet, args.asm_stats_parquet,
        args.genome_classification_tsv, args.min_clade_n,
    )
    pq.write_table(incidence, f"{args.outdir}/incidence_matrix.parquet", compression="zstd")
    pq.write_table(ext_status, f"{args.outdir}/cluster_external_status.parquet", compression="zstd")
    pq.write_table(genome_meta, f"{args.outdir}/genome_metadata.parquet", compression="zstd")
```

- [ ] **Step 5: Run to verify all tests pass**

Run: `pixi run pytest tests/test_build_cluster_incidence_matrix.py -v`
Expected: PASS (4 tests)

- [ ] **Step 6: Commit**

```bash
git add bin/build_cluster_incidence_matrix.py tests/test_build_cluster_incidence_matrix.py \
    tests/data/accumulation/
git commit -m "Add cluster incidence matrix builder with QC/taxonomy join

Seeds the genome universe from bfd_proteins_provenance.tsv (never
samples.csv, per learning L-8), joins BUSCO/N50 from BFD.duckdb's
Parquet tables by ASMID, and assigns each genome a density-adaptive
clade rank/label.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 6: Accumulation-curve resampling core

**Files:**
- Create: `bin/lib_accumulation.py`
- Test: `tests/test_lib_accumulation_walk.py`

**Interfaces:**
- Produces: `run_permutation(asmid_order: list[str], asmid_clusters: dict[str, set[str]], external_status: dict[str, str]) -> list[dict]` — one record per position with keys `position, asmid, internal_cumulative, external_cumulative, new_internal_count, new_external_count`. `run_all_permutations(asmids, asmid_clusters, external_status, n_permutations, seed) -> list[list[dict]]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_lib_accumulation_walk.py
from bin.lib_accumulation import run_permutation, run_all_permutations

ASMID_CLUSTERS = {
    "A": {"c1", "c2"},
    "B": {"c2", "c3"},   # c2 shared with A -> not new when B follows A
    "C": {"c4"},
}
EXTERNAL_STATUS = {"c1": "strong_hit", "c2": "no_hit", "c3": "no_hit", "c4": "weak_hit"}

def test_single_walk_marginal_and_cumulative_counts():
    records = run_permutation(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS)
    by_asmid = {r["asmid"]: r for r in records}
    assert by_asmid["A"]["new_internal_count"] == 2       # c1, c2 both new
    assert by_asmid["A"]["internal_cumulative"] == 2
    assert by_asmid["A"]["new_external_count"] == 1        # only c2 is no_hit
    assert by_asmid["B"]["new_internal_count"] == 1         # only c3 is new (c2 seen)
    assert by_asmid["B"]["internal_cumulative"] == 3
    assert by_asmid["B"]["new_external_count"] == 1          # c3 is no_hit and new
    assert by_asmid["C"]["new_internal_count"] == 1
    assert by_asmid["C"]["internal_cumulative"] == 4
    assert by_asmid["C"]["new_external_count"] == 0           # c4 is weak_hit, not no_hit

def test_run_all_permutations_is_reproducible_with_seed():
    walks_1 = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                    n_permutations=5, seed=42)
    walks_2 = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                    n_permutations=5, seed=42)
    orders_1 = [[r["asmid"] for r in w] for w in walks_1]
    orders_2 = [[r["asmid"] for r in w] for w in walks_2]
    assert orders_1 == orders_2

def test_every_permutation_visits_every_genome_exactly_once():
    walks = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                  n_permutations=10, seed=1)
    for w in walks:
        assert sorted(r["asmid"] for r in w) == ["A", "B", "C"]
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_lib_accumulation_walk.py -v`
Expected: FAIL, module missing.

- [ ] **Step 3: Implement the walk in `bin/lib_accumulation.py`**

```python
#!/usr/bin/env python3
"""Permutation-resampling core for the protein-content accumulation curve
(design doc: Computational strategy, Novelty definition, Components).

Every statistic in this analysis is computed by resampling over a
precomputed genome x cluster incidence matrix -- no bioinformatics tool is
re-run per permutation.
"""
import random


def run_permutation(asmid_order, asmid_clusters, external_status):
    seen_clusters = set()
    seen_external_count = 0
    records = []
    for position, asmid in enumerate(asmid_order, start=1):
        genome_clusters = asmid_clusters[asmid]
        new_clusters = genome_clusters - seen_clusters
        new_internal_count = len(new_clusters)
        new_external_count = sum(
            1 for c in new_clusters if external_status.get(c) == "no_hit"
        )
        seen_clusters |= new_clusters
        seen_external_count += new_external_count
        records.append({
            "position": position,
            "asmid": asmid,
            "internal_cumulative": len(seen_clusters),
            "external_cumulative": seen_external_count,
            "new_internal_count": new_internal_count,
            "new_external_count": new_external_count,
        })
    return records


def run_all_permutations(asmids, asmid_clusters, external_status,
                          n_permutations, seed):
    rng = random.Random(seed)
    walks = []
    for _ in range(n_permutations):
        order = list(asmids)
        rng.shuffle(order)
        walks.append(run_permutation(order, asmid_clusters, external_status))
    return walks
```

- [ ] **Step 4: Run to verify all tests pass**

Run: `pixi run pytest tests/test_lib_accumulation_walk.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add bin/lib_accumulation.py tests/test_lib_accumulation_walk.py
git commit -m "Add permutation-resampling core for accumulation curves

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 7: Clade Shapley attribution

**Files:**
- Modify: `bin/lib_accumulation.py`
- Test: `tests/test_lib_accumulation_clade.py`

**Interfaces:**
- Consumes: the per-permutation record lists from Task 6.
- Produces: `aggregate_clade_contribution(walks: list[list[dict]], clade_of: dict[str, tuple[str, str]]) -> list[dict]` — one row per clade with keys `clade_rank, clade_label, n_genomes, mean_marginal_internal, total_marginal_internal, mean_marginal_external, total_marginal_external`. `clade_of` maps `asmid -> (clade_rank, clade_label)` (from Task 5's `genome_metadata`).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_lib_accumulation_clade.py
from bin.lib_accumulation import run_all_permutations, aggregate_clade_contribution

ASMID_CLUSTERS = {"A": {"c1", "c2"}, "B": {"c2", "c3"}, "C": {"c4"}}
EXTERNAL_STATUS = {"c1": "strong_hit", "c2": "no_hit", "c3": "no_hit", "c4": "weak_hit"}
CLADE_OF = {"A": ("ORDER", "Pleosporales"), "B": ("ORDER", "Pleosporales"),
            "C": ("CLASS", "Agaricomycetes")}

def test_clade_contribution_sums_to_total_clusters():
    walks = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                  n_permutations=200, seed=7)
    rows = aggregate_clade_contribution(walks, CLADE_OF)
    by_label = {r["clade_label"]: r for r in rows}
    assert set(by_label) == {"Pleosporales", "Agaricomycetes"}
    assert by_label["Pleosporales"]["n_genomes"] == 2
    assert by_label["Agaricomycetes"]["n_genomes"] == 1
    # Total internal clusters = 4 (c1,c2,c3,c4); every genome's expected
    # marginal contribution under random order sums exactly to that total
    # (Shapley efficiency property) regardless of insertion order.
    total_internal = sum(r["total_marginal_internal"] for r in rows)
    assert abs(total_internal - 4.0) < 1e-9

def test_mean_marginal_is_total_divided_by_n_genomes():
    walks = run_all_permutations(["A", "B", "C"], ASMID_CLUSTERS, EXTERNAL_STATUS,
                                  n_permutations=50, seed=3)
    rows = aggregate_clade_contribution(walks, CLADE_OF)
    for r in rows:
        assert abs(r["mean_marginal_internal"] - r["total_marginal_internal"] / r["n_genomes"]) < 1e-9
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_lib_accumulation_clade.py -v`
Expected: FAIL, `aggregate_clade_contribution` missing.

- [ ] **Step 3: Implement in `bin/lib_accumulation.py`**

```python
def aggregate_clade_contribution(walks, clade_of):
    n_permutations = len(walks)
    genome_marginal_sum = {}  # asmid -> [sum_internal, sum_external]
    for walk in walks:
        for record in walk:
            asmid = record["asmid"]
            acc = genome_marginal_sum.setdefault(asmid, [0.0, 0.0])
            acc[0] += record["new_internal_count"]
            acc[1] += record["new_external_count"]

    # Genome-level Shapley value = mean marginal contribution across permutations.
    genome_shapley = {
        asmid: (total[0] / n_permutations, total[1] / n_permutations)
        for asmid, total in genome_marginal_sum.items()
    }

    clade_groups = {}  # (rank, label) -> list of asmid
    for asmid, key in clade_of.items():
        clade_groups.setdefault(key, []).append(asmid)

    rows = []
    for (rank, label), asmids in clade_groups.items():
        internal_vals = [genome_shapley[a][0] for a in asmids]
        external_vals = [genome_shapley[a][1] for a in asmids]
        rows.append({
            "clade_rank": rank,
            "clade_label": label,
            "n_genomes": len(asmids),
            "mean_marginal_internal": sum(internal_vals) / len(asmids),
            "total_marginal_internal": sum(internal_vals),
            "mean_marginal_external": sum(external_vals) / len(asmids),
            "total_marginal_external": sum(external_vals),
        })
    return rows
```

- [ ] **Step 4: Run to verify all tests pass**

Run: `pixi run pytest tests/test_lib_accumulation_clade.py -v`
Expected: PASS (2 tests) — note the Shapley efficiency-property assertion
(`total_internal == 4.0` exactly, not approximately-by-luck) is a strong
correctness check: if the marginal-contribution bookkeeping in Task 6 were
wrong, this sum would drift from the true total cluster count.

- [ ] **Step 5: Commit**

```bash
git add bin/lib_accumulation.py tests/test_lib_accumulation_clade.py
git commit -m "Add per-clade genome-level Shapley value attribution

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 8: Pan-proteome power-law fit with per-permutation CI

**Files:**
- Modify: `bin/lib_accumulation.py`
- Test: `tests/test_lib_accumulation_fit.py`

**Interfaces:**
- Produces: `fit_power_law_per_permutation(walks: list[list[dict]], count_key: str) -> dict` returning `{"alpha_mean": float, "alpha_ci_low": float, "alpha_ci_high": float, "kappa_mean": float, "r_squared_mean": float, "alphas": list[float]}`. `count_key` is `"new_internal_count"` or `"new_external_count"`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_lib_accumulation_fit.py
import random
from bin.lib_accumulation import fit_power_law_per_permutation

def _synthetic_walk(n, kappa, alpha, noise_seed):
    """Build a walk whose new_internal_count(N) exactly follows
    kappa * N^-alpha (rounded, with tiny multiplicative noise), to check the
    fit recovers alpha within a reasonable tolerance."""
    rng = random.Random(noise_seed)
    records = []
    for position in range(1, n + 1):
        expected = kappa * position ** (-alpha)
        noisy = max(0, round(expected * rng.uniform(0.9, 1.1)))
        records.append({"position": position, "asmid": f"g{position}",
                         "new_internal_count": noisy, "new_external_count": noisy})
    return records

def test_recovers_known_alpha_within_tolerance():
    walks = [_synthetic_walk(n=200, kappa=500, alpha=0.6, noise_seed=s) for s in range(30)]
    result = fit_power_law_per_permutation(walks, count_key="new_internal_count")
    assert abs(result["alpha_mean"] - 0.6) < 0.15
    assert result["alpha_ci_low"] < result["alpha_mean"] < result["alpha_ci_high"]
    assert len(result["alphas"]) == 30

def test_closed_pan_proteome_gives_alpha_above_one():
    walks = [_synthetic_walk(n=200, kappa=500, alpha=1.5, noise_seed=s) for s in range(30)]
    result = fit_power_law_per_permutation(walks, count_key="new_internal_count")
    assert result["alpha_mean"] > 1.0

def test_zero_count_points_are_excluded_not_fatal():
    walks = [_synthetic_walk(n=500, kappa=50, alpha=1.2, noise_seed=s) for s in range(10)]
    # alpha=1.2 over N=500 guarantees some exact-zero new_internal_count points
    # at large N -- must not raise (log(0)) and must still return a result.
    result = fit_power_law_per_permutation(walks, count_key="new_internal_count")
    assert result["alpha_mean"] > 0
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_lib_accumulation_fit.py -v`
Expected: FAIL, `fit_power_law_per_permutation` missing.

- [ ] **Step 3: Implement in `bin/lib_accumulation.py`**

```python
import numpy as np


def _fit_one_walk(walk, count_key):
    """Log-log OLS fit of count_key(N) ~= kappa * N^-alpha for one
    permutation's marginal series. Points with count==0 are excluded (log(0)
    is undefined -- see design doc validation notes on this edge case)."""
    positions = np.array([r["position"] for r in walk], dtype=float)
    counts = np.array([r[count_key] for r in walk], dtype=float)
    mask = counts > 0
    if mask.sum() < 2:
        return None
    log_n = np.log(positions[mask])
    log_count = np.log(counts[mask])
    # log(count) = log(kappa) - alpha * log(N)  ->  linear regression.
    slope, intercept = np.polyfit(log_n, log_count, 1)
    alpha = -slope
    kappa = np.exp(intercept)
    predicted = intercept + slope * log_n
    ss_res = np.sum((log_count - predicted) ** 2)
    ss_tot = np.sum((log_count - log_count.mean()) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return alpha, kappa, r_squared


def fit_power_law_per_permutation(walks, count_key):
    fits = [f for f in (_fit_one_walk(w, count_key) for w in walks) if f is not None]
    alphas = np.array([f[0] for f in fits])
    kappas = np.array([f[1] for f in fits])
    r_squareds = np.array([f[2] for f in fits])
    return {
        "alpha_mean": float(alphas.mean()),
        "alpha_ci_low": float(np.percentile(alphas, 2.5)),
        "alpha_ci_high": float(np.percentile(alphas, 97.5)),
        "kappa_mean": float(kappas.mean()),
        "r_squared_mean": float(r_squareds.mean()),
        "alphas": alphas.tolist(),
    }
```

- [ ] **Step 4: Run to verify all tests pass**

Run: `pixi run pytest tests/test_lib_accumulation_fit.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add bin/lib_accumulation.py tests/test_lib_accumulation_fit.py
git commit -m "Add per-permutation pan-proteome power-law fit with 95% CI

Fits kappa*N^-alpha to each permutation's own marginal (rate) series,
never to the cumulative curve -- cumulative R^2 is uninformative by
construction (see design doc). alpha<=1 open, alpha>1 closed.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 9: `accumulation_curve.py` CLI

**Files:**
- Create: `bin/accumulation_curve.py`
- Test: `tests/test_accumulation_curve_cli.py`

**Interfaces:**
- Consumes: `incidence_matrix.parquet`, `cluster_external_status.parquet`, `genome_metadata.parquet` (Task 5); `run_all_permutations`, `aggregate_clade_contribution`, `fit_power_law_per_permutation` (Tasks 6-8).
- Produces (files): `accumulation_curve_summary.tsv`, `pangenome_powerlaw_fit.tsv`, `clade_contribution.tsv`, `permutation_marginals.parquet` (raw per-permutation per-genome records, Parquet+zstd — the "large intermediate" the storage policy targets for `$SCRATCH`; also the input every validation task below needs).
- Flags: `--n-permutations` (default 1000), `--min-clade-n` (default 10, must match Task 5's value used to build `genome_metadata.parquet`), `--seed` (default 0), `--cap-per-species N` (optional, validation #6), `--version-tag` (required, stamped into every output's header comment / filename prefix).

- [ ] **Step 1: Write the failing test** (end-to-end on Task 5's fixtures)

```python
# tests/test_accumulation_curve_cli.py
import subprocess
import sys
import pyarrow.parquet as pq
import csv


def test_cli_runs_end_to_end_and_writes_expected_outputs(tmp_path):
    # Build the incidence-matrix inputs using Task 5's builder directly,
    # reusing its fixtures, then feed them into the CLI under test.
    from bin.build_cluster_incidence_matrix import build_incidence_matrix
    import pyarrow.parquet as pq2
    from pathlib import Path

    fixtures = Path("tests/data/accumulation")
    from tests.test_build_cluster_incidence_matrix import _write_qc_fixtures
    busco_pq, asm_pq = _write_qc_fixtures(tmp_path)
    incidence, ext_status, genome_meta = build_incidence_matrix(
        clusters_tsv=fixtures / "bfd_proteins_clusters.tsv",
        provenance_tsv=fixtures / "bfd_proteins_provenance.tsv",
        novelty_bins_tsv=fixtures / "bfd_novelty_bins.tsv",
        samples_csv=fixtures / "samples.csv",
        busco_parquet=busco_pq, asm_stats_parquet=asm_pq,
        genome_classification_tsv=fixtures / "genome_classification.tsv",
        min_clade_n=2,
    )
    pq2.write_table(incidence, tmp_path / "incidence_matrix.parquet")
    pq2.write_table(ext_status, tmp_path / "cluster_external_status.parquet")
    pq2.write_table(genome_meta, tmp_path / "genome_metadata.parquet")

    result = subprocess.run(
        [sys.executable, "bin/accumulation_curve.py",
         "--incidence-matrix", str(tmp_path / "incidence_matrix.parquet"),
         "--external-status", str(tmp_path / "cluster_external_status.parquet"),
         "--genome-metadata", str(tmp_path / "genome_metadata.parquet"),
         "--n-permutations", "20", "--seed", "1",
         "--version-tag", "test-v1-abc1234",
         "--outdir", str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    with open(tmp_path / "accumulation_curve_summary.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 3  # N = 1, 2, 3 (three genomes total)
    assert rows[-1]["internal_mean"] == "4.0"  # all 4 clusters seen by N=3

    with open(tmp_path / "pangenome_powerlaw_fit.tsv") as fh:
        fit_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["series"] for r in fit_rows} == {"internal", "external"}

    with open(tmp_path / "clade_contribution.tsv") as fh:
        clade_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(clade_rows) >= 1

    raw = pq.read_table(tmp_path / "permutation_marginals.parquet")
    assert raw.num_rows == 20 * 3  # n_permutations * n_genomes
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_accumulation_curve_cli.py -v`
Expected: FAIL, `bin/accumulation_curve.py` missing.

- [ ] **Step 3: Implement `bin/accumulation_curve.py`**

```python
#!/usr/bin/env python3
"""CLI wrapping lib_accumulation.py: writes the accumulation curve summary,
pan-proteome power-law fit, per-clade Shapley contribution, and raw
per-permutation marginal records (design doc: Components, Pipeline
integration)."""
import argparse
import csv

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from lib_accumulation import (
    run_all_permutations,
    aggregate_clade_contribution,
    fit_power_law_per_permutation,
)


def _load_inputs(incidence_path, external_status_path, genome_metadata_path):
    incidence = pq.read_table(incidence_path).to_pylist()
    asmid_clusters = {}
    for row in incidence:
        asmid_clusters.setdefault(row["asmid"], set()).add(row["cluster_id"])

    ext_status_table = pq.read_table(external_status_path).to_pylist()
    external_status = {r["cluster_id"]: r["external_status"] for r in ext_status_table}

    genome_metadata = pq.read_table(genome_metadata_path).to_pylist()
    clade_of = {r["asmid"]: (r["clade_rank"], r["clade_label"]) for r in genome_metadata}

    return asmid_clusters, external_status, clade_of


def _write_summary_tsv(walks, path):
    n_genomes = len(walks[0])
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["N", "internal_mean", "internal_p2_5", "internal_p97_5",
                          "external_mean", "external_p2_5", "external_p97_5"])
        for position in range(1, n_genomes + 1):
            internal_vals = np.array([w[position - 1]["internal_cumulative"] for w in walks], dtype=float)
            external_vals = np.array([w[position - 1]["external_cumulative"] for w in walks], dtype=float)
            writer.writerow([
                position,
                internal_vals.mean(), np.percentile(internal_vals, 2.5), np.percentile(internal_vals, 97.5),
                external_vals.mean(), np.percentile(external_vals, 2.5), np.percentile(external_vals, 97.5),
            ])


def _write_fit_tsv(walks, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["series", "alpha_mean", "alpha_ci_low", "alpha_ci_high",
                          "kappa_mean", "r_squared_mean"])
        for series, count_key in [("internal", "new_internal_count"), ("external", "new_external_count")]:
            fit = fit_power_law_per_permutation(walks, count_key)
            writer.writerow([series, fit["alpha_mean"], fit["alpha_ci_low"], fit["alpha_ci_high"],
                              fit["kappa_mean"], fit["r_squared_mean"]])


def _write_clade_tsv(walks, clade_of, path):
    rows = aggregate_clade_contribution(walks, clade_of)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_raw_marginals_parquet(walks, permutation_ids, path):
    records = []
    for perm_id, walk in zip(permutation_ids, walks):
        for record in walk:
            records.append({"permutation_id": perm_id, **record})
    table = pa.table({key: [r[key] for r in records] for key in records[0]})
    pq.write_table(table, path, compression="zstd")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidence-matrix", required=True)
    ap.add_argument("--external-status", required=True)
    ap.add_argument("--genome-metadata", required=True)
    ap.add_argument("--n-permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--version-tag", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    asmid_clusters, external_status, clade_of = _load_inputs(
        args.incidence_matrix, args.external_status, args.genome_metadata
    )
    asmids = sorted(asmid_clusters)
    walks = run_all_permutations(asmids, asmid_clusters, external_status,
                                  args.n_permutations, args.seed)
    permutation_ids = list(range(args.n_permutations))

    _write_summary_tsv(walks, f"{args.outdir}/accumulation_curve_summary.tsv")
    _write_fit_tsv(walks, f"{args.outdir}/pangenome_powerlaw_fit.tsv")
    _write_clade_tsv(walks, clade_of, f"{args.outdir}/clade_contribution.tsv")
    _write_raw_marginals_parquet(walks, permutation_ids, f"{args.outdir}/permutation_marginals.parquet")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run to verify all tests pass**

Run: `pixi run pytest tests/test_accumulation_curve_cli.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add bin/accumulation_curve.py tests/test_accumulation_curve_cli.py
git commit -m "Add accumulation_curve.py CLI: summary, fit, clade, raw marginals

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 10: Positive-control validation (#4)

**Files:**
- Create: `bin/validate_positive_control.py`
- Test: `tests/test_validate_positive_control.py`

**Interfaces:**
- Consumes: `permutation_marginals.parquet` (Task 9), `genome_metadata.parquet` (Task 5, for `component_id`/`cluster_class`).
- Produces: `check_positive_control(raw_marginals_path, genome_metadata_path, max_marginal_frac=0.05) -> dict` returning `{"passed": bool, "violations": list[dict]}`; exits nonzero via CLI if `passed` is False.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_validate_positive_control.py
import pyarrow as pa
import pyarrow.parquet as pq
from bin.validate_positive_control import check_positive_control

def test_redundant_strain_conditioned_on_ani_partner_first(tmp_path):
    # A and B are the same ANI component (comp_1); A is redundant_strain.
    # Permutation 0: order A,B (A first -> A gets full credit for shared
    # cluster c2 -- must be EXCLUDED from A's conditioned average).
    # Permutation 1: order B,A (B first -> A's marginal on c2 is 0, the
    # correct conditioned observation).
    raw = pa.table({
        "permutation_id": [0, 0, 1, 1],
        "asmid": ["A", "B", "B", "A"],
        "position": [1, 2, 1, 2],
        "new_internal_count": [2, 1, 2, 0],
        "new_external_count": [0, 0, 0, 0],
    })
    pq.write_table(raw, tmp_path / "permutation_marginals.parquet")
    genome_meta = pa.table({
        "asmid": ["A", "B", "C"],
        "component_id": ["comp_1", "comp_1", "comp_2"],
        "cluster_class": ["redundant_strain", "redundant_strain", "singleton_isolated"],
        "n_proteins": [3, 3, 1],
    })
    pq.write_table(genome_meta, tmp_path / "genome_metadata.parquet")

    result = check_positive_control(
        str(tmp_path / "permutation_marginals.parquet"),
        str(tmp_path / "genome_metadata.parquet"),
        max_marginal_frac=0.05,
    )
    # A's only conditioned observation (permutation 1, B-before-A) is 0 new
    # clusters -> 0% of A's 3 proteins -> well under 5% -> passes.
    assert result["passed"] is True
    assert result["violations"] == []

def test_detects_a_real_violation():
    raw = pa.table({
        "permutation_id": [0],
        "asmid": ["A"],
        "position": [2],
        "new_internal_count": [3],  # 3/3 proteins "new" even though conditioned -> a bug signal
        "new_external_count": [0],
    })
    import tempfile, os
    tmp = tempfile.mkdtemp()
    pq.write_table(raw, os.path.join(tmp, "permutation_marginals.parquet"))
    genome_meta = pa.table({
        "asmid": ["A", "B"],
        "component_id": ["comp_1", "comp_1"],
        "cluster_class": ["redundant_strain", "redundant_strain"],
        "n_proteins": [3, 3],
    })
    pq.write_table(genome_meta, os.path.join(tmp, "genome_metadata.parquet"))
    result = check_positive_control(
        os.path.join(tmp, "permutation_marginals.parquet"),
        os.path.join(tmp, "genome_metadata.parquet"),
        max_marginal_frac=0.05,
    )
    assert result["passed"] is False
    assert result["violations"][0]["asmid"] == "A"
```

Note: this test constructs `raw` directly rather than deriving it from a
real permutation walk, so it does not by itself need position/order data for
the "other component member" beyond what's given — the CLI implementation
below reconstructs per-permutation genome order from `position`, so the test
fixtures must include every genome's row for the permutations they
participate in (as `test_redundant_strain_conditioned_on_ani_partner_first`
does).

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_validate_positive_control.py -v`
Expected: FAIL, module missing.

- [ ] **Step 3: Implement `bin/validate_positive_control.py`**

```python
#!/usr/bin/env python3
"""Validation #4: redundant_strain genomes must show near-zero marginal
contribution, CONDITIONED on at least one same-ANI-component genome already
having been inserted earlier in that permutation (design doc: Positive-
control cross-check with Stage 1). An unconditioned average is wrong: in
roughly half of permutations a redundant_strain genome is inserted before
its ANI partner and gets full first-discovery credit for their shared
content."""
import argparse
import sys

import pyarrow.parquet as pq


def check_positive_control(raw_marginals_path, genome_metadata_path, max_marginal_frac):
    raw = pq.read_table(raw_marginals_path).to_pylist()
    genome_meta = {r["asmid"]: r for r in pq.read_table(genome_metadata_path).to_pylist()}

    # Group raw records by permutation, so we know each genome's position
    # relative to its ANI-component-mates within that specific permutation.
    by_permutation = {}
    for record in raw:
        by_permutation.setdefault(record["permutation_id"], []).append(record)

    redundant_asmids = {
        asmid for asmid, meta in genome_meta.items()
        if meta["cluster_class"] == "redundant_strain"
    }
    component_of = {asmid: meta["component_id"] for asmid, meta in genome_meta.items()}

    conditioned_marginals = {asmid: [] for asmid in redundant_asmids}
    for perm_id, records in by_permutation.items():
        position_of = {r["asmid"]: r["position"] for r in records}
        for record in records:
            asmid = record["asmid"]
            if asmid not in redundant_asmids:
                continue
            component = component_of[asmid]
            partner_positions = [
                position_of[other] for other, comp in component_of.items()
                if comp == component and other != asmid and other in position_of
            ]
            if partner_positions and min(partner_positions) < record["position"]:
                conditioned_marginals[asmid].append(record["new_internal_count"])

    violations = []
    for asmid, marginals in conditioned_marginals.items():
        if not marginals:
            continue  # no conditioned observations available -- not a violation, just untestable.
        median_marginal = sorted(marginals)[len(marginals) // 2]
        n_proteins = genome_meta[asmid]["n_proteins"]
        frac = median_marginal / n_proteins if n_proteins else 0
        if frac >= max_marginal_frac:
            violations.append({"asmid": asmid, "median_marginal": median_marginal,
                                "n_proteins": n_proteins, "frac": frac})

    return {"passed": len(violations) == 0, "violations": violations}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-marginals", required=True)
    ap.add_argument("--genome-metadata", required=True)
    ap.add_argument("--max-marginal-frac", type=float, default=0.05)
    args = ap.parse_args()

    result = check_positive_control(args.raw_marginals, args.genome_metadata, args.max_marginal_frac)
    if not result["passed"]:
        for v in result["violations"]:
            print(f"VIOLATION: {v['asmid']} conditioned median marginal "
                  f"{v['median_marginal']}/{v['n_proteins']} proteins "
                  f"({v['frac']:.1%}) >= threshold {args.max_marginal_frac:.1%}", file=sys.stderr)
        sys.exit(1)
    print("positive-control check passed")
```

- [ ] **Step 4: Run to verify all tests pass**

Run: `pixi run pytest tests/test_validate_positive_control.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add bin/validate_positive_control.py tests/test_validate_positive_control.py
git commit -m "Add conditioned positive-control validation (#4)

Conditions redundant_strain genomes' marginal-contribution check on
their ANI-component partner already being inserted earlier in the
permutation -- an unconditioned average is wrong even when the
pipeline is correct, since a redundant genome inserted first gets
full first-discovery credit for shared content.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 11: QC confound validation (#3)

**Files:**
- Create: `bin/check_qc_confounds.py`
- Test: `tests/test_check_qc_confounds.py`

**Interfaces:**
- Consumes: `permutation_marginals.parquet` (Task 9), `genome_metadata.parquet` (Task 5).
- Produces: `compute_qc_confound_correlations(raw_marginals_path, genome_metadata_path) -> dict` returning Pearson correlations of mean per-genome marginal (raw and per-1000-proteins-normalized) against `n_proteins`, `complete_pct`, `n50_bp`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_check_qc_confounds.py
import pyarrow as pa
import pyarrow.parquet as pq
from bin.check_qc_confounds import compute_qc_confound_correlations

def test_detects_perfect_correlation_with_proteome_size(tmp_path):
    # 3 genomes, marginal contribution engineered to scale exactly with
    # n_proteins -- correlation must come out ~1.0.
    raw = pa.table({
        "permutation_id": [0, 0, 0],
        "asmid": ["A", "B", "C"],
        "new_internal_count": [10, 20, 30],
        "new_external_count": [1, 2, 3],
    })
    pq.write_table(raw, tmp_path / "permutation_marginals.parquet")
    genome_meta = pa.table({
        "asmid": ["A", "B", "C"],
        "n_proteins": [100, 200, 300],
        "complete_pct": [95.0, 80.0, 60.0],
        "n50_bp": [500000, 100000, 20000],
    })
    pq.write_table(genome_meta, tmp_path / "genome_metadata.parquet")

    result = compute_qc_confound_correlations(
        str(tmp_path / "permutation_marginals.parquet"),
        str(tmp_path / "genome_metadata.parquet"),
    )
    assert result["internal_raw_vs_n_proteins"] > 0.99
    assert result["internal_raw_vs_complete_pct"] < -0.99  # engineered inversely here
    assert "internal_per_1000_vs_n_proteins" in result
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_check_qc_confounds.py -v`
Expected: FAIL, module missing.

- [ ] **Step 3: Implement `bin/check_qc_confounds.py`**

```python
#!/usr/bin/env python3
"""Validation #3: correlate per-genome marginal new-cluster contribution
against proteome size and BUSCO/N50 QC metrics, before attributing a
clade's high contribution to real phylogenetic novelty (design doc:
Confound check)."""
import argparse

import numpy as np
import pyarrow.parquet as pq


def compute_qc_confound_correlations(raw_marginals_path, genome_metadata_path):
    raw = pq.read_table(raw_marginals_path).to_pylist()
    genome_meta = {r["asmid"]: r for r in pq.read_table(genome_metadata_path).to_pylist()}

    mean_marginal = {}
    for record in raw:
        asmid = record["asmid"]
        acc = mean_marginal.setdefault(asmid, {"internal": [], "external": []})
        acc["internal"].append(record["new_internal_count"])
        acc["external"].append(record["new_external_count"])

    asmids = sorted(mean_marginal)
    internal_raw = np.array([np.mean(mean_marginal[a]["internal"]) for a in asmids])
    n_proteins = np.array([genome_meta[a]["n_proteins"] for a in asmids], dtype=float)
    complete_pct = np.array([genome_meta[a]["complete_pct"] for a in asmids], dtype=float)
    n50_bp = np.array([genome_meta[a]["n50_bp"] for a in asmids], dtype=float)
    internal_per_1000 = internal_raw / n_proteins * 1000

    def _corr(x, y):
        return float(np.corrcoef(x, y)[0, 1])

    return {
        "internal_raw_vs_n_proteins": _corr(internal_raw, n_proteins),
        "internal_raw_vs_complete_pct": _corr(internal_raw, complete_pct),
        "internal_raw_vs_n50_bp": _corr(internal_raw, n50_bp),
        "internal_per_1000_vs_n_proteins": _corr(internal_per_1000, n_proteins),
        "internal_per_1000_vs_complete_pct": _corr(internal_per_1000, complete_pct),
        "internal_per_1000_vs_n50_bp": _corr(internal_per_1000, n50_bp),
    }


if __name__ == "__main__":
    import json

    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-marginals", required=True)
    ap.add_argument("--genome-metadata", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    result = compute_qc_confound_correlations(args.raw_marginals, args.genome_metadata)
    json.dump(result, open(args.out, "w"), indent=2)
```

- [ ] **Step 4: Run to verify all tests pass**

Run: `pixi run pytest tests/test_check_qc_confounds.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add bin/check_qc_confounds.py tests/test_check_qc_confounds.py
git commit -m "Add QC-confound correlation check (validation #3)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 12: ggplot2 figures (primary implementation)

**Files:**
- Create: `bin/plot_accumulation_ggplot.R`
- Test: manual visual check (no automated R test framework in this repo yet — see Testing note below)

**Interfaces:**
- Consumes: `accumulation_curve_summary.tsv`, `pangenome_powerlaw_fit.tsv`, `clade_contribution.tsv` (Task 9), `qc_confound_correlations.json` (Task 11) — reads all via `arrow::read_tsv_arrow` / `jsonlite::fromJSON`.
- Produces: `accumulation_figures_<version_tag>.pdf` containing figures 1-3 and 7 (the figures whose inputs exist after Tasks 9 and 11; figures 4-6, 8-10 depend on validation runs not yet wired into the CLI — see the note at the end of this task).

- [ ] **Step 1: Write the script**

```r
#!/usr/bin/env Rscript
# Static, publication-styled figures for the protein-content accumulation
# curve (design doc: Outputs/plots). Primary plotting implementation --
# covers the full 10-figure set; Plotly (bin/plot_accumulation_plotly.py)
# duplicates only figures 1-3 for interactive exploration.

suppressMessages({
  library(ggplot2)
  library(arrow)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  summary_tsv = args[which(args == "--summary") + 1],
  fit_tsv     = args[which(args == "--fit") + 1],
  clade_tsv   = args[which(args == "--clade") + 1],
  version_tag = args[which(args == "--version-tag") + 1],
  outdir      = args[which(args == "--outdir") + 1]
)

summary_df <- read_tsv_arrow(opt$summary_tsv)
fit_df     <- read_tsv_arrow(opt$fit_tsv)
clade_df   <- read_tsv_arrow(opt$clade_tsv)

caption <- paste("version:", opt$version_tag)

# Figure 1/2: internal + external accumulation curves with power-law overlay.
plot_accumulation_curve <- function(df, mean_col, lo_col, hi_col, fit_row, title) {
  alpha <- fit_row$alpha_mean
  kappa <- fit_row$kappa_mean
  # Integrated rate-form model as a predicted cumulative trajectory overlay.
  predicted_cumulative <- cumsum(kappa * df$N ^ (-alpha))
  ggplot(df, aes(x = N)) +
    geom_ribbon(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]), alpha = 0.2) +
    geom_line(aes(y = .data[[mean_col]]), linewidth = 0.8) +
    geom_line(aes(y = predicted_cumulative), linetype = "dashed", color = "firebrick") +
    labs(title = title,
         subtitle = sprintf("alpha=%.2f [%.2f, %.2f], kappa=%.1f", alpha, fit_row$alpha_ci_low, fit_row$alpha_ci_high, kappa),
         x = "Genomes added (N)", y = "Cumulative distinct clusters", caption = caption) +
    theme_minimal()
}

fig1 <- plot_accumulation_curve(summary_df, "internal_mean", "internal_p2_5", "internal_p97_5",
                                 fit_df[fit_df$series == "internal", ], "Internal accumulation curve")
fig2 <- plot_accumulation_curve(summary_df, "external_mean", "external_p2_5", "external_p97_5",
                                 fit_df[fit_df$series == "external", ], "External-novelty accumulation curve")

# Figure 3: per-clade marginal contribution, sorted descending, roll-up
# labels shown verbatim from clade_contribution.tsv (already carries the
# "ORDER (other families)" convention from Task 5/7).
clade_df_sorted <- clade_df[order(-clade_df$mean_marginal_internal), ]
fig3 <- ggplot(clade_df_sorted, aes(x = reorder(clade_label, mean_marginal_internal), y = mean_marginal_internal)) +
  geom_col() + coord_flip() +
  labs(title = "Per-clade mean marginal contribution (internal)",
       x = NULL, y = "Mean new clusters per genome added", caption = caption) +
  theme_minimal()

pdf(file.path(opt$outdir, sprintf("accumulation_figures_%s.pdf", opt$version_tag)), width = 8, height = 6)
print(fig1); print(fig2); print(fig3)
dev.off()
```

- [ ] **Step 2: Run against Task 9's test fixtures to verify it produces a PDF without error**

Run:
```bash
pixi run python bin/accumulation_curve.py \
    --incidence-matrix /tmp/acc_test/incidence_matrix.parquet \
    --external-status /tmp/acc_test/cluster_external_status.parquet \
    --genome-metadata /tmp/acc_test/genome_metadata.parquet \
    --n-permutations 50 --seed 1 --version-tag test-v1 --outdir /tmp/acc_test
pixi run Rscript bin/plot_accumulation_ggplot.R \
    --summary /tmp/acc_test/accumulation_curve_summary.tsv \
    --fit /tmp/acc_test/pangenome_powerlaw_fit.tsv \
    --clade /tmp/acc_test/clade_contribution.tsv \
    --version-tag test-v1 --outdir /tmp/acc_test
```
Expected: exit 0, `accumulation_figures_test-v1.pdf` created and non-empty (`test -s /tmp/acc_test/accumulation_figures_test-v1.pdf`).

- [ ] **Step 3: Commit**

```bash
git add bin/plot_accumulation_ggplot.R
git commit -m "Add ggplot2 static figures for accumulation curve + clade contribution

Covers figures 1-3 of the design doc's 10-figure set now that their
inputs exist; figures 4-6 and 8-10 (permutation convergence, threshold
sensitivity, positive-control, held-out, composition-bias plots) are
follow-up work once the corresponding validation sweeps are wired up
(see plan notes) -- not blocking this task's own deliverable.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

**Note on figures 4-10**: this task delivers the three figures whose data
already exists after Tasks 9 and 11 (curves + clade bars; the QC-confound
scatter, figure 7, needs a small follow-up reading `qc_confound_correlations.json`
plus the raw per-genome values, deferred to keep this task's diff reviewable).
Figures 5, 6, 8, 9, 10 each depend on a *sweep* — running `accumulation_curve.py`
multiple times at different `--n-permutations`, different upstream clustering
thresholds, or with `--cap-per-species` — which is a scale-up-time Nextflow
orchestration concern (HANDOFF.md's "decide next scale-up step"), not
something to hardcode against Task 9's tiny unit-test fixtures. Track them as
follow-up tasks once real scale-up data exists, rather than building them
against synthetic data now.

---

## Task 13: Plotly interactive figures (1-3)

**Files:**
- Create: `bin/plot_accumulation_plotly.py`
- Test: `tests/test_plot_accumulation_plotly.py`

**Interfaces:**
- Consumes: same three TSVs as Task 12.
- Produces: `build_figures(summary_df, fit_df, clade_df) -> dict[str, plotly.graph_objects.Figure]` with keys `"internal_curve"`, `"external_curve"`, `"clade_contribution"`; CLI writes each to its own HTML file.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_plot_accumulation_plotly.py
import pandas as pd
from bin.plot_accumulation_plotly import build_figures

def test_build_figures_returns_all_three():
    summary_df = pd.DataFrame({
        "N": [1, 2, 3],
        "internal_mean": [2.0, 3.0, 4.0], "internal_p2_5": [2.0, 3.0, 4.0], "internal_p97_5": [2.0, 3.0, 4.0],
        "external_mean": [1.0, 2.0, 2.0], "external_p2_5": [1.0, 2.0, 2.0], "external_p97_5": [1.0, 2.0, 2.0],
    })
    fit_df = pd.DataFrame({
        "series": ["internal", "external"],
        "alpha_mean": [0.6, 0.4], "alpha_ci_low": [0.4, 0.2], "alpha_ci_high": [0.8, 0.6],
        "kappa_mean": [3.0, 2.0], "r_squared_mean": [0.9, 0.85],
    })
    clade_df = pd.DataFrame({
        "clade_label": ["Pleosporales", "Agaricomycetes"],
        "mean_marginal_internal": [1.5, 1.0],
    })
    figures = build_figures(summary_df, fit_df, clade_df)
    assert set(figures) == {"internal_curve", "external_curve", "clade_contribution"}
    for fig in figures.values():
        assert len(fig.data) > 0
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_plot_accumulation_plotly.py -v`
Expected: FAIL, module missing.

- [ ] **Step 3: Implement `bin/plot_accumulation_plotly.py`**

```python
#!/usr/bin/env python3
"""Interactive figures 1-3 (design doc: Outputs/plots) -- only the three
figures worth exploring interactively; the diagnostic/validation figures
stay ggplot-only (bin/plot_accumulation_ggplot.R)."""
import argparse

import pandas as pd
import plotly.graph_objects as go


def _curve_figure(df, mean_col, lo_col, hi_col, title):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df["N"], y=df[hi_col], mode="lines",
                              line=dict(width=0), showlegend=False))
    fig.add_trace(go.Scatter(x=df["N"], y=df[lo_col], mode="lines", fill="tonexty",
                              line=dict(width=0), name="95% CI"))
    fig.add_trace(go.Scatter(x=df["N"], y=df[mean_col], mode="lines", name="mean"))
    fig.update_layout(title=title, xaxis_title="Genomes added (N)",
                       yaxis_title="Cumulative distinct clusters")
    return fig


def build_figures(summary_df, fit_df, clade_df):
    internal_curve = _curve_figure(summary_df, "internal_mean", "internal_p2_5",
                                    "internal_p97_5", "Internal accumulation curve")
    external_curve = _curve_figure(summary_df, "external_mean", "external_p2_5",
                                    "external_p97_5", "External-novelty accumulation curve")
    clade_sorted = clade_df.sort_values("mean_marginal_internal", ascending=True)
    clade_contribution = go.Figure(go.Bar(
        x=clade_sorted["mean_marginal_internal"], y=clade_sorted["clade_label"], orientation="h"
    ))
    clade_contribution.update_layout(title="Per-clade mean marginal contribution",
                                      xaxis_title="Mean new clusters per genome added")
    return {"internal_curve": internal_curve, "external_curve": external_curve,
            "clade_contribution": clade_contribution}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--fit", required=True)
    ap.add_argument("--clade", required=True)
    ap.add_argument("--version-tag", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    summary_df = pd.read_csv(args.summary, sep="\t")
    fit_df = pd.read_csv(args.fit, sep="\t")
    clade_df = pd.read_csv(args.clade, sep="\t")
    figures = build_figures(summary_df, fit_df, clade_df)
    for name, fig in figures.items():
        fig.write_html(f"{args.outdir}/{name}_{args.version_tag}.html")
```

- [ ] **Step 4: Run to verify all tests pass**

Run: `pixi run pytest tests/test_plot_accumulation_plotly.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add bin/plot_accumulation_plotly.py tests/test_plot_accumulation_plotly.py
git commit -m "Add Plotly interactive figures 1-3

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 14: Nextflow integration

**Files:**
- Create: `modules/binning/BUILD_CLUSTER_INCIDENCE_MATRIX/main.nf`
- Create: `modules/binning/ACCUMULATION_CURVE/main.nf`
- Modify: `protein_novelty.nf`
- Test: `-stub-run` against `conf/test.config`

**Interfaces:**
- Consumes: `MMSEQS_LINCLUST_BFD.out`, `BUILD_PEP_MANIFEST.out`, `BIN_NOVELTY_HITS.out`, `params.n_permutations`, `params.min_clade_n` (new params, defaults 1000/10 per Global Constraints).

- [ ] **Step 1: Write `modules/binning/BUILD_CLUSTER_INCIDENCE_MATRIX/main.nf`**

Follow the exact structural convention of `modules/binning/BIN_NOVELTY_HITS/main.nf`
(label, `publishDir "${params.outdir}", mode: 'copy'` — per learning L-4, every
module gets one — script block calling the Task 5 CLI by `${projectDir}/bin/...`
absolute path, and a `stub:` block):

```groovy
// Build the genome x cluster incidence matrix + QC/taxonomy join (design
// doc: BUILD_CLUSTER_INCIDENCE_MATRIX). Joins BFD.duckdb's busco_genome/
// asm_stats Parquet tables by ASMID -- see that design doc section for why
// samples.csv alone is not enough (no BUSCO/N50 columns there).

process BUILD_CLUSTER_INCIDENCE_MATRIX {
    label 'build_cluster_incidence_matrix'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(clusters_tsv)
        path(provenance_tsv)
        path(novelty_bins_tsv)
        path(samples_csv)
        path(busco_parquet)
        path(asm_stats_parquet)
        path(genome_classification_tsv)

    output:
        path("incidence_matrix.parquet"), emit: incidence_matrix
        path("cluster_external_status.parquet"), emit: external_status
        path("genome_metadata.parquet"), emit: genome_metadata

    script:
    """
    ${projectDir}/bin/build_cluster_incidence_matrix.py \\
        --clusters-tsv ${clusters_tsv} \\
        --provenance-tsv ${provenance_tsv} \\
        --novelty-bins-tsv ${novelty_bins_tsv} \\
        --samples-csv ${samples_csv} \\
        --busco-parquet ${busco_parquet} \\
        --asm-stats-parquet ${asm_stats_parquet} \\
        --genome-classification-tsv ${genome_classification_tsv} \\
        --min-clade-n ${params.min_clade_n} \\
        --outdir .
    """

    stub:
    """
    touch incidence_matrix.parquet cluster_external_status.parquet genome_metadata.parquet
    """
}
```

- [ ] **Step 2: Write `modules/binning/ACCUMULATION_CURVE/main.nf`**

```groovy
// Resample the incidence matrix into accumulation curves, pan-proteome
// power-law fit, and per-clade Shapley contribution (design doc:
// ACCUMULATION_CURVE). params.n_permutations defaults to 1000, provisional
// pending the permutation-count convergence check (validation #1) -- see
// design doc's "Choosing production P".

process ACCUMULATION_CURVE {
    label 'accumulation_curve'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(incidence_matrix)
        path(external_status)
        path(genome_metadata)
        val(version_tag)

    output:
        path("accumulation_curve_summary.tsv"), emit: summary
        path("pangenome_powerlaw_fit.tsv"), emit: fit
        path("clade_contribution.tsv"), emit: clade_contribution
        path("permutation_marginals.parquet"), emit: raw_marginals

    script:
    """
    ${projectDir}/bin/accumulation_curve.py \\
        --incidence-matrix ${incidence_matrix} \\
        --external-status ${external_status} \\
        --genome-metadata ${genome_metadata} \\
        --n-permutations ${params.n_permutations} \\
        --seed 0 \\
        --version-tag ${version_tag} \\
        --outdir .
    """

    stub:
    """
    touch accumulation_curve_summary.tsv pangenome_powerlaw_fit.tsv \\
          clade_contribution.tsv permutation_marginals.parquet
    """
}
```

- [ ] **Step 3: Add `params.n_permutations` / `params.min_clade_n` defaults to `nextflow.config`**

```groovy
// inside the existing params {} block in nextflow.config
n_permutations = 1000
min_clade_n    = 10
```

- [ ] **Step 4: Wire both processes into `protein_novelty.nf`**

Read the current end of `protein_novelty.nf` first (the file as shown in
this session's earlier exploration ends around `DIAMOND_MAKEDB(uniprot_nr_ch)`
followed by the BFD-side chain through `BIN_NOVELTY_HITS`). Add, after the
existing `BIN_NOVELTY_HITS(...)` call:

```groovy
include { BUILD_CLUSTER_INCIDENCE_MATRIX } from './modules/binning/BUILD_CLUSTER_INCIDENCE_MATRIX/main.nf'
include { ACCUMULATION_CURVE }             from './modules/binning/ACCUMULATION_CURVE/main.nf'
```

(added to the top `include` block alongside the other `modules/binning` includes), and in the `workflow {}` block:

```groovy
    BUILD_CLUSTER_INCIDENCE_MATRIX(
        MMSEQS_LINCLUST_BFD.out.clusters_tsv,
        COMBINE_BFD_PROTEINS.out.provenance_tsv,
        BIN_NOVELTY_HITS.out.bins,
        samples_ch,
        Channel.fromPath(params.bfd_duckdb_busco_parquet, checkIfExists: true),
        Channel.fromPath(params.bfd_duckdb_asm_stats_parquet, checkIfExists: true),
        BUILD_PEP_MANIFEST.out.genome_classification_tsv,
    )

    ACCUMULATION_CURVE(
        BUILD_CLUSTER_INCIDENCE_MATRIX.out.incidence_matrix,
        BUILD_CLUSTER_INCIDENCE_MATRIX.out.external_status,
        BUILD_CLUSTER_INCIDENCE_MATRIX.out.genome_metadata,
        params.bfd_version_tag,
    )
```

Note: the exact upstream channel names (`MMSEQS_LINCLUST_BFD.out.clusters_tsv`,
`COMBINE_BFD_PROTEINS.out.provenance_tsv`, `samples_ch`, and whichever channel
actually carries `genome_classification.tsv` from Stage 1) must be checked
against the real, current `protein_novelty.nf` and `genome_classify.nf` output
emit names at implementation time — read both files in full immediately
before this step, since they may have changed since this plan was written,
and match the exact emit labels rather than guessing. `params.bfd_duckdb_busco_parquet`,
`params.bfd_duckdb_asm_stats_parquet`, and `params.bfd_version_tag` are new
params to add to `nextflow.config` alongside `n_permutations`/`min_clade_n`
(paths/tag for the real run belong in DeltaGain_Fungi's `-c` profile config
per this repo's established split, not hardcoded here — add empty/placeholder
defaults here and real values there).

- [ ] **Step 5: Run the stub pipeline to verify it wires up without error**

Run: `nextflow run protein_novelty.nf -profile stub -stub-run`
Expected: exit 0, both new processes appear in the run log as `COMPLETED`.

- [ ] **Step 6: Commit**

```bash
git add modules/binning/BUILD_CLUSTER_INCIDENCE_MATRIX modules/binning/ACCUMULATION_CURVE \
    protein_novelty.nf nextflow.config
git commit -m "Wire accumulation-curve analysis into protein_novelty.nf

Adds BUILD_CLUSTER_INCIDENCE_MATRIX and ACCUMULATION_CURVE after
MMSEQS_LINCLUST_BFD + BIN_NOVELTY_HITS, per design doc's Pipeline
integration section. Real BFD.duckdb paths and the run's version tag
are supplied via DeltaGain_Fungi's profile config, not hardcoded here.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Deliberately out of scope for this plan

- **Threshold-sensitivity, permutation-convergence, held-out, and
  composition-bias sweeps (figures 5, 6, 9, 10; validations #1, #2, #5, #6)**:
  each requires running `accumulation_curve.py` (and, for #2, upstream
  clustering) multiple times at real scale-up data volumes. Building the
  Nextflow orchestration for these sweeps against synthetic unit-test
  fixtures would produce untested, unrealistic code — do this once the
  500-1000-genome scale-up run (HANDOFF.md item 1) exists to sweep against.
- **`--cap-per-species` / ANI-component capping flag** on `accumulation_curve.py`
  (validation #6): same reasoning — worth adding once there's a real
  composition-bias question to answer at scale-up scale.
- **External-novelty self-inclusion check (validation #7)**: the spec itself
  says to record it only "if obtainable" — no ASMID-to-UniProt-accession
  mapping is known to exist yet; resolve as a research question before
  writing code against it.
- **DeltaGain_Fungi's `profile_protein_novelty.config` additions** (real
  `BFD.duckdb` table paths, SLURM resource directives for the two new
  labels, real `bfd_version_tag`): out of scope for the `DeltaGain` repo
  per the project's established pipeline/data split — a follow-up task in
  that sibling repo once this plan's code lands.

---

## Self-Review

**Spec coverage**: Purpose (n/a, framing only) — Scope/re-run trigger (Task 5
docstring + Global Constraints) — Computational strategy (Tasks 5-9) —
Novelty definition (Task 6) — Data provenance & versioning (Tasks 2-4) —
Components/genome universe edge cases (Task 5) — Clade Shapley (Task 7) —
Power-law fit + CI (Task 8) — Pipeline integration (Task 14) — Figures
1-3 built (Task 12-13); 4-6, 8-10 explicitly deferred with reasoning (see
"Deliberately out of scope") rather than stubbed with fake code — Validation
#3 (Task 11), #4 (Task 10) built; #1, #2, #5, #6, #7 deferred with reasoning
— Storage format (Parquet+zstd used throughout Tasks 5/9; `$SCRATCH`
placement is a Nextflow-process-resource-directive concern that lives in
DeltaGain_Fungi's profile config, out of scope here) — Testing (every task
has unit tests against realistic fixtures; Task 14 adds a stub-run check).

**Placeholder scan**: no TBD/TODO; every code block is complete and runnable
against the given fixtures. Task 14's Step 4 flags real upstream-channel-name
uncertainty explicitly (with an instruction to verify against the live files)
rather than guessing wrong names silently.

**Type consistency**: `asmid_clusters: dict[str, set[str]]` and the record
shape `{position, asmid, internal_cumulative, external_cumulative,
new_internal_count, new_external_count}` from Task 6 are used identically in
Tasks 7, 8, 9, 10, 11. `clade_of: dict[str, tuple[str, str]]` from Task 5's
`genome_metadata` (`clade_rank`, `clade_label` columns) matches Task 7's
`aggregate_clade_contribution` signature. `count_key` string values
(`"new_internal_count"` / `"new_external_count"`) are consistent between
Task 6's record keys and Task 8's `fit_power_law_per_permutation` parameter.
