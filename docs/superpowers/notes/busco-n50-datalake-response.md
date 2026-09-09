# Response: BUSCO/N50 data access under the BFD data-lake redesign

Date: 2026-09-09
Re: `busco-n50-datalake-requirements.md` (2026-09-09)
From: Fungi_BFD — full design at
`Fungi_BFD/docs/superpowers/specs/2026-09-09-bfd-duckdb-datalake-design.md`

Answers below follow your doc's numbering. Short version: your current approach (reading
`tables/busco_genome.parquet` / `tables/asm_stats.parquet` directly via `read_parquet()`) stays
correct and unchanged — nothing in this redesign requires you to change how you open the files.
What changes is what's available *alongside* them, and one gap you should close on your side.

## 1. Discovery mechanism

No glob needed, and none is coming. `tables/<name>.parquet` stays one fixed, unpartitioned file
per logical table — T-014 explicitly rejected bucket-partitioned Parquet (it reintroduces a
small-file metadata cost on our shared storage) and this redesign doesn't reopen that. Your exact
file paths keep working as-is.

## 2. Join key stability

Real gap found on the write side, being fixed: neither `summarize_busco_stats.py` nor
`summarize_asm_stats.py` currently de-duplicates or fails on a repeated `ASMID` before writing
Parquet. We found a concrete path to it — `BUSCO_GENOME`'s output is keyed by
`{ASMID}.BUSCO_summary.{lineage}.txt`; if a genome's BUSCO lineage is ever changed, the old file
isn't cleaned up and both get merged into one row-set. We're adding a hard
`assert_unique_key(df, "ASMID")` (fail loud, exit 1) to both scripts before they write. 0
duplicates exist in the data today, so this is a preventive fix, not a sign anything is
currently wrong.

**On your side**: this makes the "ASMID is unique" contract enforced upstream, but we'd still
recommend keeping your own fail-loud guard rather than removing it once ours lands — defense in
depth, same reasoning as your own L-8 convention, and it's a two-line check regardless of
whether the upstream guarantee is in place.

**Resolution rule if it ever does happen upstream despite the guard**: reject at write time
(pipeline fails, nothing publishes), never latest-write-wins. There is no scenario where a
silently-merged duplicate ASMID should reach a Parquet file you read.

## 3. Schema stability

`sql/schema.sql` (the old hand-maintained doc) was confirmed drifted from what's actually built
and is being retired in favor of `sql/table_schema.yaml` — one machine-readable file per table
(columns, types, key column, per-table `schema_version` int, changelog). A new per-table JSON
manifest (next section) carries that `schema_version` and the current column list, so you can
assert your expected schema version/columns before trusting a read, rather than only discovering
a rename via a downstream `KeyError`.

## 4. Snapshot/freshness identity — this is the main new thing for you

A sidecar JSON file per table, written immediately after (and validated against) its Parquet
file:

```
tables/_manifest/busco_genome.json
tables/_manifest/asm_stats.json
```

```json
{
  "table": "busco_genome",
  "path": "tables/busco_genome.parquet",
  "row_count": 23125,
  "built_at": "2026-09-09T14:22:03Z",
  "built_by": "MERGE_BUSCO_GENOME",
  "merge_run_id": "<nextflow session/task hash>",
  "schema_version": 1,
  "columns": ["ASMID", "complete_pct", "single_pct", "duplicated_pct", "fragmented_pct",
              "missing_pct", "n_markers", "lineage"],
  "parquet_sha256": "<sha256>"
}
```

Read this (plain JSON, no DuckDB dependency needed) and log `built_at`, `row_count`, and
`schema_version` into your own run-provenance manifest alongside the join results — that's
exactly the auditability gap your doc flagged. `merge_run_id` additionally lets you confirm two
tables you're joining (e.g. `busco_genome` + `asm_stats`) came from generations you're comfortable
mixing, if that ever matters for your analysis.

One honest note: we found the live `busco_genome.parquet` is currently stale versus disk right
now (23,125 rows vs. 23,692 actual BUSCO summary files) — this predates the redesign and is
exactly the kind of drift the manifest is meant to surface going forward instead of leaving
silent.

## 5. Update semantics

Confirmed: whole-table rebuild on each merge run, not append-only. A `MERGE_*` process globs
every relevant file on disk and rewrites the entire table's Parquet file from scratch — a row for
a given `ASMID`, once written, is not stable across separate merge runs the way an append-only
log would be. Practically: don't cache/reuse a past read's join result across sessions without
re-checking the manifest's `built_at`/`row_count` first, same as your doc's own caveat already
assumed.

One additional change on our side you don't need to do anything about: today's `MERGE_*` publish
step is not atomic (delete-then-copy), so a read landing mid-rewrite could hit a transient
truncated-file error. We're switching to temp-file + atomic rename (our storage is GPFS, where
this is atomic) — after that lands, a read either sees the complete old file or the complete new
one, never a partial one. This closes a real (if rare) failure window; it doesn't change
anything about how you read.

## Not changing

Nothing about your CLI flags, file paths, or the fact that you read Parquet directly without
opening `db/BFD.duckdb`. If you'd ever rather query via SQL (e.g. to join more tables than just
these two), `db/BFD.duckdb` will still work with the same table names — it's becoming a thin view
catalog over these same Parquet files rather than a copy, transparent to a `SELECT`.
