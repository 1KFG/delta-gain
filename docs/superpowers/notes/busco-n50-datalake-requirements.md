# BUSCO/N50 data access: requirements for a data-lake redesign

Date: 2026-09-09
From: DeltaGain (protein-content accumulation-curve analysis)
Re: proposal to treat `Fungi_BFD_runs/tables/` as data-lake access
rather than reads through a single monolithic `BFD.duckdb` file

## Context

DeltaGain's accumulation-curve analysis (`bin/build_cluster_incidence_matrix.py`)
needs per-genome BUSCO completeness and assembly N50 to run a required
confound check (does a clade's apparent novelty come from real biology,
or from noisier/lower-quality genome assembly and annotation?).

## Current state (verified 2026-09-09)

- `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/busco_genome.parquet`
  and `.../asm_stats.parquet` are each a single flat Parquet file.
- Both currently hold 23,125 rows, dated 2026-08-04, written by DuckDB
  v1.1.3 — one day before `Fungi_BFD_runs/db/BFD.duckdb` itself
  (2026-08-05), i.e. these Parquet files are themselves an export
  *from* that DuckDB build, not an independent primary source yet.
- `busco_genome.parquet` columns: `ASMID, complete_pct, single_pct,
  duplicated_pct, fragmented_pct, missing_pct, n_markers, lineage`.
- `asm_stats.parquet` columns: `ASMID, SPECIES, STRAIN, contig_count,
  total_length_bp, min_contig_bp, max_contig_bp, median_contig_bp,
  mean_contig_bp, L50, N50_bp, L90, N90_bp, gc_pct, n_gap_count,
  total_n_bases, masked_bases, masked_pct, t2t_scaffolds,
  telomere_fwd, telomere_rev`.
- DeltaGain already reads these two files **directly** via DuckDB's
  `read_parquet('<path>')` — it never opens `BFD.duckdb` itself. At
  the code level this is already Parquet-native, not
  monolithic-database-native. The open question below is what changes
  if the *storage layout* becomes a real multi-file/partitioned lake
  instead of one flat file per table.

## What DeltaGain needs, to stay correct under a data-lake redesign

1. **Discovery mechanism.** If a logical table (e.g. `busco_genome`)
   becomes multiple Parquet files — partitioned by batch, date, or
   species — DeltaGain needs either a stable glob pattern (DuckDB's
   `read_parquet()` accepts globs natively, e.g.
   `tables/busco_genome/*.parquet`) or a small manifest file listing
   the current file paths per logical table. Right now
   `build_cluster_incidence_matrix.py` takes one exact file path per
   `--busco-parquet`/`--asm-stats-parquet` CLI flag; this would need
   to become a glob-aware path or a manifest read.

2. **Join key stability.** `ASMID` must remain the unique join key
   across every file/partition. DeltaGain's current code builds a
   plain Python dict from these two tables with **no duplicate-ASMID
   guard** — unlike every other file this pipeline reads, which fails
   loudly on a duplicate key (a project convention established after
   a real false-singleton bug, see `.living/learnings.md` L-8). If
   incremental/appended writes could ever produce a duplicate ASMID
   across files or partitions, DeltaGain needs to know the intended
   resolution rule (latest-write-wins? reject at write time?) so a
   matching fail-loud or dedup guard can be added on the read side.

3. **Schema stability.** Column names/types above are hardcoded in
   DeltaGain's reader. Any rename, type change, or removed column
   needs to be a documented, versioned schema change — ideally with
   enough lead time to update the reader — not a silent change that
   only surfaces as a downstream `KeyError`.

4. **Snapshot/freshness identity.** DeltaGain's own design doc already
   flags a real caveat: a BUSCO/N50 snapshot may not match the genome
   set a given accumulation-curve run actually uses (annotation
   coverage keeps growing). A data-lake layout should expose something
   DeltaGain can record in its own run provenance manifest — a build
   timestamp, a row count, or an explicit version string per
   snapshot — so "which snapshot produced this result" is auditable
   from the output alone, not inferred from file mtimes the way it is
   today.

5. **Update semantics.** Is a partition/file append-only (a row, once
   written for a given ASMID, never changes) or is the whole table
   rebuilt from scratch on each refresh (today's apparent behavior,
   since it's exported from a rebuilt DuckDB)? This determines whether
   DeltaGain can treat a past run's join results as still valid without
   re-pulling, or must always re-read fresh before trusting them.

## Not required, just useful context

DeltaGain's join is currently read-only and infrequent (once per
accumulation-curve run, not once per genome or per task) — it does not
need low-latency point lookups, streaming reads, or write access. A
solution optimized for occasional bulk reads of two tables by a stable
join key is sufficient; nothing here requires a general-purpose query
engine beyond what DuckDB's `read_parquet()` already provides.
