# Protein-content accumulation curves for BFD fungal novelty screening

Status: approved design, not yet implemented
Date: 2026-09-09
Repo: DeltaGain (pipeline/code). Applies to data staged in DeltaGain_Fungi.

## Purpose

Quantify how much new protein content is revealed as more fungal genomes are
sequenced and annotated into BFD, in a way that is:

- A **saturation/completeness claim**: is discovery of new protein content
  still rising steeply, or flattening (open vs. closed pan-proteome, in
  pan-genomics terms), across the currently-annotated subset?
- **Sequencing-effort guidance**: which taxonomic clades are still
  contributing disproportionately much new content per genome added, as a
  signal for where more sequencing/annotation effort would pay off?
- A **pipeline/QC diagnostic**: does the novelty-screen behave sensibly as N
  grows, or does something (thresholds, reference-set gaps, annotation
  artifacts) look wrong?
- A **publication figure**: presentable, defensible quantification, not just
  an internal diagnostic.
- A demonstration of the **BFD annotation project's own value**: how much
  protein-content information was invisible before genomes had gene
  predictions, made visible now that annotation coverage is growing.

## Scope and re-run trigger

As of 2026-09-08, only 10,922 / 23,683 (46%) of `samples.csv` rows have a
matching pep file (`Fungi_BFD_runs/input/pep/`) — see HANDOFF.md's annotation
coverage note. This design targets the next scale-up run (HANDOFF item 1,
500–1000 genomes), reusing its clustering output. **This analysis will need
a full re-run — rebuilding the incidence matrix, not redesigning anything —
once BFD annotation coverage reaches roughly the full ~22k genome set.**
The design itself does not need to change for that re-run; only the input
scale does.

## Computational strategy

Cluster BFD's proteins **once** at the target scale using MMseqs2 linclust
(staying with MMseqs2 per project decision on 2026-09-08 — see the
DIAMOND-cluster vs. MMseqs2-linclust runtime/sensitivity comparison recorded
earlier in this project's history). From that single clustering run, build a
**genome × cluster incidence matrix**, plus each cluster's existing external
reference-hit status from `bin_novelty_hits.py`'s output (strong/weak/no-hit
vs. UniProt-fungi + the existing nr reference).

Every curve, permutation, and per-clade statistic below is computed by
**resampling over this precomputed matrix** — set arithmetic over
already-known cluster memberships, not repeated DIAMOND/MMseqs2 runs. This
makes permutation count (hundreds to thousands) computationally cheap
regardless of how many times the analysis is rerun, and means a future
re-derivation (e.g. at full ~22k-genome scale) only requires rebuilding the
matrix once from that run's clustering output, not repeating this whole
analysis design.

## Novelty definition: two curves, not one

- **Internal accumulation** — a protein counts as new if it doesn't cluster
  with any protein already accumulated from prior genomes in the
  permutation order. This is BFD's own pan-proteome growth, independent of
  external databases (classic pan-genome/Heaps'-law framing).
- **External novelty accumulation** — a protein counts as new only if it is
  a no-hit against UniProt-fungi + the existing nr reference (today's Stage
  3 novelty-bin definition), tracked as it accumulates across genomes.

Both curves are tracked side by side on the same genome-addition axis, from
the same incidence matrix.

## Data provenance & versioning

This closes a real gap: `FETCH_UNIPROT_FUNGI_FTP` currently pulls from
UniProt's `current_release/` FTP path with no release version or fetch date
captured anywhere in its output, and `DIAMOND_CLUSTER`'s container version is
pinned in a `params` value but never written into any output manifest.
Because BFD clustering will likely be re-run multiple times with different
parameter sets and QC filters, both reference and BFD-side datasets need an
explicit, unambiguous version tag from here on.

**Version tag format**: a short hash-based tag combining human-readable
identifying fields with a hash for collision-safety, e.g.
`bfd-v20260909-n11024-mmseqs95c90-busco70n50-a1b2c3d` — date, taxon count,
key clustering params, key QC-threshold params, then a short hash (first 7
hex chars of a SHA-256) as a tiebreaker and integrity check. The hash **must
be computed over the manifest with the `version_tag` field itself excluded**
(hash-then-embed, not hash-of-something-containing-itself) to avoid a
circular definition. Same scheme, adapted per dataset, for the UniProt
reference side.

**UniProt reference versioning** — `FETCH_UNIPROT_FUNGI_FTP` must download
`https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt` into the same
output folder as the fetched divisions, at fetch time, and parse the current
release version and date from it (confirmed real content as of 2026-09-09:
release `2026_03`, made 02-Sept-2026). This is the authoritative release
identifier — use it as the release version field rather than the fetch date,
falling back to the fetch date only if `relnotes.txt` is unreachable or its
format can't be parsed. The resulting
`uniprot_fungi_nr.fasta` (and any DIAMOND db built from it) gets a manifest
sidecar (e.g. `uniprot_fungi_nr.manifest.json`) recording:
- UniProt release version and/or fetch date
- MMseqs2 version and the exact linclust parameters used to derive the nr set
  (`--min-seq-id`, `-c`, `--cov-mode`)
- Sequence count in and out of dereplication
- The short version tag itself

**BFD dataset versioning** — each BFD clustering run gets a manifest
(`bfd_cluster_manifest.json`) recording:
- Date of the run
- Number of taxa/genomes included (and the annotation-coverage fraction at
  that time, e.g. "10,922 / 23,683 annotated")
- Clustering tool + version + full parameter set (MMseqs2 linclust
  `--min-seq-id`/`-c`/`--cov-mode`, or DIAMOND cluster's equivalent if that
  path is ever revisited)
- QC filter thresholds applied (BUSCO completeness cutoff, assembly N50
  cutoff, BUSCO-of-predicted-proteins cutoff if used) and, explicitly, the
  **list or count of genomes screened out by each filter** — not just the
  ones kept, since which genomes got excluded and why matters for
  reproducing or auditing a specific run
- Genomes excluded because they had zero clustered proteins after QC
  filtering (see Genome universe below) — recorded separately from
  QC-threshold screen-outs, since "excluded by threshold" and "produced
  nothing after filtering" are different failure modes worth telling apart
- The short version tag itself

**External reference identity** — this pipeline screens against two
distinct, separately-versioned external references, and both must be
recorded per run, not just "the nr reference": (1) NCBI nr, specifically
`nr_cluster_seq.dmnd` at `/srv/projects/db/ncbi/diamond/20260128/`
(`modules/search/DIAMOND_NOVELTY_SCREEN/main.nf`) — record its path's
existing datestamp (`20260128`) in the manifest; (2) the in-house
dereplicated UniProt-fungi set (`uniprot_fungi_nr.fasta`), versioned per the
UniProt reference versioning above. "External novelty" in this design always
means no-hit against **both**.

**BUSCO/N50 source (resolved 2026-09-09)**: no BUSCO completeness score or
assembly N50 column exists in `samples.csv` (only `BUSCO_LINEAGE`, the
lineage *name*, not a percentage) — but both metrics already exist,
ASMID-keyed, in `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/db/BFD.duckdb`,
backed by Parquet tables under `Fungi_BFD_runs/tables/`: `asm_stats.parquet`
(`N50_bp`, plus `contig_count`, `total_length_bp`, `gc_pct`, and other
assembly-QC fields) and `busco_genome.parquet` (`complete_pct`,
`single_pct`, `duplicated_pct`, `fragmented_pct`, `missing_pct`, `lineage`).
`BUILD_CLUSTER_INCIDENCE_MATRIX` should join against these two tables
(by `ASMID`) to attach QC metrics to the incidence matrix, rather than
against `samples.csv` alone.
**Caveat**: this DuckDB was last built 2026-08-05 and will need rebuilding
against BFD's latest genome-processing state before a real run at
scale-up or full ~22k scale — verify its ASMID coverage against the
annotation set actually being used each time this analysis runs, don't
assume it's current.

All downstream outputs from a given run (incidence matrix, accumulation
curves, clade contribution tables, figures) should be named or namespaced
with that run's version tag, so which dataset/parameter-set produced a given
plot is never ambiguous once several iterations exist side by side.

## Components

- **`bin/build_cluster_incidence_matrix.py`** (new). Inputs:
  `bfd_proteins_clusters.tsv` (cluster ID <-> protein ID),
  `bfd_proteins_provenance.tsv` (protein ID <-> ASMID),
  `bfd_novelty_bins.tsv` (cluster representative's external-hit status),
  `samples.csv` (ASMID -> taxonomy ranks: PHYLUM/SUBPHYLUM/CLASS/SUBCLASS/
  ORDER/FAMILY/GENUS/SPECIES). Output: one row per (genome, cluster) pair
  present, plus a per-cluster external-status lookup table.

**Genome universe and edge cases** — the permutation set for
`accumulation_curve.py` must be seeded from the ASMIDs actually present in
`bfd_proteins_provenance.tsv` (genomes that produced at least one clustered
protein), never from `samples.csv` directly — this repeats, at the
accumulation-curve layer, the exact false-singleton bug already found and
fixed once in `find_components.py` (a genome absent from the real input set
must never silently appear as if it contributed zero, or worse, get treated
as present with undefined behavior). Explicitly:
- Genomes with zero clustered proteins after QC filtering are excluded from
  the permutation universe; their count and ASMIDs are recorded in the run
  manifest as a distinct exclusion reason (see Data provenance above).
- Duplicate ASMIDs in the input are a data-integrity error, not something to
  silently deduplicate or average over — `build_cluster_incidence_matrix.py`
  must assert ASMID uniqueness and fail loudly if violated.
- A genome with a missing/blank value at the density-adaptive taxonomic rank
  rolls up to the next coarser rank automatically (same rollup logic used
  for low-N ranks); if unresolvable all the way to phylum, it is bucketed as
  `unclassified` rather than dropped from the per-clade analysis.

- **`bin/accumulation_curve.py`** (new). The resampling engine. For P random
  permutations of the annotated-genome set:
  - Walk the permutation; at each step record cumulative distinct clusters
    seen so far (internal curve) and cumulative distinct clusters seen so
    far whose external status is no-hit (external-novelty curve).
  - Record, per genome in the walk, how many *new* clusters (internal) and
    new *no-hit* clusters (external) that specific genome introduced — this
    per-genome marginal count feeds the per-clade view at zero extra
    permutation cost.

- **Clade attribution method** — for each genome's recorded marginal
  contribution, tag it with that genome's clade at the finest taxonomic rank
  meeting a density threshold (e.g. family if the family has >= `min_clade_n`
  annotated genomes, else roll up to order, then class, then phylum).
  Aggregate: mean per-clade marginal contribution per genome added, and mean
  per-clade total contribution, averaged across all P permutations. This is
  the *expected marginal contribution under random insertion order* — this
  is exactly the **genome-level Shapley value** for the "new clusters"
  characteristic function, computed efficiently via random-permutation
  averaging rather than exact combinatorics (the standard practical
  approximation to Shapley values). Clade-level labels built by rolling
  low-N genomes up to a coarser rank must be shown as e.g. `ORDER (other
  families)`, not as the bare order name, so a bar in the per-clade plot is
  never misread as "the whole order" when it actually represents only the
  families too sparse to report individually.

- **Pan-proteome power-law fit** (Tettelin et al. 2005 form; earlier drafts
  of this spec called this "Heaps' law fit" and conflated the cumulative and
  rate forms — corrected here). Fit
  `n_new(N) ~= kappa * N^(-alpha)` — where `n_new(N)` is the **mean marginal
  (per-genome, per-permutation-position) new-cluster count**, not the
  cumulative curve — via nonlinear least squares, separately for the
  internal and external series. Decision rule: `alpha <= 1` indicates an
  open pan-proteome (still discovering at a sustained rate); `alpha > 1`
  indicates it is closing/closed. Fitting this rate form directly to the
  *cumulative* curve was considered and rejected: cumulative curves are
  monotonic and smooth by construction, so R^2 there is close to 1
  regardless of whether the underlying process is open or closed and is not
  informative. The cumulative accumulation curve (figures 1-2) remains the
  primary descriptive visualization; the fitted rate-form model may be
  integrated and overlaid on it as a predicted trajectory for visual
  comparison, but `kappa`/`alpha`/R^2 are reported from the rate-curve fit,
  not claimed from the cumulative overlay.
  **Uncertainty**: fit this model per permutation (not once to the
  across-permutation mean), giving a distribution of `alpha` across the P
  permutations. Report the mean `alpha` with its 2.5/97.5 percentile
  interval (see CI-band definition below) rather than a single point
  estimate with no uncertainty.

## Pipeline integration

Added to `protein_novelty.nf`, after `MMSEQS_LINCLUST_BFD` + `BIN_NOVELTY_HITS`:

```
MMSEQS_LINCLUST_BFD.out + BUILD_PEP_MANIFEST.out (taxonomy) + BIN_NOVELTY_HITS.out
        |
        v
BUILD_CLUSTER_INCIDENCE_MATRIX  ->  incidence_matrix.parquet, cluster_external_status.parquet
        |
        v
ACCUMULATION_CURVE  (params.n_permutations, params.min_clade_n)
        |
        +-> accumulation_curve_summary.tsv   (per-N mean/CI, internal + external)
        +-> pangenome_powerlaw_fit.tsv        (kappa, alpha + 95% CI, R^2 per curve)
        +-> clade_contribution.tsv           (per-clade mean marginal + total contribution, rank used)
```

Plotting is deliberately kept out of the Nextflow DAG: it consumes the three
TSVs above as a standalone step, so regenerating a figure doesn't require
re-running the resampling, and the underlying numbers are always available
as plain tables for review before anyone trusts a plot.

## Outputs / plots

**`bin/plot_accumulation_ggplot.R`** is the primary implementation, covering
the full figure set below as static, publication-styled PDFs. A second,
lightweight **`bin/plot_accumulation_plotly.py`** covers only figures 1-3 —
the ones people actually want to explore interactively (hover per-N values,
zoom, drill into a clade bar) — rather than duplicating all figures in both
stacks; the diagnostic/validation figures (4-10 below) are one-shot checks
read once per run, where a static PDF is sufficient and keeping only one
implementation to maintain is the better trade for a research pipeline.

**Bias toward more figures, not fewer.** Each required validation step
should produce its own diagnostic plot, not just a pass/fail number — both
to catch problems visually and because these are exactly the figures that
demonstrate the dataset's dimensions, dynamics, and sensitivity to
parameter/cutoff choices for a paper or report. Minimum figure set:

1. *(ggplot + plotly)* Internal accumulation curve (mean + permutation-band
   ribbon), with the integrated pan-proteome power-law fit overlaid as a
   predicted trajectory, `alpha`/`kappa` annotated.
2. *(ggplot + plotly)* External-novelty accumulation curve, same treatment.
3. *(ggplot + plotly)* Per-clade marginal-contribution bar chart (at the
   density-adaptive rank chosen per clade, roll-up bars labeled e.g. `ORDER
   (other families)`), sorted descending.
4. *(ggplot)* Pan-proteome power-law rate-curve fit plot: mean marginal
   new-cluster count vs. N with the fitted `kappa * N^(-alpha)` curve and R^2
   annotated — this is where the fit quality is actually shown, since it is
   not informative on the cumulative curve (see the fit definition above).
5. *(ggplot)* Permutation-count convergence plot (validation #1): `alpha`
   estimate and CI width vs. P.
6. *(ggplot)* Clustering-threshold sensitivity plot (validation #2, at
   scale-up scale): `alpha` vs. min-seq-id/coverage threshold.
7. *(ggplot)* QC-confound scatter (validation #3): per-genome marginal
   new-cluster count (raw and per-1,000-proteins-normalized) vs. proteome
   size, vs. BUSCO `complete_pct`, and vs. `N50_bp`.
8. *(ggplot)* Positive-control check plot (validation #4): conditioned
   marginal-contribution distribution for `redundant_strain`-flagged genomes
   vs. the unconditioned distribution for `singleton_isolated`/
   `novel_strain_same_species` genomes.
9. *(ggplot)* Held-out predictive-check plot (validation #5): distribution
   of out-of-sample fit quality across permutations, first-half-fitted curve
   overlaid on an example second-half trajectory.
10. *(ggplot)* Composition-bias sensitivity plot (validation #6): curve/
    `alpha` with vs. without genomes-per-species (or ANI-component) capping.

Every figure's title/caption should include the run's version tag (see Data
provenance & versioning) so a figure pulled out of context still identifies
which dataset and parameter set produced it.

## Validation (required before trusting a result from this analysis)

All of the following are required steps, not optional sensitivity checks,
before any output from this pipeline is used to support a saturation claim,
sequencing-effort recommendation, or publication figure:

1. **Permutation-count convergence check** — run at increasing P (e.g. 50,
   100, 500, 1000) and confirm the mean curve and the per-permutation
   `alpha` distribution (mean and 2.5/97.5 percentile width, not just the
   mean) stabilize as P grows, before picking a production P. Cheap, since
   resampling doesn't re-run any bioinformatics tool. This is also the
   source of the CI reported in the pan-proteome power-law fit above — the
   fit's uncertainty comes from this per-permutation distribution, not from
   a single fit to the mean curve.

2. **Clustering-threshold sensitivity** — MMseqs2's min-seq-id/coverage
   thresholds are still placeholders per HANDOFF.md. Recompute the curve and
   `alpha` at a couple of alternate thresholds **at the scale-up run
   (500-1000 genomes or larger), not at the n=100 pilot** — the Testing
   section below already notes n=100 is too small for the fit itself to be
   meaningful, so running this check at that scale would just add noise to
   noise. If "open vs. closed" flips depending on threshold choice, that is
   itself a finding to report, not a result to average away.

3. **Confound check: annotation/assembly quality vs. genuine novelty** —
   correlate each genome's marginal new-cluster count against (a) that
   genome's BUSCO `complete_pct` and `N50_bp` (from `BFD.duckdb` /
   `busco_genome.parquet` + `asm_stats.parquet`, joined by ASMID — see Data
   provenance & versioning above; verify the DuckDB's ASMID coverage matches
   the genome set actually in use, since it will need rebuilding before
   scale-up) and (b) that genome's total predicted-protein count (proteome size),
   which *is* already available and is an obvious, currently-missing
   confound: a genome with more predicted proteins has a higher chance of
   containing something novel by sheer count, independent of real
   phylogenetic novelty. Report marginal new-cluster count normalized per
   1,000 proteins alongside the raw count for this reason. Do not attribute
   a clade's high marginal contribution to real phylogenetic novelty until
   both confounds have been checked.

4. **Positive-control cross-check with Stage 1** — naively checking that
   `redundant_strain`-flagged genomes have near-zero marginal contribution
   *averaged over all random permutations* is wrong: in roughly half of
   permutations such a genome is inserted *before* its ANI partner and
   receives full first-discovery credit for their shared content, so its
   unconditional average will not be near zero even when everything is
   correct. The check must instead be **conditioned**: for each
   `redundant_strain` genome, compute its marginal contribution only across
   permutations where it is inserted *after* at least one other member of
   its same ANI component (from `genome_classify.nf`'s component/cluster
   assignment, not just the `redundant_strain` label alone). Define "near
   zero" concretely as a configurable threshold, e.g. median
   conditioned-marginal-contribution `< 5%` of that genome's own total
   distinct protein count (`redundant_strain_max_marginal_frac`, default
   0.05) — and compare that conditioned distribution against the
   unconditioned distribution for `singleton_isolated`/`novel_strain_same_species`
   genomes as a sanity contrast. Implement as an automated assertion in
   `bin/accumulation_curve.py`'s output validation; a violation indicates a
   real bug in the ANI calls or the incidence matrix, not a novel finding.

5. **Held-out predictive check** — fit the pan-proteome power-law model on
   the first half of a permutation's genome-addition order, then measure how
   well it predicts the observed curve on the second half, **repeated across
   permutations** (not a single permutation), reporting the distribution of
   out-of-sample fit quality. Report this out-of-sample distribution
   alongside (not instead of) in-sample R^2, since a model that only fits
   well in-sample is weaker evidence of real saturation behavior than one
   with genuine out-of-sample predictive accuracy across many permutations.

6. **Composition-bias sensitivity** — the annotated 46% (and, later, the
   full ~22k set) is not a random sample of fungal diversity; some
   species/genera likely have many more annotated strains than others.
   Re-run the curve with genomes-per-species capped (e.g. at 1-2) and check
   whether `alpha` changes meaningfully. The capped subset must be **redrawn
   independently within each permutation** (not fixed once up front), so the
   check reflects the sensitivity of the random-order process itself rather
   than one arbitrary choice of which strains to keep. Consider using Stage
   1's ANI-based components (`genome_classify.nf`) as the capping unit
   instead of the `SPECIES` string field, since ANI components are a more
   direct measure of redundancy than a taxonomy label. Report alongside the
   uncapped result, not as a replacement for it.

7. **External-novelty self-inclusion confound** — a BFD genome's proteins
   may already be present in UniProt/TrEMBL (e.g. if that genome or a close
   relative was previously deposited and annotated by UniProt), in which
   case "no-hit against UniProt-fungi" partly reflects sequence-matching
   noise rather than genuine absence from the public record. If per-ASMID
   UniProt-presence/accession-mapping is obtainable, record it and treat
   such genomes as a separate stratum in the external-novelty curve.
   Otherwise, this caveat must be stated explicitly wherever the external
   curve is used to support the "value of the annotation project" framing
   from the Purpose section.

## Known caveat (record, do not attempt to solve here)

The annotated subset is not a random sample of fungal diversity — it
reflects whatever has been sequenced/annotated so far. The accumulation
curve therefore describes discovery-within-the-annotated-subset, not a
calibrated estimate of "true fungal pan-proteome size." This caveat must be
stated explicitly wherever any result from this analysis is reported
(figure caption, report text, etc.), so the saturation claim is not
overstated.

## Storage format & compute-storage policy

- **File format**: large tabular outputs (incidence matrix, per-permutation
  raw results before summarization, cluster tables) should be written as
  **Parquet with zstd compression** rather than plain TSV — Parquet already
  compresses at the column level, so layering a separate gzip pass on top of
  Parquet is redundant and can hurt (compressing already-compressed data).
  The three small flat summary TSVs (`accumulation_curve_summary.tsv`,
  `pangenome_powerlaw_fit.tsv`, `clade_contribution.tsv`) stay **plain,
  uncompressed TSV** — they are kept specifically for direct human
  inspection and diffing between runs, and compressing them buys negligible
  space at the cost of that convenience.
- **Compute-storage placement**: per this project's established HPCC
  convention (`$SCRATCH` is node-local NVMe, not shared), large intermediate
  data during a SLURM job — the raw permutation walk output before it's
  summarized, any uncompressed intermediate matrix — should be written to
  `$SCRATCH` and compressed there before copying the final Parquet/TSV
  outputs back to `/bigdata`. This isn't an absolute rule: small-to-medium
  files (the three summary TSVs, manifests, figures) are fine written and
  compressed directly on `/bigdata`; it's specifically the large cluster/
  incidence-matrix-scale files where scratch-then-compress-then-copy matters.

## Testing

- Unit-testable on a small synthetic incidence matrix (a handful of genomes,
  a handful of clusters, known expected new-cluster counts per permutation
  order) — verifies the resampling logic and pan-proteome power-law fit math
  without needing real bioinformatics tools.
- Smoke test against the existing n=100 pilot data
  (`bfd_proteins_clusters.tsv` etc. from job 28206550) before trusting output
  at the 500-1000-genome scale-up — confirms the incidence-matrix builder
  correctly parses real cluster/provenance/taxonomy files, though N is too
  small for the power-law fit itself to be meaningful yet.
