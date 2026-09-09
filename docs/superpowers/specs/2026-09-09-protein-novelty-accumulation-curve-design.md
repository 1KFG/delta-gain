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
identifying fields with a hash of the full manifest for collision-safety,
e.g. `bfd-v20260909-n11024-mmseqs95c90-busco70n50-a1b2c3d` — date, taxon
count, key clustering params, key QC-threshold params, then a short hash
(e.g. first 7 hex chars of a SHA-256 of the full manifest below) as a tiebreaker
and integrity check. Same scheme, adapted per dataset, for the UniProt
reference side.

**UniProt reference versioning** — `FETCH_UNIPROT_FUNGI_FTP` must record, at
fetch time, either the UniProt release version string (e.g. `2026_02`) if
available from the FTP `reldate.txt` alongside the divisions, or the fetch
date as a fallback when a release string can't be confirmed. The resulting
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
- The short version tag itself

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
  the *expected marginal contribution under random insertion order* — a
  standard, cheap, unbiased way to split "first discovery" credit across
  clades, deliberately simpler than a combinatorially expensive block-Shapley
  computation over clades, which was considered and rejected as not worth
  the added cost here.

- **Heaps' law fit** — nonlinear least-squares fit of
  `new_clusters(N) ~= K * N^(-gamma)` to the mean internal curve and
  separately to the mean external curve, reported with `K`, `gamma`, and R^2.
  Reported as a supporting quantitative summary alongside the descriptive
  curve, not as a replacement for it.

## Pipeline integration

Added to `protein_novelty.nf`, after `MMSEQS_LINCLUST_BFD` + `BIN_NOVELTY_HITS`:

```
MMSEQS_LINCLUST_BFD.out + BUILD_PEP_MANIFEST.out (taxonomy) + BIN_NOVELTY_HITS.out
        |
        v
BUILD_CLUSTER_INCIDENCE_MATRIX  ->  incidence_matrix.tsv, cluster_external_status.tsv
        |
        v
ACCUMULATION_CURVE  (params.n_permutations, params.min_clade_n)
        |
        +-> accumulation_curve_summary.tsv   (per-N mean/CI, internal + external)
        +-> heaps_law_fit.tsv                (K, gamma, R^2 per curve)
        +-> clade_contribution.tsv           (per-clade mean marginal + total contribution, rank used)
```

Plotting is deliberately kept out of the Nextflow DAG: it consumes the three
TSVs above as a standalone step, so regenerating a figure doesn't require
re-running the resampling, and the underlying numbers are always available
as plain tables for review before anyone trusts a plot.

## Outputs / plots

Two consumers of the same three TSVs (plus the validation outputs below), so
interactive and static figures never drift apart:

- **`bin/plot_accumulation_plotly.py`** — interactive HTML.
- **`bin/plot_accumulation_ggplot.R`** — static, publication-styled PDF.

**Bias toward more figures, not fewer.** Each required validation step
(below) should produce its own diagnostic plot, not just a pass/fail number
— both to catch problems visually and because these are exactly the figures
that demonstrate the dataset's dimensions, dynamics, and sensitivity to
parameter/cutoff choices for a paper or report. Minimum figure set, all
generated by both scripts above from the same underlying tables:

1. Internal accumulation curve (mean + permutation-band ribbon), with Heaps'
   law fit line overlaid and gamma/K/R^2 annotated.
2. External-novelty accumulation curve, same treatment.
3. Per-clade marginal-contribution bar chart (at the density-adaptive rank
   chosen per clade), sorted descending.
4. Permutation-count convergence plot (validation #1): gamma estimate and CI
   width vs. P.
5. Clustering-threshold sensitivity plot (validation #2): gamma vs.
   min-seq-id/coverage threshold, small multiples or overlaid curves.
6. QC-confound scatter (validation #3): per-genome marginal new-cluster
   count vs. BUSCO completeness and vs. assembly N50.
7. Positive-control check plot (validation #4): marginal contribution
   distribution for `redundant_strain`-flagged genomes vs. all others.
8. Held-out predictive-check plot (validation #5): first-half-fitted Heaps'
   curve overlaid on the observed second-half trajectory.
9. Composition-bias sensitivity plot (validation #6): curve/gamma with vs.
   without genomes-per-species capping.

Every figure's title/caption should include the run's version tag (see Data
provenance & versioning) so a figure pulled out of context still identifies
which dataset and parameter set produced it.

## Validation (required before trusting a result from this analysis)

All of the following are required steps, not optional sensitivity checks,
before any output from this pipeline is used to support a saturation claim,
sequencing-effort recommendation, or publication figure:

1. **Permutation-count convergence check** — run at increasing P (e.g. 50,
   100, 500, 1000) and confirm the mean curve and gamma estimate stabilize,
   and the CI band narrows roughly as expected with more permutations, before
   picking a production P. Cheap, since resampling doesn't re-run any
   bioinformatics tool.

2. **Clustering-threshold sensitivity** — MMseqs2's min-seq-id/coverage
   thresholds are still placeholders per HANDOFF.md. Recompute the curve and
   gamma at a couple of alternate thresholds at pilot scale. If "open vs.
   closed" flips depending on threshold choice, that is itself a finding to
   report, not a result to average away.

3. **Confound check: annotation/assembly quality vs. genuine novelty** —
   correlate each genome's marginal new-cluster count against that genome's
   BUSCO completeness and assembly N50 (or equivalent QC metrics already
   available from the genome-classification stage) before attributing a
   clade's high marginal contribution to real phylogenetic novelty rather
   than noisier gene prediction (fragmented models, contamination).

4. **Positive-control cross-check with Stage 1** — genomes Stage 1
   (`genome_classify.nf`) already flagged `redundant_strain` must show
   near-zero marginal new-cluster contribution in this analysis. Implement
   as an automated assertion/check in `bin/accumulation_curve.py` output
   validation, not just an eyeballed sanity check — a violation indicates a
   real bug in either the ANI calls or the incidence matrix, not a novel
   finding.

5. **Held-out predictive check** — fit Heaps' law on the first half of a
   permutation's genome-addition order, then measure how well it predicts
   the observed curve on the second half. Report out-of-sample fit quality
   alongside (not instead of) the in-sample R^2, since a model that only
   fits well in-sample is weaker evidence of real saturation behavior.

6. **Composition-bias sensitivity** — the annotated 46% (and, later, the
   full ~22k set) is not a random sample of fungal diversity; some
   species/genera likely have many more annotated strains than others.
   Re-run the curve with genomes-per-species capped (e.g. at 1-2) and check
   whether gamma changes meaningfully — this tells you whether the
   "openness" estimate is being driven by real diversity or by oversampled
   strains, and must be reported alongside the uncapped result.

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
  Plain gzip or zstd remains appropriate for the smaller flat TSV/CSV summary
  files (`accumulation_curve_summary.tsv`, `heaps_law_fit.tsv`,
  `clade_contribution.tsv`) kept around for direct human inspection.
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
  order) — verifies the resampling logic and Heaps' fit math without needing
  real bioinformatics tools.
- Smoke test against the existing n=100 pilot data
  (`bfd_proteins_clusters.tsv` etc. from job 28206550) before trusting output
  at the 500-1000-genome scale-up — confirms the incidence-matrix builder
  correctly parses real cluster/provenance/taxonomy files, though N is too
  small for the Heaps' fit itself to be meaningful yet.
