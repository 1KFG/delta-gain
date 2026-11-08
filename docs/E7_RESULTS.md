# E7 — lineage-specific expansions and secretome (SSP) novelty

Run 2026-11-08 (job 31865426; scripts `bin/e7_stage1.py`, `bin/e7_cys.py`, `bin/e7_analyze.py`, launcher `launch/e7.sbatch`;
results `DeltaGain_Fungi/results/e7/`: `e7_summary.txt`, `e7_ssp_within_order.tsv`, `e7_survivor_ssp.tsv`,
`e7_family_copy_number.tsv`, `e7_genome_counts.tsv`).

## Method

* **Units.** All statistics are species-level: genomes of one species in one class (baseline / new) are averaged into one unit, so
  strains do not inflate n. New units are split into `N_known` (species also has a baseline genome, 833 units) and `N_newsp`
  (new species, 4,540 units); 2,682 baseline units. Only genomes of class baseline/new enter (34 genomes in other classes are excluded).
* **Within-order contrast.** diff = sum_o w_o (mean_new_o − mean_base_o) / sum_o w_o, w_o = n_b·n_n/(n_b+n_n), over orders with
  ≥ 3 units of each class (71 orders for SSP, 71 for copy number). SSP/SP: permutation of class labels inside orders
  (3,000 draws) and bootstrap CI over units (1,000). Copy number: bootstrap over units inside (order, class) (1,000), two-sided
  p, Benjamini–Hochberg q over 323 families.
* **Signal peptide (SP).** SignalP probability ≥ 0.5 (best row per protein).
* **SSP candidate.** SP + length ≤ 300 aa + ≥ 4 Cys (counted from the BFD protein FASTAs) + no GPI-anchor (PredGPI
  `GPI-anchor` feature, 3.3 M proteins) + no TM helix beyond a signal anchor (TMHMM PredHel = 0, or PredHel = 1 whose helix
  starts ≤ SP end + 10; proteins absent from the TMHMM table are treated as 0 helices — the table contains no PredHel = 0 rows).
  It is a heuristic, **not** a functional effector call.
* **Copy number.** Per genome, number of distinct proteins with a dbCAN HMM hit of the family (`HMM_id` without `_N` subfamily
  suffix) or a MEROPS hit whose MEROPS id maps to the family (`merops_lib.families.tab`); zeros kept. Families seen in ≥ 50
  genomes (219 CAZy, 104 MEROPS). Two scales: raw copies per genome, and copies per 1,000 proteins (controls for proteome size;
  new proteomes average 10.2 k proteins vs 11.1 k for baseline).
* **Survivors.** Survivor representatives (`survivor_annot.tsv`, novel at 50 %) by genome-count class; expected SP/SSP counts use the same
  genome's overall rate (observed/expected), which removes clade and genome-quality composition.

## Results

### Signal peptides and SSP candidates, new vs baseline (within order, species level)

Whole proteomes: SP fraction 7.64 % (new) vs 7.71 % (baseline); SSP candidates 1.369 % vs 1.359 % (2.51 M and 0.84 M proteins).

| metric | contrast | within-order diff (95 % CI) | perm p |
|---|---|---|---|
| SP fraction | new − baseline | −0.144 pts (−0.211, −0.077) | 0.0003 |
| SSP fraction | new − baseline | −0.037 pts (−0.056, −0.018) | 0.0003 |
| SSP among SP proteins | new − baseline | −0.17 pts (−0.33, −0.01) | 0.034 |
| SSP fraction | known-species new − baseline | +0.001 pts (−0.026, +0.027) | 0.94 |
| SSP fraction | new-species − baseline | −0.045 pts (−0.064, −0.026) | 0.0003 |

New genomes do **not** carry more SSP candidates than same-order baseline genomes; the effect is small (about −3 % relative) and comes
from new species, not from new strains of known species.

### Survivors (novel at 50 %) are not enriched for secreted or small-secreted proteins

| genome-count class | reps | SP frac | SP obs/exp | SSP frac | SSP obs/exp |
|---|---|---|---|---|---|
| 1 | 12.23 M | 5.1 % | 0.82 | 0.93 % | 0.88 |
| 2–4 | 3.60 M | 6.0 % | 0.89 | 1.09 % | 0.98 |
| 5–19 | 1.07 M | 6.1 % | 0.90 | 1.13 % | 1.01 |
| ≥ 20 | 0.18 M | 3.9 % | 0.72 | 0.84 % | 0.79 |

(all new proteins: SP 7.64 %, SSP 1.37 %). Survivor representatives are, if anything, depleted in signal peptides relative to their own
genomes; SSP candidates occur at about the genome rate (obs/exp 0.98–1.01) for multi-genome clusters and below it for singletons and
very widespread clusters. The novel fraction is therefore **not** driven by a secretome/effector-like class.

### CAZy / MEROPS copy number, new vs baseline within order (71 orders, 8,055 units, 323 families)

| scale | tested | expanded (q < 0.05) | contracted (q < 0.05) |
|---|---|---|---|
| raw copies per genome | 323 (219 CAZy, 104 MEROPS) | 12 | 174 |
| copies per 1,000 proteins | 323 | 19 (12 CAZy, 7 MEROPS) | 23 (18 CAZy, 5 MEROPS) |

The raw scale is dominated by proteome size (new genomes have about 8 % fewer proteins, many are drafts), which makes about half of all
families look contracted; the per-1,000 scale is the fair comparison and gives a balanced, small result. Absolute effects are tiny
(≤ 0.014 copies per 1,000 proteins); relative effects are 4–15 % for common families (GH33, GH36, M14, GT34) and up to several-fold
for rare ones (GH19, CE14, S41). Per-1,000 expansions: CE14, GH19, S41, A39, GH42, GH33, GH36, GH188, AA5, M14, GT2_Glyco_trans_2, GT34.
Contractions: MEROPS C100 and M100, PL20, PL38, GH30, GH135, GT25, CE17 (raw scale additionally CBM13, PL35, CBM52, GH44, I66).
Full table: `e7_family_copy_number.tsv`.

## Reading

* No evidence that new genomes add a secretome or effector-like (SSP) novelty beyond same-order baseline genomes, and the survivor pool is not enriched for signal peptides.
* Lineage-specific CAZy/MEROPS expansion is not a general feature of new genomes: after proteome-size normalisation 19 families expand
  and 23 contract, with small effects. This matches the family-level nulls (new genomes add about 1–3 points beyond same-order baseline genomes).

## Limits (honest)

* SSP calling is a standard heuristic with fixed thresholds (SP ≥ 0.5, ≤ 300 aa, ≥ 4 Cys); no sensitivity analysis over thresholds, no
  effector-prediction tool, and TMHMM absence is read as "no helix".
* BFD annotated both groups with the same pipeline, but new genomes include more drafts and metagenome-derived assemblies, so
  proteome size and gene-model quality differ; per-1,000 normalisation does not remove incompleteness of individual families.
* Copy number counts proteins with a hit, not independent duplications (no phylogeny, no synteny); "expanded" means a higher mean count only.
* Only the 71 orders with ≥ 3 units in each class contribute; clades without baseline representatives are not tested.
* Survivor SP/SSP labels refer to the cluster representative only.
* Family-level bootstrap treats species within an order as exchangeable and does not account for phylogenetic relatedness inside orders.
