# Family-level diversity plan (Pfam · CAZy · MEROPS)

Status: **draft v1, 2026-11-06** · Companion to [`DELTA_GAIN_PLAN.md`](DELTA_GAIN_PLAN.md) · Living document: update §10 as stages complete.

## 1. Question

> For a defined protein family F (a Pfam domain, a CAZy family/subfamily, a MEROPS peptidase family),
> do the **newly annotated genomes** add sequence diversity that the **baseline genomes** did not already contain?

The whole-proteome funnel (`DELTA_GAIN_PLAN.md`) asks this for *entire proteins* at 95/70/50 % identity. This plan asks it
**inside one family at a time**, at finer resolution, using family-defined *regions* (domains / catalytic modules) rather than whole proteins.
The two views are complementary:

| | Funnel (whole protein) | Family walk (this plan) |
|---|---|---|
| Unit | whole protein, 80 % bidirectional coverage | family region (HMM alignment span / MEROPS hit span) |
| Multi-domain proteins | fail the coverage test if architecture differs → counted "novel" | domain extracted → architecture ignored, then analysed separately |
| Resolution | 95 / 70 / 50 % identity | 90 / 70 / 50 / 30 % identity **plus phylogeny** |
| Annotation source for baseline | NCBI proteins (different predictor) | **BFD re-annotation for both groups** (same pipeline) |
| Interpretation | "is the protein new to the baseline?" | "is the family's sequence/phylogenetic space larger?" |

The BFD-vs-NCBI predictor difference is a confounder for the funnel but **not** here, because Pfam/CAZy/MEROPS tables
were computed on the BFD proteins of all 23.7 k genomes. Baseline vs new is defined by *genome class*
(`genome_classes.tsv`), not by where the protein sequence came from.

## 2. Definitions

* **Baseline genome** – NCBI assembly in the Fungi_BFD set with a prior NCBI annotation and protein file (5,703).
* **New genome** – no prior NCBI annotation (17,864 matched in the manifest; +196 `not_annotated`, +2 `ambiguous_tag` kept out of headline numbers).
* **Family member** – a protein with ≥ 1 hit to F passing the family-specific threshold (§4).
* **Family region** – the aligned span of the member to the family model (Pfam `ali_from..ali_to`; CAZy `s_start..s_end`; MEROPS `q_start..q_end`).
* **Region cluster at identity *t*** – MMseqs2 cluster of family regions within F (`--min-seq-id t -c 0.8 --cov-mode 0`).
* **New-only cluster** – a region cluster at tier *t* that contains **no** baseline member. Counted by member, not by representative, to avoid rep-ownership bias.
* **Added diversity (ΔD)** – gain in a diversity metric when new genomes are added to the baseline, measured against a **size-matched null** (§5).

## 3. Inputs (all on HPCC, `/bigdata/stajichlab/shared/projects/BFD/`)

Tables in `Fungi_BFD_runs/tables/` (parquet files dated 2026-11-05 21:41 … 11-06 09:09; **use the parquet files — `db/BFD.duckdb` is a stale Sep-9 82 MB file**). Row counts verified by review 2026-11-06:

| Table | Rows | Key columns | Notes |
|---|---:|---|---|
| `pfam.parquet` (20 GB) | 368,209,721 | `protein_id, pfam_id, pfam_acc, domain_i_evalue, domain_score, hmm_from/to, ali_from/to, env_from/to, domain_num(_of)` | Pfam 38.2, `hmmsearch --cut_ga` (verified; 0 rows below GA). **Do not add an extra e-value cut** (drops 8.4 % of GA-passing hits) |
| `cazy.cazymes_hmm.parquet` | 9,164,068 | `species_prefix, HMM_id, protein_id, evalue, coverage, s_start/s_end` | dbcanlight HMM; filters (E ≤ 1e-15, cov ≥ 0.35) already applied. **Family = `HMM_id` with `_N` subfamily suffix stripped** (`GH5_7` → `GH5`) |
| `cazy.overview.parquet` | 9,310,657 | `protein_id, EC, cazyme_fam, sub_fam, diamond_fam, substrate, toolcount` | `diamond_fam` is `-` in all rows (no DIAMOND run); `toolcount` ∈ {1,2} = HMM ± dbCAN-sub agreement |
| `merops.parquet` | 43,462,876 | `protein_id, merops_id, percent_identity, aln_length, q_start/q_end, evalue, bitscore` | BLASTP vs `merops_scan.lib` release 124, `-evalue 1e-10 -max_target_seqs 10 -seg yes` (~6.4 hits/protein → choose best hit per protein per family). Family map: `/srv/projects/db/MEROPS/2023-01/merops_lib.families.tab` (all 3,664 MER ids map) |
| `species.parquet` | 23,680 | `ASMID, LOCUSTAG, PHYLUM…SPECIES, NCBI_TAXONID` | `protein_id` prefix = `LOCUSTAG` (join ≈ 100 %; CAZy has 17–22 unmatched rows from prefix `F0B62757`) |
| `asm_stats`, `busco_genome` | 23,584 / 23,582 | assembly + BUSCO | covariates for confounder control |

**Usable genome counts for this analysis** (genomes with BFD proteins): baseline ≈ **5,537** (not 5,703 — that is `baseline_ok` over all NCBI assemblies; `genome_classes.tsv` has 5,584 baseline rows, 5,549 matched); new ≈ **17,838** with Pfam rows. Genomes in classes `annotated_pending` (10), `not_in_ncbi_table` (24) and `new/not_annotated` (11 with CAZy rows) must be excluded explicitly. Ratio new:baseline ≈ 3.2×.

**Taxonomic structure of the new set (review finding):** 10,673 of 17,864 new genomes (60 %) belong to species that already have a baseline genome (e.g. *S. cerevisiae* 149 base / 1,148 new; *C. auris* 21 / 930; *F. oxysporum* 61 / 750); only 64 new genomes lie in the 23 orders with no baseline genome; 57 orders have baseline but no new genomes.
| `swissprot_annot.parquet` | 575,503 | `accession, ec_numbers, pfam, go_terms` | characterised-protein anchors |

Other: `DeltaGain_Fungi/results/delta_funnel/genome_classes.tsv` (baseline/new), funnel tier clusters (63.8 M BASE + 183.1 M NEW input proteins), BFD protein FASTAs (`combine_bfd_proteins.py`).

Join path (all steps are DuckDB, no sequences needed until §6 stage F1):

```
family hit table.protein_id ──prefix──► species.LOCUSTAG ──► species.ASMID ──► genome_classes.class {baseline,new}
                                                         └──► taxonomy (PHYLUM … SPECIES)
```

## 4. Workflow

```
                    ┌─────────────────────────────────────────────────────────────┐
 pfam / cazy /      │ F0  GLOBAL CENSUS  (DuckDB, one pass, minutes–hours)        │
 merops tables  ───►│  per family: n_genomes_base / n_genomes_new, copy number,   │
 + genome classes   │  prevalence, new-only families, taxonomic spread, effect    │
 + taxonomy         │  sizes vs size-matched null                                 │
                    └───────────────┬─────────────────────────────────────────────┘
                                    │ rank families: (a) strongest signal  (b) controls  (c) biology interest
                                    ▼
                    ┌─────────────────────────────────────────────────────────────┐
                    │ F1  EXTRACT family regions (per chosen family)              │
                    │  ids from table ► seqkit/mmseqs createsubdb on BFD proteome │
                    │  trim to ali span (± 10 aa), drop < 80 % of model length    │
                    └───────────────┬─────────────────────────────────────────────┘
                                    ▼
                    ┌─────────────────────────────────────────────────────────────┐
                    │ F2  REGION CLUSTERING within family (MMseqs2, nested)       │
                    │  90 → 70 → 50 → 30 % id, 80 % bidirectional coverage        │
                    │  label each cluster: n_base, n_new, n_genomes, taxa         │
                    └───────────────┬─────────────────────────────────────────────┘
                                    ▼
              ┌─────────────────────┴──────────────────────┐
              ▼                                            ▼
 ┌─────────────────────────────┐              ┌──────────────────────────────────┐
 │ F3  RAREFACTION / NULLS     │              │ F4  PHYLOGENY (per family)       │
 │  new→base size-matched      │              │  reps (50–70 %) ► MAFFT ► trim   │
 │  baseline-vs-baseline null  │              │  ► FastTree/IQ-TREE ► Faith PD   │
 │  taxon-stratified           │              │  new-only clades, branch length  │
 └──────────────┬──────────────┘              └───────────────┬──────────────────┘
                └──────────────────┬──────────────────────────┘
                                   ▼
                    ┌─────────────────────────────────────────────────────────────┐
                    │ F5  ARCHITECTURE + DIVERGENCE + ANCHORS                      │
                    │  domain architectures new/old; score-near-threshold members; │
                    │  catalytic residues; place on SwissProt/CAZy characterised   │
                    └───────────────┬─────────────────────────────────────────────┘
                                    ▼
                    ┌─────────────────────────────────────────────────────────────┐
                    │ F6  JOIN WITH FUNNEL: which survivors fall in which families │
                    │  (funnel survivors_clusters.tsv ⨝ family labels)             │
                    └───────────────┬─────────────────────────────────────────────┘
                                    ▼
                    ┌─────────────────────────────────────────────────────────────┐
                    │ F7  PRIORITISE for structure (ColabFold/Foldseek) & curation │
                    └─────────────────────────────────────────────────────────────┘
```

### Family-type specifics

| | Pfam | CAZy (dbCAN) | MEROPS |
|---|---|---|---|
| Family unit | `pfam_acc` | `HMM_id` family; refine to `sub_fam` (e-family) and `substrate`/`EC` | MEROPS family (S8, M14 …) and clan, from MER id |
| Member filter | **GA-gated already; no extra cut.** Optional sensitivity: model coverage ≥ 50 % and ≥ 80 % | already filtered (E ≤ 1e-15, cov ≥ 0.35); strict set = `toolcount = 2` (HMM and dbCAN-sub agree) | best hit per protein **per family**; 22 k proteins hit > 1 family; `aln_length ≥ 50 % of subject` needs subject lengths from `merops_scan.lib` (add) |
| Region | `env_from..env_to` (one region per domain instance; repeats kept separate, report merged-envelope sensitivity) | `s_start..s_end` | `q_start..q_end` |
| Divergence signal | `domain_score` relative to the family GA | HMM `coverage`; sub-family assignment vs family-only | `percent_identity` of best MEROPS hit (low = distant) |
| Extra | architecture string from `domain_num` ordering | substrate / EC diversity per subfamily | catalytic-triad / active-site motif per family |

## 5. Metrics and null models

Raw "new genomes contain X more clusters" is **meaningless on its own**: there are ~3.2× more new genomes than baseline
genomes, they are taxonomically skewed, and **60 % of them are additional strains of species already in the baseline** — so at the 90/70 % tiers
much "novelty" is strain/population variation or paralog splitting, not new biology. Every claim must be against a null, and
**every metric is reported at three dereplication levels: all genomes, one genome per species (headline), one per genus; and split into
"new taxon" (species/genus absent from baseline) vs "within known species".**

1. **Prevalence delta** – families present in ≥ *k* new genomes and ≤ *j* baseline (k=3, j=0 as "new-only families"); FDR-controlled across families (Benjamini–Hochberg in F0).
2. **Rarefied cluster richness** – rarefy baseline and new to the **same n within the orders they share** (a 5,703-genome order-matched draw from new is impossible: Σ min(base,new) over orders is only 4,873); use species-collapsed accumulation curves for baseline alone vs baseline+new. Report ratio and 95 % interval.
3. **Leave-k-out null (baseline vs baseline)** – hold out *k* baseline genomes taxon-matched to a random *k*-subset of new genomes, compare their novelty against the remaining baseline at equal *k* (a single-genome leave-one-out does not model adding thousands of genomes). Run species-collapsed. New genomes are interesting only if they exceed this distribution.
4. **Per-genome marginal gain** – order genomes randomly; accumulate; fit saturation curve (re-use `accumulation_curve.py` / `lib_accumulation.py`); baseline curve vs baseline+new.
5. **Phylogenetic diversity added** – Faith PD (sum of branch length) contributed by new members on the family tree vs the same for the size-matched null; plus fraction of **new-only clades** (≥ 3 new members, no baseline member).
6. **Distant-member fraction** – share of new members whose score lies in the lowest decile of the baseline score distribution (Pfam/CAZy) or whose best-MEROPS identity < 40 %.
7. **Architecture novelty** – number of distinct domain architectures observed only in new genomes.
8. **Subfamily / substrate novelty (CAZy)** – e-families, EC numbers and substrate labels seen only in new genomes; placement of new clusters relative to SwissProt/characterised members.

Counting rule: members, never representatives. Identity-tier results always reported for all tiers (a signal only at 30 % means distant homologs, only at 90 % means near-duplicates).

## 6. Pilot sequence

Run in this order so each stage can fail cheaply.

| Step | Family | Why | Expected size |
|---|---|---|---|
| **P0 (control, −)** | **Proteasome `PF00227`** (F0: 5,534 base / 17,789 new genomes, ~14 copies/genome) and **Tubulin `PF00091`** (4.6 copies; species-level q = 1). *(The review named `PF00246` as the proteasome — that is Peptidase_M14, i.e. P3's contrast family.)* | calibrates the null: strain variation will show "novelty" here too, so the pass criterion is meaningful only after species collapse | PF00227 ≈ 10⁵ proteins |
| **P1 (CAZy)** | **GH18** (`GH18` dbCAN, strip subfamily suffix, + `PF00704`) | large, diverse; cross-checks CAZy vs Pfam views of one family | GH18 253,602; PF00704 266,000 |
| **P2 (CAZy, +)** | **AA9** (LPMO; `PF03443`) | known diversified/expanded family; positive control | AA9 192,879; PF03443 197,003 (~20× the first estimate) |
| **P3 (MEROPS)** | **S8** subtilisin-like (and **M14** as a contrast) | MEROPS BLAST-hit path; S8 varies with lifestyle | S8 265,246; M14 57,804 |
| **P4 (Pfam, +)** — **deferred until after shutdown** | `PF00172` (3.21 M proteins), `PF00067` (2.02 M) | exceed the 2 M linclust trigger; stress tests only | 10⁶ |
| **P5 (Pfam, dark)** | A set of **DUF** families concentrated in new genomes (from F0) | links to "dark proteome" / structure goal | varies |

A family qualifies as a pilot only after F0 gives its size; if > 2 M regions, cluster at 90 % first (linclust) and carry reps forward. **Preregister the pilot list** (above) before looking at F0 results to avoid selecting on the outcome.

## 7. Measures of success

| Level | Criterion |
|---|---|
| Pipeline valid | P0 control: new-vs-null difference within the baseline-vs-baseline null at all tiers (|z| < 2) |
| Method useful | ≥ 1 family where new genomes exceed the baseline-vs-baseline null at ≥ 2 tiers **and** gain PD on the tree |
| Biologically real | Those signals persist after (a) removing assemblies in the lowest BUSCO quartile, (b) requiring ≥ 80 % model coverage, (c) taxon stratification |
| Reproducible | Whole table F0 rerunnable in one command; per-family stages F1–F5 driven by a single `family_id` arg |
| Useful to the funnel | Family labels attached to ≥ 90 % of funnel survivors that carry a Pfam hit; list of families that concentrate novelty |

A **negative** result (new genomes add little beyond the null) is a valid outcome and informs how much of the 18 k-genome annotation effort adds family-level novelty versus redundancy.

## 8. Risks and confounders

| ID | Risk | Mitigation |
|---|---|---|
| R1 | New genomes dominate counts (3.1×) | rarefaction + leave-one-out null (§5) |
| R2 | Taxonomic skew (new = environmental / under-sampled lineages) | taxon-stratified sampling; report per-order tables; treat "novelty from new taxa" as an *expected* source, separate from "within-taxon" novelty |
| R3 | ~~MEROPS family map missing~~ **Closed**: `/srv/projects/db/MEROPS/2023-01/merops_lib.families.tab` maps all 3,664 MER ids (5,097 rows; headers also carry `#S01A#`) | join on `merops_id`; best hit per protein per family |
| R10 | **60 % of new genomes are additional strains of baseline species** → 90/70 % tiers measure population variation | species/genus dereplication as headline; "new taxon" vs "within known species" split (§5) |
| R11 | Baseline genome count differs by table (5,703 vs 5,537 usable) | always compute n from the join, never hard-code |
| R12 | Multiple testing across ~16 k Pfam + CAZy + MEROPS families | BH-FDR in F0; preregistered pilot list |
| R13 | Gene-model quality may differ between baseline and new (e.g. RNA-seq availability) | covariate check; stratify by BUSCO/assembly stats |
| R4 | Fragmented/partial genes inflate "distinct" sequence | require ≥ 80 % of model length; drop `gene model partial` flags; compare with and without |
| R5 | Contamination / HGT-like bacterial proteins | taxonomy filter via best DIAMOND/`taxonomy` columns; inspect new-only clades manually |
| R6 | ~~Pfam run not TC-gated~~ **Closed**: Pfam 38.2 `hmmsearch --cut_ga` (`Fungi_BFD/nextflow/modules/BFD/PFAM/main.nf`); 0 rows below GA | no extra filtering |
| R7 | Multi-copy families → paralog inflation | report per-genome copy number; cluster-level counts use distinct genomes as well as members |
| R8 | Assembly quality heterogeneity | covariates from `busco_genome`, `asm_stats`; stratify by BUSCO completeness |
| R9 | ID mapping between BFD ids and funnel ids (`NEW__<ASMID>__<orig>`) | check on 1 genome before F6 |

## 9. Compute and time budget

**Time is the binding constraint: the cluster is scheduled to shut down 2026-11-10 08:00.** Anything that must read the 20 GB Pfam table
or the BFD proteomes has to run before then, or the inputs must be copied out first.

| Stage | Compute | Needs cluster? |
|---|---|---|
| F0 global census | DuckDB over parquet; ~minutes–2 h, 32 cores, ≤ 200 GB, `highclock` (AVX2) | yes (data are here) |
| F1 region extraction (pilots P0–P4) | seconds–minutes per family | yes |
| F2 clustering | MMseqs2; tiny for most, hours for P4 | yes (or later, if regions are exported) |
| F3 null models | DuckDB/numpy, small | no (after F0 export) |
| F4 trees | MAFFT+FastTree on ≤ 20 k reps, hours | after export OK |
| F5–F7 | small / structure prediction elsewhere | no |

**Priority before shutdown (revised after review):**

| Day | Work |
|---|---|
| Day 1 (Nov 6–7) | F0 as one DuckDB pass over the parquet files (GA-gated Pfam; CAZy with subfamily suffix stripped; MEROPS via `families.tab`), output per-family genome counts per class **and per species-collapsed set** plus copy number; export immediately |
| Day 2 | **One** extraction job that reads the 23 k per-genome BFD protein files once (`Fungi_BFD_runs/input/pep/<tag>.proteins.fa.gz`; ~1–2 h I/O) and writes regions + membership TSVs for P0–P3 and the top F0 DUFs. **Do not reuse `funnel_input.faa`** — its baseline sequences are NCBI proteins, not BFD. Copy regions, memberships, `species`, `genome_classes`, `asm_stats`, `busco_genome`, `merops_lib.families.tab` off-cluster |
| Day 3 | MMseqs2 tiers on pilot families (minutes each); optionally extract regions for the top 50–100 F0 families as an archive |
| Also export | funnel `tier_*_clusters.tsv` + provenance (without them F6 cannot be done later; check the BASE-id (NCBI) → family-table mapping first) |

**After the shutdown:** nulls/rarefaction, trees and Faith PD (not feasible at 2×10⁵ members before Nov 10; PD on a FastTree of 50 % reps is fragile anyway), F5–F7, funnel join (F6). **Cut before shutdown:** P4, tree building, Faith PD. Disk: /bigdata 94 % full, 19 TB free.

HPC rule: **never abu_dhabi nodes** (no AVX2; DuckDB/MMseqs die with "Illegal instruction"). Use `-p highclock` or `-p epyc,intel` with `--constraint="ryzen|broadwell|cascade|sapphire"`. `highmem` h01–h06 are Ivy Bridge (no AVX2) as well.

## 10. Status log (update as you go)

| Date | Stage | Result |
|---|---|---|
| 2026-11-06 | Plan drafted; table schemas verified (§3) | — |
| 2026-11-06 | Independent (Fable) review + HPCC verification | Corrected: Pfam GA-gated (R6 closed), MEROPS family map exists (R3 closed), CAZy family unit/`diamond_fam`/thresholds, baseline n ≈ 5,537, pilot sizes, nulls redesigned for 60 % same-species new genomes, FDR + preregistration, Day 1–3 priorities |
| 2026-11-06 | **F0 census done** (`bin/family_census.py`, job 31779947, 1m29s) | 5,549 baseline / 17,864 new genomes; species with baseline 2,755, with new 5,454, **new-taxon species 4,632 (85 % of new species)** — while 60 % of new *genomes* are strains of known species (a few species hold very many genomes). Families: 16,053 Pfam, 6,799 CAZy (family/subfamily/e-family), 412 MEROPS. "New-only" (≥3 new, 0 baseline genomes): 1,597 Pfam, 17 CAZy families, 299 e-families, 16 MEROPS families. **Species-level q < 0.05 for 5,721 Pfam families** — most have tiny effect sizes (e.g. proteasome, ribosomal: log2 ratio ≈ −0.01 but q ≈ 10⁻³) so significance alone is uninformative; use effect size + the F3 nulls. GH18: genome-level p = 7e-10 but species-level q = 0.07 (illustrates the strain effect) |
| 2026-11-06 | **F1 extraction done** (`bin/family_extract.py`, job 31780033, 8 min, 738 MB in `results/family_regions/`) | 18 family files: PF00227 363,233 regions / 333,095 proteins (5,534 base / 17,789 new genomes — identical to F0); GH18 253,791 / 253,194; AA9 193,304; S08 264,977; M14 57,710; PF00091 104,159; PF00704 268,146; PF03443 198,613; 10 DUFs (PF07799 DUF1643 … PF28466 DUF8407; 156–6,893 regions). 23,366 genome files read, 0 missing, 0 unmatched ids |
| | P0 control | not started |
| | P1 GH18 | not started |

## 11. Next actions (handoff checklist)

1. `bin/family_census.py` — DuckDB script producing `results/family_census/{pfam,cazy,merops}_census.tsv` (per family: genomes base/new, copies, rarefied expectation, new-only flag). Verify the Pfam TC gating first (R6).
2. Obtain and load the MEROPS id→family map (R3).
3. `bin/family_extract.py <family>` — region FASTA + membership TSV (genome class, taxonomy, score).
4. `bin/family_cluster.sh <family>` — nested MMseqs2 tiers, labelled by baseline/new.
5. `bin/family_nulls.py <family>` — rarefaction and leave-one-out null.
6. Pilot P0 → P1 → P2 → P3 and record results in §10.
7. When the funnel summary completes (`results/delta_funnel/summary/`), run F6.
