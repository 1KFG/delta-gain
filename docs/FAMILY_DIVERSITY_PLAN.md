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
| 2026-11-06 | **F2 clustering done** (`bin/family_cluster.py`, array 31780205, 18 tasks ≤ 72 s each, `results/family_clusters/`, `tier_summary_all.tsv`) | Fraction of new regions whose cluster has **no baseline member** (overall / within-known-species / new-taxon), 90 → 70 → 50 → 30 %: **PF00227 proteasome (control)** 9.9/1.3/22 % → 0.7/0.3/1.3 % → 0.2/0.2/0.4 % → 0.1/0.1/0.2 %; **PF00091 tubulin** 8.8/1.5/20 % → 1.0/0.3/2.0 % → 0.2 % → 0.1 %; **GH18** 21.4/5.1/43.7 % → 7.0/1.3/14.9 % → 1.5/0.3/3.0 % → 0.15 %; **AA9** 15.6/2.7/34.8 % → 3.2/0.4/7.2 % → 0.6/0.1/1.4 %; **S08** 23.1/6.1/50.9 % → 8.6/1.8/19.6 % → 3.3/0.9/7.1 % → 0.7/0.3/1.3 %; **M14** 13.7/2.2/33 % → 3.9/0.6/9.4 % → 1.1/0.3/2.4 %. Species-collapsed pair fractions follow the same ordering. **Reading:** novelty concentrates in the new-taxon stratum (most of it is simply taxa absent from the baseline); within known species it is small for the control (≈ 1 % at 90 %, ≈ 0.3 % at 70 %) and 2–5× larger for GH18/S08/AA9 at 90–70 %. S08 and GH18 stay above the control down to 50 %. **Not yet established:** whether these exceed what an equal number of *baseline* genomes would add (F3 baseline-vs-baseline null) — until then these are descriptive. Control also shows that 90 %-tier "novelty" in a conserved family is mostly taxon sampling |
| 2026-11-06 | **F3 baseline-vs-baseline null done** (`bin/family_null.py`, array 31781060, 18 tasks, `results/family_null/null_all.tsv`) | Order-matched paired draws (H = held-out baseline genomes, N = equally many new genomes, both scored against baseline∖H; 114 shared orders, ~490 genomes each at f = 0.2; 200 draws; genome and species-collapsed variants; 18 families × 4 tiers = 72 rows per setting). **New genomes add essentially the same novelty as equally many same-order baseline genomes**: no row has a diff CI wholly below 0; CI wholly above 0 in only 0–3 of 72 rows per setting; mean diff is positive in ~70–80 % of rows but small (about +0.5 to +2.5 points). Example (species-collapsed, f = 0.2, c90): S08 new 50.8 % vs null 46.3 % (diff +4.5 [−0.7, +9.1]); GH18 43.0 vs 40.7; AA9 37.0 vs 35.9; **control PF00227 20.0 vs 18.2 (+1.8 [0.0, +3.5])**. **The family ranking seen in F2 (S08 > GH18 > AA9 > control) is already present in the null (S08 46 % vs control 18 % at c90) — it reflects inherent family diversity, not anything specific to the new genomes.** Known-species stratum: diff positive in ~50–63 % of rows (genome variant: CI > 0 in 11/68 at f = 0.2; species-collapsed: 0/68) — consistent with a small strain/population effect that disappears under species collapse. **Scope limit:** only orders present in both sets are tested (57 baseline-only orders and the many new-only orders/species are excluded), so this says nothing about taxa the baseline lacks — that needs a leave-one-order-out null (F3b). Draw-level "p" values are the fraction of paired draws with diff ≤ 0, not a formal test |
| | P0 control | done (F2 + F3): control shows the same small positive shift as pilots, so it does not separate "pilot" from "control" |
| 2026-11-06 | **F3b leave-one-order-out null done** (`bin/family_null_order.py`, array 31785571, `results/family_null_order/{ordnull_all,order_null_summary}.tsv`) | Reference for each order o = baseline without order o. `base_out` = K = 10 baseline genomes drawn from o (orders with >= 20 baseline genomes, 31-43 per family); `new` = K = 10 new genomes from o (>= 10 new genomes); `new_only_pool` = 10 new genomes drawn from the 64 genomes in the 23 orders with no baseline (too few per order, pooled; reference = whole baseline); 200 draws, genome and species-collapsed variants. **Holding out a whole order is far more novel than the within-order F3 null** (c70 median: 60-85 % for S08/GH18/AA9/M14, 1-4 % for tubulin/proteasome controls) and **new genomes in baseline orders look the same as held-out baseline orders** (e.g. S08 c70 68 % vs 71 %, GH18 71 % vs 70 %). **Genomes from baseline-less orders are NOT more novel than a typical held-out baseline order**: at c70 the pool sits 6-30 points *below* the base_out median for the diverse families (the two low-diversity controls are a little above: PF00091 +7, PF00227 +0.5) (S08 56 % vs 71 %, GH18 56 % vs 70 %, M14 58 % vs 85 %, PF03443 36 % vs 68 %); at c50 differences are within about +/-6 points, and at c30 the pool is a few points above (S08 5.4 % vs 1.2 %, M14 5.3 % vs 0 %). Interpretation: the family-level novelty of new taxa is what you get from sampling any unseen clade, not a special property of the new genomes. **Caveats:** the pool mixes 23 orders (and may include close relatives of covered orders; mixing lowers novelty relative to a single clade), uses only one pooled comparison per family, and the 64 genomes are few; the low-diversity Pfam DUF families (PF07799, PF14024 ...) are saturated at 1.0 and uninformative |
| 2026-11-06 | **Wide run (F1-F3b on 717 families)**: 340 CAZy + 174 MEROPS (all) + top 200 Pfam by F0 species-level enrichment (+ pilots/DUFs); `results/family_{regions,clusters,null,null_order}_wide` (F1 20 min, 18.0 M regions; 568 families clustered, 149 skipped for < 30 regions); summaries use families with >= 1000 regions (192 CAZy, 91 MEROPS, 192 Pfam), species-collapsed, f = 0.2 | **Unselected families (CAZy, MEROPS): new genomes add almost nothing beyond equally many same-order baseline genomes** (c70 mean diff +0.7 / +1.7 points; diff CI wholly > 0 in 1.6 % / 9.9 % of families, never < 0; c90 +2.3 / +3.1 points, 6 % / 17 % CI > 0) - the pilot conclusion holds across families. **Enrichment-selected Pfam: clearly positive** (c70 mean diff +14.6 points, 97 % of families positive, but CI > 0 in only 2.6 %) - expected given the families were chosen for being enriched in new genomes, so this is a best case, not an estimate; **needs a size-matched random-Pfam control** before it is interpreted. F3b: baseline-less-order pool vs held-out baseline orders: CAZy/MEROPS c70 -6.8 / -6.4 points (pool above in 30 % / 28 % of families), c50 +3-4 points (68-69 % above), c90 -5 points; new genomes in baseline orders equal held-out baseline orders (median diff 0 to -1.5 points). Pfam F3b rows are mostly saturated (median difference 0.0) |
| 2026-11-07 | **Random-Pfam control for the wide run** (198 families extracted, 197 clustered/nulled; `results/pfam_control_families.txt`, `family_*_ctrl`): 200 Pfam families with >= 100 new genomes and NOT significantly enriched in new genomes (pool 6,376), each matched to one enriched family by nearest (log g_new, log g_base) - match is loose (enriched mean 3,842 new / 397 baseline genomes vs control 2,588 / 612); same F3 settings, families >= 1000 regions | **The +15-point effect is specific to the enriched families**: c70 mean F3 diff **+14.9 points (98 % positive)** for the enriched 200 vs **+2.2 points (66 % positive)** for the controls; c90 +15.1 vs +2.8; c50 +7.1 vs +1.6. CI wholly > 0: 2.7 % vs 0.5 % (c70). So an average Pfam family behaves like the CAZy/MEROPS families (new genomes add ~1-3 points beyond same-order baseline genomes), while the families that are over-represented in new genomes also carry markedly more novel sequence variants. We cannot separate selection (families picked because they are over-represented) from a real effect of new genomes, so the enriched set is a **candidate list for follow-up**, not a genome-wide estimate. Many enriched families are small, lineage-restricted or viral-like (e.g. PF00519), so check contamination/fragmentation before interpreting |
| 2026-11-07 | **F4 artefact/quality diagnostics** (`bin/family_qc.py`, `results/family_qc{,_ctrl}/`: `family_qc.tsv`, `qc_by_quality.tsv`, `qc_within_order.tsv`) for enriched Pfam, control Pfam, CAZy, MEROPS, tiers c70/c50 | **(1) The enriched Pfam list is a Saccharomycotina yeast list**: of its new-genome regions 453k are Saccharomycetales, 421k Serinales, 70k Pichiales (Ascomycota in 198/199 families); BUSCO of those orders is 75-80 % (fungi_odb gene loss), so 99.7 % of enriched-family regions sit in BUSCO < 95 genomes - a lineage effect, not a quality flag. Enrichment is also driven by sampling depth (new genomes include ~2,000 Saccharomycetales and ~2,200 Serinales). **(2) Raw novelty is 2-3x higher in BUSCO < 95 / N50 < 50 kb genomes for CAZy/MEROPS (c70 e.g. 12.6 % vs 4.9 %), but this is mostly clade composition**: compared within the same order the low-minus-high BUSCO gap is only +2.2 points (CAZy), +4.6 (MEROPS), +2.2 (control Pfam) and N50 gap +0.2-0.4 points; enriched Pfam shows the opposite sign for N50 (low-N50 genomes less novel, -9 points) so assembly contiguity does not explain the enriched-family novelty. **(3) Other checks (c70, new regions)**: top-species share of novel regions is small (median 0.02 enriched, 0.07 CAZy, 0.06 MEROPS: not one-species driven); ~38-50 % of novel regions are in clusters private to one species (so a species-level count is still needed); truncated regions are not enriched among novel regions for Pfam (0.16 novel vs 0.18 known) or MEROPS (0.17 vs 0.16) and only mildly for CAZy (0.050 vs 0.034); model coverage of novel Pfam/CAZy regions is lower (0.84 vs 0.90 enriched; 0.82 vs 0.90 CAZy), consistent with partial or divergent domains |
| 2026-11-07 | **F5 deep dive of S08, GH18, AA9, M14 at c70** (`bin/family_deepdive.py`, `results/family_deepdive/` log; new-genome regions, novel = no baseline member in the cluster) | Novelty is **clade-specific, not family-wide**: GH18 7.0 % overall (5,253 baseline-free clusters, 3,281 species) but 2.9 % in Hypocreales and 2.1 % in Eurotiales vs 72.8 % Ascoideales, 42.6 % Dipodascales, 34.4 % Pichiales, 23.0 % Boletales; S08 5.0 % in Hypocreales vs 31.1 % Agaricales, 38.2 % Phaffomycetales, 22.0 % Dipodascales; AA9 3.1 % overall, 8.2 % Agaricales; M14 2.6-4.3 % in the main subfamilies, 22.2 % Agaricales. Per MEROPS subfamily (MER id of the best hit) S08 novelty runs from 14-41 % (small subfamilies, 200-250 regions) to 17.6 % for the largest (MER0000358, 6,921 regions, 656 species); M14 from 1 to 20 %. Baseline-free clusters are small: mean 1.5-1.7 species per cluster (S08 1.46, GH18 1.70, AA9 1.54, M14 1.67). **Interpretation: new diversity concentrates where baseline sampling is thin (Agaricales/Boletales, early-diverging Saccharomycotina) and is almost absent in well-sampled Hypocreales/Eurotiales.** Not done: catalytic-residue and architecture analysis (would need alignments and domain-architecture tables) |
| | F3c (optional) | per-order novelty ranked against the null, and a clade-matched null for the 23 baseline-less orders (closest-covered-order distance) |
| | P1 GH18 | not started |

## 11. Next actions (handoff checklist)

1. `bin/family_census.py` — DuckDB script producing `results/family_census/{pfam,cazy,merops}_census.tsv` (per family: genomes base/new, copies, rarefied expectation, new-only flag). Verify the Pfam TC gating first (R6).
2. Obtain and load the MEROPS id→family map (R3).
3. `bin/family_extract.py <family>` — region FASTA + membership TSV (genome class, taxonomy, score).
4. `bin/family_cluster.sh <family>` — nested MMseqs2 tiers, labelled by baseline/new.
5. `bin/family_nulls.py <family>` — rarefaction and leave-one-out null.
6. Pilot P0 → P1 → P2 → P3 and record results in §10.
7. When the funnel summary completes (`results/delta_funnel/summary/`), run F6.


## Extension ideas from the review panel (2026-11-07) and what we are doing with them

Three read-only consultants proposed work plans (statistical design, biology, artefact attack). Merged, ranked by value/cost before the 2026-11-10 shutdown. "Primary" = the planned 1-6 questions below.

| # | Idea | Source | Decision |
|---|---|---|---|
| E1 | **Predictor-mismatch control for the funnel**: same T baseline genomes annotated by NCBI vs BFD, same reference; novelty difference = predictor effect (`bin/predictor_control_prep.py`, `launch/predictor_control.sbatch`) | artefact (top risk) | **done**: predictor effect +1.5-3 points (see DELTA_GAIN_PLAN §11) |
| E2 | **Rank decomposition of funnel novelty**: tag each new genome by the lowest rank shared with a baseline genome (strain / species-in-genus / genus-in-family / higher) and report novel proteins per genome per stratum, species-collapsed | design, artefact | **done** (DELTA_GAIN_PLAN §11) |
| E3 | Rarefy new genomes to k per species (k = 1, 2, 5) and re-score novelty vs baseline | design | do after E2 |
| E4 | Species accumulation curves, baseline species first, then new species (marginal gain per species) | design | queue |
| E5 | Contamination triage of funnel survivors (length, intron density, GC/codon outliers, genomes with extreme survivor counts) | biology, artefact | queue (part of primary item 2) |
| E6 | Orphan vs taxon-restricted survivors (species / genus / order restricted), >= 2 species required | biology | queue |
| E7 | Lineage-specific expansions (copy number) and secretome/SSP novelty; TE-domain flags | biology | queue, lower priority |
| E8 | Split-sample test of Pfam enrichment (define on half the species, measure on the other half); permutation FDR | artefact | queue |
| E9 | Parameter sensitivity of the family runs (cov-mode, -c, ali vs env) on pilot families | artefact | optional |
| E10 | Chao/Good-Turing, capture-recapture, effective genome number, HGT calling, antiSMASH | design, biology | **not planned**: assumptions fail (non-random, clonal samples), or no table/compute before shutdown |
