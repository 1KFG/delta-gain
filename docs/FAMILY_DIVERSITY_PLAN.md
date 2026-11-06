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

Tables in `Fungi_BFD_runs/tables/` (rebuilt 2026-11-05/06; the consolidated `db/BFD.duckdb` is being rebuilt – **re-verify names/counts when it lands**):

| Table | Rows | Key columns | Notes |
|---|---:|---|---|
| `pfam.parquet` (20 GB) | 368,209,721 | `protein_id, pfam_id, pfam_acc, domain_i_evalue, domain_score, hmm_from/to, ali_from/to, domain_num(_of)` | **verify** the run used `--cut_tc`; if not, apply per-family TC from Pfam-A.hmm.dat |
| `cazy.cazymes_hmm.parquet` | 9,164,068 | `species_prefix, HMM_id, protein_id, evalue, coverage, s_start/s_end` | dbCAN HMM family (e.g. `GH27`) |
| `cazy.overview.parquet` | 9,310,657 | `protein_id, EC, cazyme_fam, sub_fam, diamond_fam, substrate, toolcount` | `sub_fam` = dbCAN-sub e-family (e.g. `GH27_e9`); `toolcount` = tool consensus |
| `merops.parquet` | 43,462,876 | `protein_id, merops_id, percent_identity, aln_length, q_start/q_end, evalue, bitscore` | BLAST vs MEROPS library; **needs MER-id → family/clan map** (§8, risk R3) |
| `species.parquet` | 23,680 | `ASMID, LOCUSTAG, PHYLUM…SPECIES, NCBI_TAXONID` | `protein_id` prefix = `LOCUSTAG` |
| `asm_stats`, `busco_genome` | 23,693 | assembly + BUSCO | covariates for confounder control |
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
| Member filter | `domain_i_evalue ≤ 1e-5` **and** score ≥ TC (verify) ; drop AntiFam-like hits if present | `evalue ≤ 1e-15` and `coverage ≥ 0.35` (dbCAN default) **and** `toolcount ≥ 2` for the strict set | best hit `evalue ≤ 1e-5`, `aln_length ≥ 50 % of subject` |
| Region | `ali_from..ali_to` | `s_start..s_end` | `q_start..q_end` |
| Divergence signal | `domain_score` / TC | HMM `coverage`, `diamond_fam` vs HMM disagreement | `percent_identity` of best MEROPS hit (low = distant) |
| Extra | architecture string from `domain_num` ordering | substrate / EC diversity per subfamily | catalytic-triad / active-site motif per family |

## 5. Metrics and null models

Raw "new genomes contain X more clusters" is **meaningless on its own**: there are 3.1× more new genomes than baseline
genomes, and they are taxonomically skewed. Every claim must be against a null.

1. **Prevalence delta** – families present in ≥ *k* new genomes and ≤ *j* baseline (k=3, j=0 as "new-only families").
2. **Rarefied cluster richness** – for tier *t*, sample *n* = 5,703 genomes from the new set (1,000 draws, taxon-stratified to mirror the baseline's order-level composition); compare cluster count with the baseline's. Report ratio and 95 % interval.
3. **Leave-one-out null (baseline vs baseline)** – hold out a random baseline genome (or clade); how many *new-only clusters* does it add to the remaining baseline? New genomes are interesting only if they exceed this distribution. This is the key control for "any extra genome adds something".
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
| **P0 (control, −)** | A conserved housekeeping family (e.g. a BUSCO-linked Pfam such as a ribosomal or proteasome subunit) | delta ≈ 0 expected; validates the null and the pipeline | ~1 copy/genome |
| **P1 (CAZy)** | **GH18** (chitinases; `GH18` dbCAN + `PF00704`) | large, diverse, ecologically interesting (mycoparasites, insect pathogens), cross-checks CAZy vs Pfam views of one family | 10⁴–10⁵ members |
| **P2 (CAZy, +)** | **AA9** (LPMO) | known diversified/expanded family; positive control | 10⁴ |
| **P3 (MEROPS)** | **S8** subtilisin-like (and **M14** as a contrast) | MEROPS query-based path; S8 is large and varies with lifestyle | 10⁵ |
| **P4 (Pfam, +)** | Zn(II)₂Cys₆ `PF00172` or cytochrome P450 `PF00067` | very large lineage-specific families; stress test for scale | 10⁵–10⁶ |
| **P5 (Pfam, dark)** | A set of **DUF** families concentrated in new genomes (from F0) | links to "dark proteome" / structure goal | varies |

A family qualifies as a pilot only after F0 gives its size; if > 2 M regions, cluster at 90 % first (linclust) and carry reps forward.

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
| R3 | **MEROPS family map missing** (table has `MER#` ids only) | obtain MEROPS `domain.txt`/`pepunit.lib` family assignment (MEROPS FTP) and load as `merops_family` table; **blocker for P3** |
| R4 | Fragmented/partial genes inflate "distinct" sequence | require ≥ 80 % of model length; drop `gene model partial` flags; compare with and without |
| R5 | Contamination / HGT-like bacterial proteins | taxonomy filter via best DIAMOND/`taxonomy` columns; inspect new-only clades manually |
| R6 | Pfam run not TC-gated | verify; re-threshold with `pfam_len` + per-family TC table before F0 |
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

**Priority before shutdown (in order):** (1) F0 for Pfam + CAZy + MEROPS; (2) F1 FASTA export of regions for P0, P1, P2, P3 (and P4 if time); (3) copy `genome_classes.tsv`, `species.parquet`, the F0 outputs and the region FASTAs to storage that survives the shutdown.

HPC rule: **never abu_dhabi nodes** (no AVX2; DuckDB/MMseqs die with "Illegal instruction"). Use `-p highclock` or `-p epyc,intel` with `--constraint="ryzen|broadwell|cascade|sapphire"`. `highmem` h01–h06 are Ivy Bridge (no AVX2) as well.

## 10. Status log (update as you go)

| Date | Stage | Result |
|---|---|---|
| 2026-11-06 | Plan drafted; table schemas verified (§3) | — |
| | F0 census | not started |
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
