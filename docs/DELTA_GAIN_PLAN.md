# Delta Gain — plan, algorithm, assumptions and measures of success

Status: living document. Last updated 2026-11-05.
Convention: **[M]** = measured on the cluster, **[E]** = estimate (replace when measured), **[TBD]** = not yet known.

---------------------------------------------------------------------------------------------------

## 1. The question

> How much NEW protein-level information do newly annotated fungal genomes add,
> compared with what was already known from previously annotated fungal genomes?

"Delta gain" is the novelty contributed by a set of **new** genomes, measured against an explicit
**baseline** set of genomes. It is *not* "how much is missing from public databases" (that was the
first design, see §12), and it is not a property of any single genome in isolation.

Three levels of answer, in increasing biological content:

1. **Cluster level** — what fraction of new proteins fall in sequence clusters with no baseline member,
   at several identity tiers (95 / 70 / 50 %)?
2. **Family level** — do the unmatched proteins form real families (multiple genomes, multiple genera),
   do they carry known Pfam domains (new architectures) or no known domain (dark)?
3. **Case studies** — for a single Pfam / CAZy family, do the new genomes add new *subfamilies*
   beyond what the baseline genomes show, more than a random-order null predicts?

## 2. Definitions

| Term | Definition |
|---|---|
| **Baseline genome** | Fungal assembly with NCBI-submitted/RefSeq gene models (`annotation_info` present, ≥1 protein-coding gene) **and** a `*_protein.faa.gz` on disk in `NCBI_ASM`. |
| **New genome** | Assembly with **no** prior NCBI annotation; the only gene models are BFD's (funannotate pipeline). |
| **Baseline proteins** | The NCBI protein FASTA of the baseline genomes (their own gene models). |
| **Query proteins** | BFD-annotated proteins of the new genomes. |
| **Novel at tier t** | A query protein whose tier-t cluster contains no baseline protein. Counted **per protein** (member-weighted), never per representative. |
| **Survivor** | A query cluster still baseline-free at the last tier; input to the sensitive search stage. |
| **Pending** | NCBI-annotated assembly whose proteins are not on disk yet. Excluded from both sets so it can never be counted as "new". |
| **Sibling-flagged** | A "new" strain for which a *different assembly of the same strain* has NCBI annotation (dropped by the NCBI_fungi species+strain dedup). Reported separately, excluded from the headline new set. |

## 3. Inputs and expected sizes

```
 NCBI datasets snapshot 2026-10-23  ........................ 27,401 fungal assemblies            [M]
 NCBI_fungi download (1KFG/2026/NCBI_fungi/source/NCBI_ASM) . 23,838 folders (23,115 GCA, 723 GCF)  [M]
 BFD protein manifest (full run) ........................... 23,447 matched, 229 not_annotated,
                                                              4 ambiguous_tag                       [M]
```

Prior-annotation status of the 27,401 assemblies (`DeltaGain_Fungi/prior_annotation_status.tsv`) [M]:

| Status | Assemblies | Role |
|---|---:|---|
| annotated, protein file on disk (`baseline_ok`) | **5,703** | **baseline** (after gap download 2026-11-05; was 5,511) |
| annotated, folder present, protein file missing | 10 | NCBI lists an annotation but ships no protein file; excluded from both sets |
| annotated, never downloaded | 1,184 | same-strain duplicates dropped by the project's dedup (1,118 already have an annotated sibling in the baseline; 57 strains flagged) |
| no prior annotation, downloaded | 18,145 | **new** |
| no prior annotation, not downloaded | 2,359 | not in our set |

Protein counts: **measured** by the full funnel input build of 2026-11-06 [M] (the earlier estimates were 55M / 135-195M):

| Set | Genomes | Proteins/genome | Proteins |
|---|---:|---:|---:|
| Baseline (NCBI proteins) | 5,703 | ~11.2k | **63,766,960** [M] |
| Query (BFD, new genomes) | 17,864 | ~10.2k | **183,095,650** [M] |
| Funnel input (both) | 23,567 | | **246,862,610** [M] (17,715 sequences < 30 aa dropped; 0 duplicate ids; 0 missing files) |
| After tier-0 linclust 95/90 | | | ~85–110M [E; pilot 2.0×, full run 2.25× reduction] |

Reference points already measured:
- 2000-genome pilot: 21M proteins → **10,615,362** representatives (95/90) [M].
- Full 23.4k-genome BFD set, 95/90: 1,307 chunks × 60,000 ≈ **78M** representatives [M]; linclust 1h25m, 224 GB RAM [M].

## 4. Workflow

```
                    NCBI datasets: 27,401 fungal assemblies  [M]
                                      │
              annotation_info + protein.faa on disk?   (prior_annotation_status.tsv)
          ┌───────────────────────────┴────────────────────────────┐
          ▼ yes                                                    ▼ no
  BASELINE genomes  ~5.7k [E]                             NEW genomes  18,145 [M]
  NCBI protein.faa  ~55M proteins [E]                     BFD funannotate proteins ~135-195M [E]
          │                                                        │
          │         build_funnel_input.py  (len ≥ 30 aa, ids BASE__… / NEW__…)
          └───────────────────────────┬────────────────────────────┘
                                      ▼
            ┌──────────────────────────────────────────────┐
  TIER 0    │ MMseqs2 linclust, baseline + new together    │   in : ~190-250M [E]
            │ id 95 %, cov 90 %                            │   out: ~85-110M reps [E]
            └──────────────────────┬───────────────────────┘
                                   ▼  representatives only
  TIER 1    │ linclust id 70 %, cov 80 % (bidirectional)   │   baseline-free clusters = novel@70
                                   ▼
  TIER 2    │ linclust id 50 %, cov 80 %                   │   baseline-free clusters = novel@50
                                   ▼
            delta_funnel_summary.py   (member-weighted; per genome, per tier)
                │                          │
                ▼                          ▼
      tier_summary.tsv            survivors_clusters.tsv   (baseline-free at 50 %)
      genome_tier.tsv                      │
                                           ▼
  STAGE 2  (not yet written) sensitive search of survivors vs baseline reps
            DIAMOND --more-sensitive, bidirectional-coverage rule, then profile
            search (Pfam / MMseqs2-profile).  Survivor count decides the cost.   [TBD]
                                           │
                                           ▼
  STAGE 3  what are they?   Pfam on survivors  →  known domain / new architecture / dark
            clade restriction, gene-neighbourhood, secretion signal
            structure (ESMFold → Foldseek) on a small, high-confidence family subset          [TBD]

  Parallel investigations (can start as soon as the baseline list is final):
  A. one Pfam family walk   (pfam.parquet, 666 MB, all proteomes already scanned)
  B. one CAZy family        (needs dbCAN at scale; current cazy parquet is a test file)
        both:  cluster domain-bearing proteins (~50 % id) → add genomes in order
               (baseline first, then new) → new subfamilies per genome vs random-order null
```

Everything in Tiers 0–2 is MMseqs2 CPU time; no DIAMOND until the survivor set is small.

## 5. Stage specifications

| Stage | Script | Key parameters | Output |
|---|---|---|---|
| Prior-annotation table | (ad hoc; `dg_prior_table.py`) | NCBI `annotation_info` + protein file ≥ 1 KB | `prior_annotation_status.tsv` |
| Gap download | `fetch_gap.sh` (NCBI_fungi `datasets` v18) | `--include gff3,protein,cds` | proteins added to `NCBI_ASM` |
| Classify BFD genomes | `bin/classify_genomes_by_prior_annotation.py` | match on ASMID, else `GC[AF]_n.n` prefix | `genome_classes.tsv`, `new_manifest.tsv` |
| Build funnel input | `bin/build_funnel_input.py` | min length 30 aa; drops duplicate ids | `funnel_input.faa`, `funnel_provenance.tsv` |
| Tier clustering | `launch/delta_funnel.sbatch` | linclust; tier 0: 0.95/0.90; later 0.80 cov, cov-mode 0; 64 cpu / 400 GB / 2 d | `tier_*_clusters.tsv`, `tier_*_reps.faa` |
| Summary | `bin/delta_funnel_summary.py` | DuckDB; known-flag propagated through nesting | `summary/` |
| Pangenomes (later) | `bin/species_pangenome.py`, `bin/compare_species_pangenomes.py` | see HANDOFF | species tables |

Run it: `sbatch launch/delta_funnel.sbatch` (defaults point at the full manifest and status table).

## 6. Assumptions (each can be wrong — see §9 for how we test them)

1. NCBI `annotation_info` presence means usable prior gene models. Submitter annotation counts, not only RefSeq.
2. BFD proteins of a new genome and NCBI proteins of a baseline genome are comparable. Different gene
   predictors give different models; some "novelty" will be gene-model difference, not biology.
3. Identity tiers are a sound proxy for "same family" at 95 / 70 / 50 %; below ~50 % linclust loses
   sensitivity, which is why survivors get a sensitive search.
4. A cluster with no baseline member at a tier is novel *relative to this baseline*, not absolutely.
5. Baseline choice defines the result. Baseline = genomes with prior annotation is a convenience
   split, not a claim that those genomes are phylogenetically representative.
6. Species-strain dedup in NCBI_fungi is acceptable (1,118 of 1,187 dropped annotated assemblies have an
   annotated sibling in the baseline; 57 strains are flagged because the kept sibling is unannotated).
7. Pfam hits from the existing annotation represent each proteome fairly.

## 7. Expected results and hypotheses

- Novel fraction falls steeply with lower identity tier: most new proteins have baseline relatives at 50 %.
- Novelty is **concentrated** — a minority of genomes (and clades) carries most of it. The 2000-genome pilot
  vs public nr-fungi: top 5 % of genomes hold 40 % of no-hit proteins; median genome 1.8 % no-hit [M].
- Many survivors carry known Pfam domains (new architectures), not new domains.
- Novelty per genome depends on taxonomic distance to the baseline: new genera/families add more than
  extra strains of a well-sampled species (testable with the genus/order breakdown).
- Within a single Pfam/CAZy family, new genomes add subfamilies at a rate above zero but below the
  per-protein novelty rate (diverse forms inside known families, not whole new families).

## 8. Measures of success

The project succeeds if each of these holds:

1. **Pipeline:** the funnel runs end to end on the full set within the stated budget (≤ 2 days, 400 GB) and the
   tier table reproduces on re-run.
2. **Self-consistency:** novel@95 ≥ novel@70 ≥ novel@50 for every genome (nesting never violated).
3. **Controls behave:**
   - *positive:* withholding a whole known order (`validate_positive_control.py`) makes that order's novelty
     jump; a baseline-vs-baseline split gives near-zero novelty at 95 %.
   - *null:* shuffled or reversed proteins give ~100 % novel at every tier.
4. **Robustness:** the headline fractions change by less than a stated margin when (a) low-quality genomes
   (BUSCO/OMArk bottom decile) are removed, (b) contaminant-flagged proteins are removed, (c) sibling-flagged
   strains are included/excluded.
5. **Interpretability:** survivors resolve into a reportable split — known domain / new architecture / dark —
   and families are called only with ≥ 3 genomes in ≥ 2 genera.
6. **Case study:** for one Pfam family, a delta-vs-genomes-added curve exists with a random-order null band,
   and we can say whether new genomes add subfamilies beyond the null.
7. **Decision value:** we can state, per clade, whether more sequencing/annotation of that clade is expected to add
   substantial new information (the practical purpose of delta gain).

## 9. Confounders, checks and known weaknesses

| Risk | Check |
|---|---|
| Self-hits: public nr already contains many query genomes (inflates "known" in the old design) | **Removed by design** (query genomes have no NCBI annotation), and the old 2000-run hits are being re-examined (§11 step 2b) |
| Gene-model differences, not biology | Compare novelty of BFD-annotated *baseline* genomes vs their own NCBI proteins (a calibration of annotation noise) |
| Contamination / over-prediction | Taxonomy of no-hit proteins; per-genome no-hit rate vs gene count and mean protein length; drop bottom BUSCO decile |
| Fragments scored as novel | Bidirectional coverage in the sensitive stage; weak hits in the old run were coverage-driven (median E ≈ 0, subject cov 50 %) |
| Rep-ownership bias (which genome "owns" a cluster rep) | Counting is per member protein; ownership never enters the metric |
| Uneven taxon sampling | Per-order/genus tables; permutation order tests |
| linclust at 50 % misses homologs | Survivor stage uses DIAMOND/profile search; report novelty at each tier, not just the last |
| Placeholder thresholds (min length 30, cov 0.9/0.8, tiers) | Sensitivity table at ±1 step |

## 10. Compute budget

| Step | Resources | Status |
|---|---|---|
| Funnel (3 linclust tiers on ~190-250M proteins) | 64 cpu, 400 GB, ≤ 2 d; reference: full BFD linclust 1h25m / 224 GB [M] | smoke test passed (80 genomes, 2m13s) [M]; **full run done 2026-11-06: job 31780270, 1h47m, 64 cpu** [M] |
| Pfam family walk | duckdb + mmseqs on a subset | not started |
| Survivor search | TBD from survivor count | not written |
| Structure prediction | TBD; only on a small family subset | not planned in detail |

## 11. Status log

- 2026-09/10: 100-genome and 2000-genome pilots of the universal DIAMOND design [M].
- 2026-11-04: full-set linclust done (78M reps); full DIAMOND search launched.
- 2026-11-05: strategy review (Fable) + comparison with the metagenomic-families paper
  (Aplakidou et al., NAR 2026): their funnel is Pfam/AntiFam filter → fast LAST → DIAMOND on the remainder →
  linclust 30 %/80 % families (≥ 25 members) → AlphaFold + Foldseek. Decision: **stop the universal search**
  (~15 % done; 5 days to the shutdown; mostly re-confirming known strong hits) and rebuild around an explicit baseline.
- 2026-11-05: prior-annotation status built for 27,401 assemblies; gap download done (192 of 202 accessions got
  proteins; 10 have none); baseline now 5,703 genomes. Funnel scripts written and smoke-tested.
- 2026-11-05: raw hits compressed and verified, 2000-genome re-bin and self-hit check done, both `work/`
  directories deleted.
- 2026-11-06: **first funnel run crashed in the summary stage** ("CSV Error ... Expected 3 columns, found 1"). Root cause was
  not quoting: 44 baseline folder names carried a trailing carriage return (the gap-download list was read from a CRLF CSV).
  MMseqs treats `\r` as whitespace and DuckDB as a newline, so those ids split. Fixed by stripping CR in the status table and
  writer, `build_funnel_input.py` and the summary parser, and renaming 44 directories / 220 files on disk. Rerun (31780270)
  gave base_genomes 5,703, base_missing_file 0.
- 2026-11-06: **funnel result [M]** (183,095,650 new proteins; "novel" = tier cluster has no baseline member):

  | Tier | Novel proteins | Fraction |
  |---|---:|---:|
  | 95 % id / 90 % cov | 91,857,334 | 50.2 % |
  | 70 % | 53,498,892 | 29.2 % |
  | 50 % | 41,042,274 | 22.4 % |

  Per genome (17,864; `summary/genome_tier.tsv`) the spread is wide: median fraction novel 35.6 % / 14.1 % / 10.2 % at
  95 / 70 / 50 %, 90th percentile 99.5 % / 82.1 % / 62.0 % (new-taxon genomes vs strains of baseline species).
  Nested clusters from stage 4: `summary/survivors_clusters.tsv` (1.0 GB).
- 2026-11-07: **predictor-mismatch control (E1) [M]** (`bin/predictor_control_prep.py`, `launch/predictor_control.sbatch`,
  `results/predictor_control/`): 300 baseline genomes annotated by NCBI (run A) vs by BFD (run B), each scored against the same
  reference of 1,000 other baseline genomes' NCBI proteins (the test genomes are not in the reference). Novel fraction A -> B:
  95 % 57.3 -> 60.0 %, 70 % 34.5 -> 35.8 %, 50 % 25.7 -> 26.8 % (protein counts 3.37 M vs 3.32 M; per-genome median change
  +0.3 / +0.8 / +1.0 points, 10th-90th percentile at 50 % -4.6 to +6.0 points; BFD higher in 64 % of genomes). **The NCBI-vs-BFD
  gene-model difference inflates funnel novelty by only ~1.5-3 points (mean) - it does not explain the headline 50 / 29 / 22 %.**
  Absolute levels are higher than the full funnel because the reference has 1,000 not 5,703 genomes.
- 2026-11-07: **rank decomposition of funnel novelty (E2) [M]** (`results/delta_funnel/summary/genome_stratum.tsv`): new genomes
  tagged by the lowest rank they share with a baseline genome. Novel at 50 %: strain of baseline species 13.3 % (10,673 genomes,
  60 % of new proteins, 36 % of novel proteins, median genome 5 %); new species in baseline genus 27.6 % (5,087 genomes, 35 % of
  novel proteins); new genus in baseline family 57.0 % (1,629; 22 %); new family 63.2 %, new order 59.6 %, new class or higher
  60.2 % (475 genomes together, 7 % of novel proteins). Novel proteins per genome rise from 1,374 (strain) to 2,817 (new species)
  to 5,652-6,351 (new genus/family). Species-collapsed fractions at 50 %: 15.5 / 34.5 / 59.1 / 63.0 / 60.9 / 70.4 %.
- 2026-11-07: **survivor clusters and their labels (F6) [M]** (`bin/survivor_annotate.py`, `results/survivor_annot/`). Novel-at-50 % survivors:
  17,095,910 clusters / 41.0 M proteins; 12.2 M clusters are seen in one genome only, 4.86 M in >= 2, 1.25 M in >= 5, 184 k in >= 20
  (9.4 M proteins). Labelled by the cluster representative only (BFD protein id joins directly to the BFD tables). **Most "novel at 50 %"
  clusters still carry a known Pfam domain**: fraction with a Pfam hit 41.6 % (1 genome), 51.2 % (2-4), 55.7 % (5-19), 35.9 % (>= 20);
  CAZy 1-2 %, MEROPS 1-2 %; unannotated 58 % / 49 % / 44 % / 64 %. Most frequent Pfam among survivors in >= 5 genomes: protein kinase
  PF00069, MFS PF07690, fungal Zn2Cys6 PF00172 / PF04082, cytochrome P450 PF00067 - i.e. divergent members of large known families. The >= 20-genome
  class has the largest unannotated share (64 %), consistent with strain-/population-specific or repeat-derived clusters; needs the
  contamination/TE triage (E5).
- 2026-11-07: **survivor sensitive search, 20,000-query test [M]** (`bin/survivor_split.py`, `bin/survivor_search_summary.py`,
  `launch/survivor_search.sbatch`): query = survivor representatives seen in >= 5 genomes (test: the 20,000 largest), target = 16.45 M tier-50
  representatives of baseline-containing clusters, MMseqs2 search -s 6.0, e <= 1e-5, query coverage >= 0.4. Result for the test
  (proteins): homolog (e <= 1e-10, identity >= 40 %) 75.6 %, close hit with low identity 4.5 %, distant (e 1e-10..1e-5) 1.9 %,
  **no hit 18.1 %**. So most large survivor clusters are detectable homologs of baseline-known clusters that merely fail the 50 % / 80 %
  coverage rule; roughly one cluster in six has no baseline homolog even by sensitive search. The test took 3.5 min on 32 CPUs;
  the full 1.25 M-query run (job 31828362) is queued. Test queries are the largest clusters, so the full-set split may differ.
- 2026-11-07: **taxonomic breadth of survivor clusters (E6) [M]** (`bin/survivor_taxa.py`, `results/survivor_taxa/`; survivors expanded through the
  nested tiers to their NEW proteins; 17,095,910 clusters / 41,042,274 proteins recovered exactly). **82 % of survivor clusters (14.06 M)
  and 61 % of survivor proteins (25.0 M) come from ONE species**: 12.24 M single-genome clusters (13.1 M proteins), 1.41 M / 3.6 M proteins in 2-4 genomes,
  306 k / 2.6 M proteins in 5-19 genomes, 110 k / 5.7 M proteins in >= 20 genomes of one species (strain-population or accessory content of heavily
  resequenced species, e.g. the S. cerevisiae and C. auris collections). Only 3.04 M clusters (16.0 M proteins) span >= 2 species; 0.90 M clusters
  (5.3 M proteins) span >= 2 genera (0.32 M clusters span >= 2 families). Among clusters seen in >= 5 genomes (1.25 M): 33 % one species,
  42 % one genus, 16 % one family, 4 % one order, 5 % multi-order. **So most novel-at-50 % protein mass is species-restricted; broadly distributed
  novel families are a small fraction.**
- 2026-11-07: **contamination / TE triage of survivors (E5) [M]** (`bin/survivor_triage.py`, `results/survivor_triage/`). Genome flags: 1,067 genomes
  with GC > 8 points from their order median and 322 with BUSCO < 50 (together 6-7 % of new genomes) supply **8.2 %** of the novel-at-50 %
  protein mass; GC outliers are clearly more novel only for "new species in a baseline genus" (frac novel 48.5 % vs 28.8 %), and the 59 genomes
  of new class or higher contribute 33 % of that tiny stratum's novelty from flagged genomes. Transposon/retroelement Pfam domains
  (210 Pfam families by name) label only ~2 % of Pfam-annotated survivor clusters and **~1 % of survivor proteins** (0.9-1.1 % by genome-count class).
  Caveats: GC is a crude contamination proxy; unannotated TE-derived genes are not counted. **Contamination and annotated TEs do not explain the survivor mass.**
- 2026-11-07: **species-collapsed / rarefied funnel novelty (E3) [M]** (from `genome_stratum.tsv`). The headline 50.2 / 29.2 / 22.4 % is
  protein-weighted, so heavily resequenced species (60 % of new proteins are strains of baseline species) dominate it. Rarefying to k random genomes
  per species (50 draws, protein-weighted) gives 95 / 70 / 50 %: k = 1 76.9 / 49.9 / 37.6; k = 2 74.6 / 47.5 / 35.8; k = 5 70.5 / 44.0 / 33.1; k = 10 67.1 / 41.3 / 31.2;
  all genomes 50.2 / 29.2 / 22.4. One value per species (median over its genomes): 5,454 species, mean novel 78.7 / 52.4 / 39.4 %, median at 50 % 37 %;
  only the 4,632 species absent from the baseline: 85.6 / 58.3 / 43.7 %. **So the all-genome headline understates novelty per species and per lineage
  (resequenced strains add little); the per-species view (about 40 % novel at 50 % identity) is the better measure of how much sequence diversity the new
  genomes add, and is dominated by new species.**
- 2026-11-07: **full survivor sensitive search [M]** (job 31828367, highclock, 1 h 25 min; `results/survivor_search/`). Queries = the 1,254,423 survivor
  representatives seen in >= 5 genomes (18.4 M proteins); targets = 16.45 M representatives of baseline-containing 50 % clusters; MMseqs2 -s 6.0, e <= 1e-5, query coverage >= 0.4.
  Clusters (proteins): **homolog of a baseline-known cluster (e <= 1e-10, identity >= 40 %) 79.3 % (77.6 %)**; close hit with identity < 40 % 8.5 % (7.3 %); distant (e 1e-10..1e-5)
  1.6 % (1.8 %); **no hit 10.7 % = 133,599 clusters / 2.45 M proteins (13.3 %)**. Pfam labels: homolog class 61 %, low-identity 44 %, distant 9 %, **no-hit class 4.8 %**;
  CAZy/MEROPS in the no-hit class 0.2 %. No-hit clusters are mostly narrow: 62 % one species (mean 2.2 species), 5 % span >= 2 genera, 1 % >= 2 orders.
  **Conclusion: of the 41.0 M proteins that are novel at 50 % identity / 80 % coverage, the ones in clusters seen in >= 5 genomes are overwhelmingly
  detectable homologs of baseline families; only ~2.5 M proteins (about 6 % of survivor proteins, 1.3 % of all new proteins) in 134 k clusters have no baseline homolog even
  by sensitive search, and these are mostly small, lineage-restricted and unannotated (candidate orphan / lineage-specific genes, or gene-model artefacts - not yet
  triaged). Clusters in 1-4 genomes (15.8 M clusters, 22.6 M proteins) were not searched.**
- Profile search of the no-hit survivor clusters (2026-11-07): the 133,599 no-hit clusters (2,453,840 proteins) were re-clustered at 30 % identity into 107,819
  profile clusters, each turned into a profile and searched (-s 7.0, e <= 1e-3, qcov >= 0.4) against baseline-known reps. **10,976 profile clusters (10.2 %), 315,862
  proteins (12.9 %) gain a baseline hit; 96,843 clusters (89.8 %), 2.14 M proteins remain without any detectable baseline homolog.**
- Singleton sample (1 M survivor clusters seen in 1 genome, sensitive search): 63.0 % homolog, 12.8 % close low-identity, 2.9 % distant, 21.2 % no hit (by clusters).
  Singletons are therefore about twice as likely as >= 5-genome clusters (10.7 %) to lack any homolog, consistent with a mix of real orphans and gene-model artefacts.
- Lifestyle correlates (2026-11-07; `bin/lifestyle_novelty.py`, FungalTraits 1.2 genus table, Polme et al. 2020, `refs/lifestyle/`; results `results/lifestyle*`).
  1,061 of 1,204 genera with a new genome (16,085 of 17,865 genomes) have a lifestyle; unit = genus; contrast = genus mean minus its order mean (within-order),
  permutation test shuffling labels among genera inside each order (2,000 draws). Coarse genus-level labels, so only a screen.
  All new genomes, within-order deviation in the novel fraction at 50 % identity: mycorrhizal +4.4 points (p 0.0005), saprotroph +1.2 (p 0.016), plant pathogen
  -3.9 (p 0.0005), animal-associated -2.8 (n.s.); yeast growth form -4.1 vs filamentous +4.4. These shrink or vanish when relatedness is held fixed:
  restricted to new species in a baseline genus, mycorrhizal stays +7.5 (p 0.003) while plant pathogen (-0.7), saprotroph, animal and yeast are about 0; restricted to
  new genus or higher, plant pathogens stay -3.0 (p 0.0045) and mycorrhizal +2.0 (p 0.012). Reading: most lifestyle differences are sampling depth / relatedness
  (plant pathogens and yeasts have deep baseline coverage); the robust signal is that **mycorrhizal (mostly Agaricomycetes) genomes are more novel than their
  order mates at every relatedness level**, and plant pathogens are somewhat less novel among taxonomically distant genomes. Not causal: lifestyle is nested in clade.
- Sensitive search of survivors in 2-4 genomes (job 31830779, 5 h; 3,602,662 clusters of the 17.1 M searched; queries = clusters with 2-4 genomes, same -s 6.0 settings):
  by clusters 74.2 % homolog, 11.1 % close low-identity, 1.9 % distant, **12.8 % no hit (460,890 clusters, 1.22 M proteins)**. Together with the >= 5-genome (10.7 %) and
  singleton-sample (21.2 %) results, the no-hit share rises as genome count falls (10.7 -> 12.8 -> 21.2 %), as expected if no-hit clusters are enriched for gene-model
  artefacts and recent lineage-specific genes. (Note: 3.6 M queries rather than the 5.08 M planned: the 2-4 genome clusters are fewer than the earlier estimate.)
- Alignment-based catalytic-residue check (`bin/catalytic_align.py`, hmmalign to Pfam PF00704 / PF00082 profile; catalytic sites located on the HMM consensus; c70 clusters):
  GH18 (sites D1, D4, D6, E8): regions covering all sites / all canonical / canonical given coverage = baseline 97.6 / 83.5 / 85.5 %; new-known 98.2 / 84.2 / 85.8 %;
  new-novel 93.1 / 74.4 / 79.8 % (65,500 / 175,122 / 13,169 regions). S08 (Asp, His, Ser): baseline 90.9 / 66.8 / 73.5 %; new-known 93.6 / 70.3 / 75.1 %; new-novel
  70.2 / 44.8 / 63.9 % (64,787 / 183,064 / 17,126). Asp and Ser are canonical in 98-99.8 % of every group; the His position is only 74-75 % canonical even in
  baseline (the site is variable across S8 subfamilies or the HMM column is imperfect), so use the His column only for relative comparison (novel 65 % vs 74-75 %).
  Reading: novel regions are far more often truncated (S08: 30 % lack a site vs 6-9 %), which explains most of the intact-motif deficit; among spanning regions
  novel ones are 6 (GH18) to 10 (S08) points less often fully canonical. Novel S08 hits are therefore mostly fragments or divergent subtilisin-like proteins,
  consistent with the earlier regex result and the likely false-hit S08 on carboxylases; they are not evidence of new catalytically active enzymes.
- BGC-backbone proxy (user-approved substitute for antiSMASH, which cannot finish for 23.5 k genomes before the shutdown; `launch/bgc_run.sh`, results `results/family_*_bgc`):
  10 secondary-metabolism backbone / tailoring Pfams (ketosynthase PF00109, AMP-binding PF00501, condensation PF00668, acyltransferase PF00698, thioesterase PF00975,
  UbiA prenyltransferase PF01040, terpene synthase C PF03936 and PF19086, ketoreductase PF08659, dehydratase PF14765; PF01397 had no regions) through F1-F3 (regions
  from 0.01-2.0 M per family). Baseline-vs-baseline null with matched genome counts (307 and 1,010 genomes per side): the novel fraction of new genomes minus the null
  at c70 is -1.3 to +2.3 points, every 95 % interval spans 0 except ketoreductase PF08659 (+0.9 points, [0.0, 1.8], at 1,010 genomes; +0.2 at c50). **New genomes add no more
  backbone-enzyme diversity than the same number of baseline genomes**, matching the CAZy/MEROPS/random-Pfam pattern; secondary-metabolism backbones are not a hotspot of
  new protein-level diversity. Limits: Pfam backbone domains only (no cluster context, no product prediction); terpene synthases (PF03936) have wide intervals (1.4 k baseline genomes).
- HGT proxy (user-approved substitute for a phylogeny-based screen; `bin/hgt_proxy.py`, NCBI taxdump, results `results/hgt_proxy`): best Swiss-Prot hit (e <= 1e-5, qcov >= 0.5)
  of survivor representatives, grouped by taxonomy of the hit organism, metagenome-bin genomes excluded. Prokaryote best hit: survivors in 1 genome 4.8 % (2.99 M with hit),
  2-4 genomes 4.2 % (1.18 M), >= 5 genomes 4.1 % (0.45 M) vs 5.2 % for 2.72 M background proteins (hash sample, non-survivors); >= 40 % identity 1.0-1.3 % vs 1.3 %. Fungal best hit 75-79 %.
  **No enrichment of prokaryote-like hits among novel proteins** (no evidence of a large HGT / contamination pool). Limits: only proteins with a Swiss-Prot hit (< 50 % of
  survivors), no phylogeny, so individual transfers are not tested.
- Status 2026-11-07 (goal run): done: funnel, F0-F5 (717 families + subfamilies + random control), predictor control, rank decomposition, rarefaction,
  survivor labels / breadth / triage / sensitive search (>= 5 genomes, 2-4 genomes, singleton sample), profile search of no-hit clusters, F3c class-matched null,
  split-sample test, alignment-based catalytic analysis (GH18, S08), lifestyle correlates. Not done: structure prediction (cannot fit before the 2026-11-10 shutdown),
  antiSMASH itself (replaced by the BGC-backbone proxy), phylogeny-based HGT (replaced by the Swiss-Prot taxonomy proxy), per-order F3c ranking, search of the remaining 1-genome clusters beyond the 1 M sample.

## 12. What the old design showed and why it changed

The first design searched every representative (78M) against NCBI nr restricted to fungi (`nr_cluster_seq.dmnd`,
470M-sequence clustered nr, `--taxonlist 4751` at search time, so the full 195 GB index is still scanned) and
UniProt-fungi, then binned hits strong / weak / no-hit. 2000-genome result [M]:

| | strong | weak | no hit |
|---|---:|---:|---:|
| vs nr-fungi | 78.0 % | 17.6 % | 4.35 % |
| overall (nr or UniProt) | 80.0 % | 15.7 % | 4.31 % |

Re-examination of the raw hits (2026-11-05, 10.6M queries, `reports/rebin_selfhit_n2000_per_genome.tsv`) [M]:

- The old rule is reproduced exactly (37 of 10.6M queries differ, a tie-break).
- **Best-class re-bin:** 731,225 of the 1,870,422 "weak" queries (39 %) have another row that meets the strong rule.
  Strong rises 78.0 % → 84.9 %, weak falls 17.6 % → 10.7 %; no-hit is unchanged (4.35 %). The single-row rule was
  misclassifying a large part of "weak".
- **Self-hit evidence:** genomes that already have NCBI annotation (488 in the sample; their proteins are in nr)
  have 47.5 % of queries with a 99–100 % identity top hit, vs 12.6 % for genomes without prior annotation (1,502).
  No-hit is 2.6 % vs 5.2 %. Part of this gap is taxon composition, but it supports excluding prior-annotated genomes
  from the query side, as the new design does.

Problems: (1) it measures "not yet public", not delta between genome sets; (2) public nr contains many query
genomes, so strong hits can be self-matches; (3) weak hits were mostly coverage-driven, not distant homology;
(4) counts were per representative; (5) ~2,600 DIAMOND jobs to re-find the 80 % that match anyway. The
2000-genome tables are kept as a pilot; the full-set clustering (`results/protein_novelty_full`) is reused.

## 13. Open questions

- Which genomes are "new" when the same strain has several assemblies? (57 flagged strains now.)
- Is a date-ordered addition curve possible? Assembly release dates exist in the NCBI snapshot; annotation dates
  exist for baseline genomes only.
- Tier set: add 30 %? Needs a sensitivity test (linclust is weak there).
- Include UniProt for functional labels on survivors only?
- Which Pfam family and which CAZy family for the case studies?

## 14. Where the numbers come from

- `DeltaGain_Fungi/prior_annotation_status.tsv` — per-assembly status table.
- `DeltaGain_Fungi/results/protein_novelty_n2000/` — pilot tables; `…/protein_novelty_full/` — full clustering.
- `Fungi_BFD_runs/tables/{pfam,species,asm_stats}.parquet` — Pfam hits, taxonomy, assembly stats.
- `NCBI_fungi/lib/ncbi_accessions.json` — NCBI datasets snapshot (annotation info, assembly dates).
- Pangenome and accumulation scripts, their tests and the earlier decision history: `DeltaGain_Fungi/HANDOFF.md`.
