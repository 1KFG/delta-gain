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
| annotated, protein file on disk (`baseline_ok`) | 5,511 | **baseline** |
| annotated, folder present, protein file missing | 158 | being fetched (job `fetch_gap`) → baseline |
| annotated, never downloaded | 1,228 | 44 fetched; 1,184 are same-strain duplicates dropped by the project's dedup (1,118 already have an annotated sibling in the baseline) |
| no prior annotation, downloaded | 18,145 | **new** |
| no prior annotation, not downloaded | 2,359 | not in our set |

Protein counts (all **[E]** until the funnel input is built; per-genome averages from a 40+40 genome smoke test):

| Set | Genomes | Proteins/genome | Proteins |
|---|---:|---:|---:|
| Baseline (NCBI proteins) | ~5.7k | ~9.6k [M, 40 genomes] | ~55M [E] |
| Query (BFD, new genomes) | ~18.1k | 7.5k–10.7k [E] | ~135–195M [E] |
| Funnel input (both) | ~23.8k | | ~190–250M [E] |
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
| Funnel (3 linclust tiers on ~190-250M proteins) | 64 cpu, 400 GB, ≤ 2 d; reference: full BFD linclust 1h25m / 224 GB [M] | smoke test passed (80 genomes, 2m13s) [M]; full run pending baseline download |
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
- 2026-11-05: prior-annotation status built for 27,401 assemblies; gap download (1,386 accessions) submitted;
  funnel scripts written and smoke-tested.
- In progress: **2a** compress old raw hits; **2b** self-hit / best-class re-bin on the 2000-genome hits; **2c** delete
  the two `work/` directories after verification.
- Next: finish gap download → refresh status → launch the full funnel → fill §3 with measured counts →
  write the survivor stage → Pfam family walk.

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
