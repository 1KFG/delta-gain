# DeltaGain

**A Framework for Quantifying the Annotation Yield of Large-Scale Fungal Genome Sequencing**

## What this is

DeltaGain asks: when we annotate genomes that were previously unannotated (the
bulk of the ~22k-genome Fungi_BFD set), how much genuinely new information does
that bring in — beyond what UniProt (fungi, 2026_02) and RefSeq fungal
proteomes already give us?

The focus is not "new genomes exist" but **what annotating them adds**: new
domains, new domain architectures, new members of known gene families,
taxonomic/phylogenetic breadth not otherwise sampled, and lineage-specific
expansions of known gene families (secondary metabolites, effectors,
CAZymes, etc.) — five separate axes, not one novelty number. See the design
doc for the full rationale.

## This is the pipeline code repo — not where results/reports live

As of 2026-09-07, this project split into two repos, the way an nf-core-style
pipeline separates from the project instances that run it:

- **`DeltaGain`** (this repo, `1KFG/delta-gain`) — the reusable pipeline
  itself: `nextflow.config`, the entry `.nf` scripts, `modules/`, `lib/`,
  `bin/`, `conf/`, and `tests/` (stub-validation data for the pipeline logic,
  not project results).
- **[`DeltaGain_Fungi`](https://github.com/1KFG/delta-gain-fungi)**
  (`1KFG/delta-gain-fungi`) — the project instance that runs this pipeline
  against the Fungi_BFD dataset specifically. `docs/` (the results website,
  served at [dg.fungalgenomes.org](https://dg.fungalgenomes.org)),
  `reports/` (Quarto report source), and eventually `results/`/project-
  specific params all live there, not here.

**Implication for `conf/`:** the current profile configs still hardcode
Fungi_BFD-specific paths (`Fungi_BFD_runs/samples.csv`, `input_clean_genomes/`,
`nr_cluster_seq.dmnd`, etc.) directly in this repo's tracked config — a
holdover from before the split. The intended end state is for this repo to
carry generic parameter *names* with placeholder/no defaults, and for
`DeltaGain_Fungi` to supply the real paths via its own `-params-file`. That
refactor hasn't happened yet; treat the current hardcoded defaults as
Fungi_BFD-specific until it does.

`results/` and `work/` are also still here (see `.gitignore`'s note) rather
than in `DeltaGain_Fungi`, because a real pipeline run was actively writing to
them at split time — move them over once that run completes, not mid-run.

## Design doc

`/bigdata/stajichlab/shared/projects/BFD/Ideas/protein_novelty_esm2_plan.md`
— living plan, includes the Fable review and revisions.

## Two separate reference needs — don't conflate them

Working through the architecture surfaced an important distinction (see
git history / conversation for the full reasoning):

1. **Genome-level classification** (Stage 1, built): "which BFD genomes are
   redundant strains of each other?" — this is **BFD-internal only**. It does
   NOT need UniProt/RefSeq reference genomes mixed in — whether a species is
   already known elsewhere is a taxonomy lookup, not a genome-ANI question.
2. **Protein-level novelty screen** (Stage 3, partially scaffolded): "has this
   protein sequence been seen before?" — pure sequence search against a
   protein reference database. Genome/assembly provenance is irrelevant here.

An earlier draft of this pipeline mixed reference genomes into Stage 1's skani
clustering — removed as unnecessary once this distinction was clear.

## Stage 1 — genome classification (built, stub-validated)

Label-free: genomes are clustered by ANI distance, not bucketed by their
pre-existing species name first. This matters because bucketing by label has
a chicken-and-egg problem for the exact case this project cares about — a
genuinely novel or mislabeled species either becomes an outlier in the wrong
bucket, or (for genomes with no clean species match at all) never gets
compared to anything. The existing `Fungi_BFD_runs/results/ANI/skani/`
SPECIES/GENUS outputs are label-bucketed in exactly this way; Stage 1 here is
deliberately label-free instead.

Pipeline (`genome_classify.nf`, profile `genome_classify`):

1. **MASH_SKETCH** — per genome (`Fungi_BFD_runs/input_clean_genomes/<ASMID>.fa.gz`)
2. **MASH_PREFILTER_GLOBAL** — ONE all-vs-all mash distance across every BFD
   genome (~22k), no taxonomic bucketing
3. **FIND_COMPONENTS** (`bin/find_components.py`) — union-find connected
   components at `params.prefilter_ani` (default 80%); a component spanning
   more than one existing species label is exactly the kind of disagreement
   this step is meant to surface
4. **SKANI_TRIANGLE** — exact-ANI refinement, fanned out one task per
   multi-member component only (never all-vs-all at the expensive-ANI level)
5. **CLASSIFY_GENOMES** (`bin/classify_genomes.py`) — combines the above into
   `genome_classification.tsv` with `cluster_class`: `redundant_strain`,
   `novel_strain_same_species`, `label_mismatch`, `unresolved_below_ani_floor`,
   `singleton_isolated`. See that script's docstring for exact definitions.

Reuses `Fungi_BFD_runs/samples.csv` and `input_clean_genomes/` directly (not
reinvented) — see `conf/profile_genome_classify.config`.

**Provisioning:** Apptainer, using the same pre-built sandbox images already
proven working in `Fungi_BFD/nextflow/conf/profile_ANI.config` (mash, skani —
that config's comments document a squashfuse-timeout failure mode under
module-based provisioning at this pipeline's concurrency, fixed by switching
to sandboxes; reused here rather than re-solving it).

**Validated:** `nextflow run genome_classify.nf -profile test -stub-run --n_test 2`
passes end-to-end (tiny synthetic genomes in `tests/data/`).

**Known Nextflow gotcha hit while building this** (nextflow/26.04.6): the
channel operator `.take(n)` fails to parse ("Missing process or function
take(...)") when used via reassignment (`ch = ch.take(n)`) — including the
exact pattern in the nextflow-hpcc skill's own template — but works fine
chained inline onto a fresh channel definition. `genome_classify.nf` uses the
inline form; worth flagging if the skill template gets reused elsewhere.

**Not yet run for real** — only stub-validated. Before a real run: confirm
`params.singularity_cache` / `NXF_SINGULARITY_CACHEDIR` is set, and note the
fixed 90%/95% ANI thresholds are a known simplification (fungal species
boundaries vary by lineage; see design doc) worth revisiting once real
results come back.

## Stage 3 — protein novelty screen (partially scaffolded, not wired into a workflow yet)

Quick-and-dirty v0, converged on during architecture discussion: DIAMOND
blastp against the **already-staged** `nr_cluster_seq.dmnd`
(`/srv/projects/db/ncbi/diamond/20260128/`, 470M sequences, built with
`--taxonmap`/`--taxonnodes`/`--taxonnames`), restricted to Fungi at search
time via `--taxonlist 4751` — no separate fungi-only database needs to be
built or hosted; the taxon restriction is a single flag on the search command
itself.

Considered and rejected: pulling a fresh nr slice or building one from
scratch — unnecessary, since a taxonomy-aware nr-derived DIAMOND db is already
staged. Also considered: UniProt fungi via REST stream
(`https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:4751&format=fasta`)
and UniProt Reference Proteomes/Fungi (curated, one proteome per species) —
still worth pulling as a complementary curated baseline for axis 4
(taxonomic breadth), just not needed for the nr-based novelty screen itself.

**UniProt fungi reference set — two paths, both in `protein_novelty.nf`:**

1. **`FETCH_UNIPROT_FUNGI`** (REST `/stream`) — small connectivity/content
   smoke test only. **Live-tested 2026-09-06**: produced 20 genuine Swiss-Prot
   fungi records end-to-end (`--uniprot_test_records 20`, the profile's safety
   default). Two real API constraints discovered and now documented in the
   module: `/stream` has no server-side small-N option (`size=` 403s even at
   1), and it has a hard server-side result-count cap that the *unrestricted*
   (Swiss-Prot+TrEMBL) fungi query exceeds outright — 403 "too many results...
   fetch batch by batch," rejected before any data streams, so client-side
   truncation can't work around it. Only `reviewed:true` (Swiss-Prot, 38,451
   entries) is wired here; kept as the cheap smoke-test path, not the real
   comprehensive pull.

2. **`FETCH_UNIPROT_FUNGI_FTP` → `CONVERT_DAT_TO_FASTA` → `COMBINE_UNIPROT_FASTA`
   → `MMSEQS_LINCLUST`** — the real, comprehensive path, built after
   confirming (via REST total-count headers) that Swiss-Prot (38,451) +
   TrEMBL (17,534,069) = 17,572,520, matching the unrestricted total exactly —
   i.e. Swiss-Prot alone is only ~0.2% of UniProtKB fungi; TrEMBL is
   necessary for comprehensiveness. Bulk-downloads UniProt's pre-split
   taxonomic-division files from FTP instead of paginating the REST API
   (**live-verified 2026-09-06 that both files still exist**:
   `uniprot_sprot_fungi.dat.gz` ~60MB, `uniprot_trembl_fungi.dat.gz` ~13GB) —
   no query size cap, no pagination needed, one bulk download per release.
   - `.dat` flatfiles (not FASTA — UniProt doesn't publish a FASTA division)
     carry full InterPro/Pfam/GO cross-references (`DR` lines) for **every**
     entry, including TrEMBL (from UniProt's own automatic InterProScan-based
     annotation pipeline). `CONVERT_DAT_TO_FASTA` extracts both FASTA *and*
     these cross-refs (`bin/dat_to_fasta_and_annotations.py`, via
     `Bio.SwissProt.parse`) rather than discarding them — this gives Stage 5
     (domain/fold axis classification) pre-computed annotations for any
     UniProt sequence a BFD protein hits, for free.
   - **Real bug found and fixed while validating this**: the default
     `python3`'s Biopython (1.79) throws `IndexError` parsing current
     (2026-release) UniProt flatfile syntax in `Bio.SwissProt.parse`'s
     feature-table reader. Biopython 1.87 (also available as a module) parses
     the identical real record correctly — confirmed live against an actual
     downloaded sample. The module pins `module load biopython/1.87` inline,
     not the default.
   - `MMSEQS_LINCLUST` dereplicates the combined ~17.5M-sequence set (linclust,
     not cascaded `cluster` — same scaling reasoning as the design doc's
     Step 2). Uses the already-cached
     `mmseqs2-17.b804f--hd6d6fdc_1.img` Apptainer image — nothing new to build.

**Both fetch modules and the full chain are `-stub-run` validated** (7/7
processes succeed end-to-end: FETCH_UNIPROT_FUNGI, both
FETCH_UNIPROT_FUNGI_FTP divisions, both CONVERT_DAT_TO_FASTA divisions,
COMBINE_UNIPROT_FASTA, MMSEQS_LINCLUST). **Not yet run for real past small
tests** — the FTP path is a genuine 13GB download + a long, resource-heavy
linclust job on ~17.5M sequences (`epyc`, 32 cpus, 256GB, 48h in the profile
config) and hasn't been submitted yet; confirm before doing so.

**BFD-side extraction and novelty screen — wired 2026-09-08 (HANDOFF.md
items 1/3/4):**

- **`BUILD_PEP_MANIFEST`** (`bin/build_pep_manifest.py`) maps each
  samples.csv ASMID to its `Fungi_BFD_runs/input/pep/<tag>.proteins.fa`
  file. The pep filenames are NOT ASMID-keyed — they're named by
  `makeSampleTag(SPECIES, STRAIN)`, the exact same species+strain
  sanitization convention `Fungi_BFD/nextflow/modules/common/utils.nf`
  uses, ported to Python here rather than reimplemented ad hoc. As of
  2026-09-08: 10,922/23,683 samples.csv rows (46%) have a matching pep file
  (the annotation run is still in progress — most of the ~22k genomes
  aren't annotated yet); 6 tag collisions exist where two ASMIDs map to the
  same species+strain tag (e.g. two `Fusarium graminearum PH-1` accessions)
  — those are flagged `ambiguous_tag` rather than guessed at, since the
  filename alone can't disambiguate which ASMID's proteins are in the file.
  `params.n_test` caps how many `matched` rows get used downstream (see the
  profile config's SAFETY DEFAULT — 100, for an initial pilot, not the full
  annotated set yet).
- **`COMBINE_BFD_PROTEINS`** (`bin/combine_bfd_proteins.py`) concatenates
  the matched pep FASTAs, prefixing every sequence ID with its ASMID
  (`<ASMID>__<original_id>`) so provenance survives dereplication —
  otherwise funannotate's per-genome locus tags collide across genomes.
- **`MMSEQS_LINCLUST`** (generalized, see its header comment) dereplicates
  BFD's own combined protein set the same way it dereplicates UniProt-fungi
  — reused via a DSL2 `include ... as` alias
  (`MMSEQS_LINCLUST_BFD`/`MMSEQS_LINCLUST_UNIPROT`), not duplicated.
- **`SPLIT_FASTA`** (seqkit) chunks the dereplicated representative set for
  parallel DIAMOND search.
- **`DIAMOND_NOVELTY_SCREEN`** (vs. `nr_cluster_seq.dmnd`, `--taxonlist`)
  and the new **`DIAMOND_MAKEDB` → `DIAMOND_SEARCH_UNIPROT`** (vs. this
  pipeline's own dereplicated UniProt-fungi set — no `--taxonlist` needed,
  already fungi-only by construction) both run per chunk.
- **`BIN_NOVELTY_HITS`** (`bin/bin_novelty_hits.py`) combines both targets'
  hits per query into `strong_hit`/`weak_hit`/`no_hit`, E-value + coverage
  per the design doc's framing fix. **Threshold defaults are placeholders**
  (reuses `params.diamond_evalue`'s existing 1e-5 for "strong", 70%
  bidirectional coverage as a conventional rule of thumb) — not yet backed
  by this dataset's real hit-score distribution; revisit once a real pilot
  run's output is in hand.

**Full `protein_novelty.nf` graph (15 processes) is `-stub-run` validated
end-to-end.** Still stub-only past this point, though: the FTP/linclust
UniProt path (see above) and the BFD-side chain both need a real run before
any of their numbers mean anything — start with the pilot scale
(`params.n_test`), not the full annotated set.

**Not yet built:** the Stage 3b controls (artefact/contamination/
fast-evolving-family filters) and Stage 5 domain/fold sorting from the
design doc — both still design-doc-only.

## Layout

- `lib/utils.nf` — Groovy resource-scaling helpers (mash/skani cpu/mem/time by
  group size), sized for this project's actual scale rather than copied
  verbatim from the sibling project's per-species-group tiers
- `modules/` — one process per tool step, DSL2 `main.nf` per subdirectory
- `bin/` — Python helpers (`find_components.py`, `classify_genomes.py`)
- `conf/` — per-stage profile configs (`profile_genome_classify.config`) +
  `test.config` for `-stub-run`
- `tests/data/` — tiny synthetic samples.csv + genomes for stub validation
- `results/`, `logs/`, `work/` (gitignored) — created on first real run;
  temporarily still here rather than in `DeltaGain_Fungi` (see the repo-split
  section above) because a run was in flight when the split happened
