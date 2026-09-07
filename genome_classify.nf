#!/usr/bin/env nextflow
/*
 * genome_classify.nf — DeltaGain Stage 1: label-free genome classification.
 *
 * Answers "which of the 22k BFD genomes are redundant strains of each other?"
 * WITHOUT bucketing by pre-existing species name first (the chicken-and-egg
 * problem a species-name-bucketed approach has for genuinely novel/mislabeled
 * species -- see /bigdata/stajichlab/shared/projects/BFD/Ideas/protein_novelty_esm2_plan.md).
 *
 * BFD-internal only: this does NOT compare against UniProt/RefSeq reference
 * genomes (see MASH_PREFILTER_GLOBAL's header comment for why). "Is this
 * species already known elsewhere" is answered separately, at the taxonomy
 * level, not by genome ANI.
 *
 * Pipeline:
 *   1. MASH_SKETCH          per genome (input_clean_genomes/<ASMID>.fa.gz)
 *   2. MASH_PREFILTER_GLOBAL  one all-vs-all mash dist across ALL genomes
 *   3. FIND_COMPONENTS      label-free union-find clustering at prefilter_ani
 *   4. SKANI_TRIANGLE       exact ANI refinement, per multi-member component
 *   5. CLASSIFY_GENOMES     combine + classify (see bin/classify_genomes.py
 *                           docstring for the cluster_class definitions)
 */

nextflow.enable.dsl = 2

include { MASH_SKETCH }           from './modules/mash/MASH_SKETCH/main.nf'
include { MASH_PREFILTER_GLOBAL } from './modules/mash/MASH_PREFILTER_GLOBAL/main.nf'
include { FIND_COMPONENTS }       from './modules/components/FIND_COMPONENTS/main.nf'
include { SKANI_TRIANGLE }        from './modules/skani/SKANI_TRIANGLE/main.nf'
include { CLASSIFY_GENOMES }      from './modules/classify/CLASSIFY_GENOMES/main.nf'

workflow {

    samples_csv = file(params.samples, checkIfExists: true)

    // Reuses Fungi_BFD_runs/samples.csv as-is (ASMID, SPECIES, ... columns) --
    // not reinvented here. Genome files resolved via params.genome_dir using
    // the existing <ASMID>.fa.gz convention (Fungi_BFD_runs/input_clean_genomes).
    // NOTE on `.take(n)`: confirmed on nextflow/26.04.6 that the reassignment
    // form `ch = ch.take(n)` (the pattern in the nextflow-hpcc skill's own
    // template) fails to parse ("Missing process or function take(...)"),
    // while `.take(n)` chained directly inline works fine. Root cause not
    // fully pinned down (DSL2 reserves `take:` as a workflow-section keyword
    // and appears to mis-scope the reassigned form specifically) — chaining
    // inline sidesteps it either way. `.take(-1)` is Nextflow's documented
    // "take everything" sentinel, used here for the n_test=0 (no limit) case.
    genomes_ch = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def asmid = row.ASMID.trim()
            def fa    = file("${params.genome_dir}/${asmid}.fa.gz", glob: false)
            tuple(asmid, fa)
        }
        .filter { asmid, fa -> fa.exists() }
        .take(params.n_test as int > 0 ? (params.n_test as int) : -1)

    MASH_SKETCH(genomes_ch)

    MASH_PREFILTER_GLOBAL(
        MASH_SKETCH.out.map { asmid, msh -> msh }.collect()
    )

    FIND_COMPONENTS(MASH_PREFILTER_GLOBAL.out, samples_csv)

    // Fan out: one SKANI_TRIANGLE task per multi-member component.
    component_genomes_ch = FIND_COMPONENTS.out.component_genome_lists
        .flatten()
        .map { f -> tuple(f.baseName.tokenize('.')[0], f) }

    SKANI_TRIANGLE(component_genomes_ch)

    CLASSIFY_GENOMES(
        FIND_COMPONENTS.out.components,
        SKANI_TRIANGLE.out.map { comp_id, tsv -> tsv }.collect().ifEmpty([])
    )
}
