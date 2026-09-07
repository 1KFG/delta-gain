// Exact-ANI refinement within one mash-prefiltered connected component.
// Components are typically small (one candidate species/strain cluster), so
// this is the expensive-but-accurate step, run only where mash already found
// enough similarity to be worth refining -- never all-vs-all across all 22k.

include { skaniCpusFor; skaniMemoryFor; skaniTimeFor } from '../../../lib/utils.nf'

process SKANI_TRIANGLE {
    tag    { comp_id }
    label  'skani'

    // n_genomes passed as an explicit val(), not derived by reading
    // genome_list inside these closures -- same convention already
    // established in Fungi_BFD/nextflow/modules/ani/compare/SKANI_COMPARE
    // ("passed as an explicit val(n_genomes) channel element so the resource
    // directive closure can use it reliably at task-submission time").
    // Confirmed live 2026-09-07 why the alternative doesn't work: calling
    // .countLines() (a Nextflow-only Path extension) on the process input
    // inside a resource closure throws
    // java.nio.file.ProviderMismatchException, and even forcing
    // .toFile().readLines().size() then fails with "No such file or
    // directory" -- resource directives are evaluated before the task's own
    // file staging happens, so reading file content from an input path at
    // that point is unreliable regardless of API used. Computing the count
    // in Groovy at channel-construction time (see genome_classify.nf) and
    // passing it as a plain value sidesteps the whole class of problem.
    cpus   { skaniCpusFor(n_genomes) }
    memory { skaniMemoryFor(n_genomes, task.attempt) }
    time   { skaniTimeFor(n_genomes, task.attempt) }

    input:
        tuple val(comp_id), val(n_genomes), path(genome_list)

    output:
        tuple val(comp_id), path("${comp_id}.skani.tsv")

    script:
    """
    cut -f2 ${genome_list} > ${comp_id}.paths.txt
    # --sparse: pairwise (query, ref, ANI%%) triplets, not a phylip matrix --
    # matches the format already produced by the existing Fungi_BFD ANI
    # pipeline (see Fungi_BFD_runs/results/ANI/skani/SPECIES/*/*.ani.tsv) and
    # what bin/classify_genomes.py's parser expects.
    skani triangle --sparse ${params.skani_preset_flag} \\
        -l ${comp_id}.paths.txt \\
        --min-af ${params.skani_min_af} \\
        -t ${task.cpus} \\
        -o ${comp_id}.skani.tsv
    """

    stub:
    """
    touch ${comp_id}.skani.tsv
    """
}
