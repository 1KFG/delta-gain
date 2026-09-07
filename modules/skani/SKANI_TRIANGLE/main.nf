// Exact-ANI refinement within one mash-prefiltered connected component.
// Components are typically small (one candidate species/strain cluster), so
// this is the expensive-but-accurate step, run only where mash already found
// enough similarity to be worth refining -- never all-vs-all across all 22k.

include { skaniCpusFor; skaniMemoryFor; skaniTimeFor } from '../../../lib/utils.nf'

process SKANI_TRIANGLE {
    tag    { comp_id }
    label  'skani'

    cpus   { skaniCpusFor(genome_list.countLines()) }
    memory { skaniMemoryFor(genome_list.countLines(), task.attempt) }
    time   { skaniTimeFor(genome_list.countLines(), task.attempt) }

    input:
        tuple val(comp_id), path(genome_list)

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
