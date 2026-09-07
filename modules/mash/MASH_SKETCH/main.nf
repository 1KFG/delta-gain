// Per-genome mash sketch. Cheap and embarrassingly parallel — one task per
// genome, fans in for the single global MASH_PREFILTER_GLOBAL comparison.

process MASH_SKETCH {
    tag   { asmid }
    label 'mash_sketch'

    input:
        tuple val(asmid), path(genome_fa)

    output:
        tuple val(asmid), path("${asmid}.msh")

    script:
    """
    mash sketch -k ${params.mash_kmer} -s ${params.mash_sketch_size} \\
        -o ${asmid} ${genome_fa}
    """

    stub:
    """
    touch ${asmid}.msh
    """
}
