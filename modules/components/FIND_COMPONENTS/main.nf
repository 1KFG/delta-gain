process FIND_COMPONENTS {
    label 'components'
    publishDir "${params.outdir}/genome_classification", mode: 'copy'

    input:
        path(mash_dist)
        path(samples_csv)

    output:
        path("components.tsv"),                emit: components
        path("comp_*.genomes.tsv"), optional: true, emit: component_genome_lists

    script:
    """
    ${projectDir}/bin/find_components.py \\
        --mash-dist ${mash_dist} \\
        --samples ${samples_csv} \\
        --genome-dir ${params.genome_dir} \\
        --prefilter-ani ${params.prefilter_ani} \\
        --outdir .
    """

    stub:
    """
    printf "asmid\\tcomponent_id\\tcomponent_size\\texisting_species\\n" > components.tsv
    touch comp_00001.genomes.tsv
    """
}
