process CLASSIFY_GENOMES {
    label 'classify'
    publishDir "${params.outdir}/genome_classification", mode: 'copy'

    input:
        path(components_tsv)
        path("skani_results/*")   // all comp_*.skani.tsv, gathered flat

    output:
        path("genome_classification.tsv")

    script:
    """
    ${projectDir}/bin/classify_genomes.py \\
        --components ${components_tsv} \\
        --skani-dir skani_results \\
        --ani-cluster-threshold ${params.ani_cluster_threshold} \\
        --ani-outlier-threshold ${params.ani_outlier_threshold} \\
        --out genome_classification.tsv
    """

    stub:
    """
    printf "asmid\\tcomponent_id\\tcomponent_size\\texisting_species\\tmajority_species_in_component\\tbest_within_component_ani\\tbest_ani_partner\\tcluster_class\\n" > genome_classification.tsv
    """
}
