// Compares MMSEQS_LINCLUST_BFD's and DIAMOND_CLUSTER's dereplication of the
// SAME BFD protein set (bin/compare_clusterings.py) -- Adjusted Rand Index
// plus basic per-method cluster stats. See that script's header for why ARI
// (not O(n^2) pairwise enumeration) is the right metric here.

process COMPARE_CLUSTERINGS {
    label 'compare_clusterings'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(mmseqs_clusters)
        path(diamond_clusters)

    output:
        path("bfd_clustering_comparison.tsv"), emit: comparison

    script:
    """
    ${projectDir}/bin/compare_clusterings.py \\
        --a mmseqs_linclust=${mmseqs_clusters} \\
        --b diamond_cluster=${diamond_clusters} \\
        --out bfd_clustering_comparison.tsv
    """

    stub:
    """
    printf "metric\\tmmseqs_linclust\\tdiamond_cluster\\n" > bfd_clustering_comparison.tsv
    """
}
