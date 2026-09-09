// Turn raw DIAMOND hits (both targets, all chunks) into strong-hit /
// weak-hit / no-hit per BFD query protein (bin/bin_novelty_hits.py) --
// HANDOFF.md item 4. See that script's header for why the E-value/coverage
// thresholds used here are placeholder defaults, not yet empirically tuned.

process BIN_NOVELTY_HITS {
    label 'bin_novelty_hits'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(query_fasta)
        path(nr_cluster_hits, stageAs: 'nr_cluster_hits/*')
        path(uniprot_hits, stageAs: 'uniprot_hits/*')

    output:
        path("bfd_novelty_bins.tsv"), emit: bins

    script:
    // One chunk per SPLIT_FASTA piece per target -- bin_novelty_hits.py
    // accumulates repeated --hits specs for the same target label across
    // calls (see its header), so each chunk's TSV is passed as its own
    // --hits nr_cluster=... / --hits uniprot_fungi=... rather than
    // concatenated first.
    def nr_args      = nr_cluster_hits.collect { "--hits nr_cluster=${it}" }.join(' ')
    def uniprot_args = uniprot_hits.collect { "--hits uniprot_fungi=${it}" }.join(' ')
    """
    ${projectDir}/bin/bin_novelty_hits.py \\
        --queries ${query_fasta} \\
        ${nr_args} \\
        ${uniprot_args} \\
        --evalue-strong ${params.diamond_evalue} \\
        --out bfd_novelty_bins.tsv
    """

    stub:
    """
    printf "query_id\\toverall_status\\n" > bfd_novelty_bins.tsv
    """
}
