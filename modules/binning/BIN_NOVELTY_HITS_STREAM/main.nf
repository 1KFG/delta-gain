// Streaming variant of BIN_NOVELTY_HITS (bin/bin_novelty_hits_stream.py).
// Takes each chunk's FASTA together with that chunk's two hit files, in chunk
// order, and builds the manifest the script reads. Memory is bounded by the
// largest chunk instead of all queries. Selected with --bin_streaming true.
//
// All four inputs must be lists in the SAME chunk order (the workflow sorts
// by chunk id before splitting them). The script's --queries check verifies
// that chunk order equals the order of the full query FASTA and fails loudly
// if it does not.

process BIN_NOVELTY_HITS_STREAM {
    label 'bin_novelty_hits'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(query_fasta)
        path(chunk_fastas, stageAs: 'chunks/*')
        path(nr_cluster_hits, stageAs: 'nr_cluster_hits/*')
        path(uniprot_hits, stageAs: 'uniprot_hits/*')

    output:
        path("bfd_novelty_bins.tsv"), emit: bins

    script:
    def fa = [chunk_fastas].flatten()
    def nr = [nr_cluster_hits].flatten()
    def up = [uniprot_hits].flatten()
    assert fa.size() == nr.size() && fa.size() == up.size() :
        "BIN_NOVELTY_HITS_STREAM: ${fa.size()} chunk FASTAs, ${nr.size()} nr hit files, ${up.size()} uniprot hit files"
    def rows = (0..<fa.size()).collect { i ->
        "printf '%s\\t%s\\t%s\\t%s\\n' ${fa[i].baseName} ${fa[i]} ${nr[i]} ${up[i]} >> manifest.tsv"
    }.join('\n')
    """
    printf 'chunk\\tqueries\\tnr_cluster\\tuniprot_fungi\\n' > manifest.tsv
    ${rows}
    ${projectDir}/bin/bin_novelty_hits_stream.py \\
        --manifest manifest.tsv \\
        --queries ${query_fasta} \\
        --evalue-strong ${params.diamond_evalue} \\
        --out bfd_novelty_bins.tsv
    """

    stub:
    """
    printf "query_id\\toverall_status\\n" > bfd_novelty_bins.tsv
    """
}
