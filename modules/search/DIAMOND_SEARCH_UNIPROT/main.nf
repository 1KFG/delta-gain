// Second Stage-3 novelty-screen target: DIAMOND blastp against the
// dereplicated UniProt-fungi representative set (DIAMOND_MAKEDB's output),
// complementing DIAMOND_NOVELTY_SCREEN's nr_cluster_seq.dmnd search --
// see HANDOFF.md item 3 ("DIAMOND search against nr_cluster_seq.dmnd
// and/or uniprot_fungi_nr.fasta").
//
// No --taxonlist here (unlike DIAMOND_NOVELTY_SCREEN): the db is already
// fungi-only by construction (built from the UniProt-fungi pull), so
// there's no broader-taxonomy nr db to restrict.
//
// Separate module rather than a parameterized single DIAMOND process:
// this db is a Nextflow-produced artifact staged normally as a `path`
// input (small enough relative to nr_cluster_seq.dmnd's ~195GB to be worth
// staging/symlinking), while nr_cluster_seq.dmnd is referenced directly by
// absolute path to avoid Nextflow ever trying to stage it -- different
// enough handling that folding both into one process would need more
// conditional branching than it's worth.

process DIAMOND_SEARCH_UNIPROT {
    tag       { chunk_id }
    label     'diamond'
    container params.diamond_container

    input:
        tuple val(chunk_id), path(query_faa)
        path(uniprot_dmnd)

    output:
        tuple val(chunk_id), path("${chunk_id}.uniprot_fungi.tsv")

    script:
    // See DIAMOND_NOVELTY_SCREEN's identical comment: DIAMOND's sensitivity
    // modes are standalone flags, not a `--sensitivity <value>` option.
    """
    diamond blastp \\
        --db ${uniprot_dmnd} \\
        --query ${query_faa} \\
        --evalue ${params.diamond_evalue} \\
        --${params.diamond_sensitivity} \\
        --max-target-seqs ${params.diamond_max_targets} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle \\
        --threads ${task.cpus} \\
        --out ${chunk_id}.uniprot_fungi.tsv
    """

    stub:
    """
    touch ${chunk_id}.uniprot_fungi.tsv
    """
}
