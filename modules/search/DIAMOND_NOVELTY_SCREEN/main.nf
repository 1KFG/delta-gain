// Stage 3 novelty screen: DIAMOND blastp against the already-staged, taxonomy-
// aware nr_cluster DIAMOND db, restricted to Fungi at search time via
// --taxonlist -- no separate "fungi-only" database needs to be built or hosted.
//
// DB: /srv/projects/db/ncbi/diamond/20260128/nr_cluster_seq.dmnd
//   Built with --taxonmap/--taxonnodes/--taxonnames (see make_nr_dmnd.sh in that
//   dir), so --taxonlist works directly. 470M sequences, already clustered
//   (less redundant than full nr.dmnd's 990M), so prefer nr_cluster over nr here.
// Taxid 4751 = Fungi (NCBI Taxonomy).
//
// Output columns (outfmt 6 + stitle for quick manual triage):
//   qseqid sseqid pident length mismatch gapopen qstart qend sstart send
//   evalue bitscore qcovhsp scovhsp stitle
//
// Binning into strong-hit / weak-hit / no-hit happens downstream (not in this
// process) -- see bin/bin_novelty_hits.py -- using E-value + coverage, per the
// design doc's framing fix (a raw %identity cutoff wrongly flags paralogs and
// fast-evolving orthologs as "novel").

process DIAMOND_NOVELTY_SCREEN {
    tag   { chunk_id }
    label 'diamond'

    input:
        tuple val(chunk_id), path(query_faa)

    output:
        tuple val(chunk_id), path("${chunk_id}.nr_cluster_fungi.tsv")

    script:
    """
    diamond blastp \\
        --db ${params.nr_cluster_dmnd} \\
        --taxonlist ${params.fungi_taxid} \\
        --query ${query_faa} \\
        --evalue ${params.diamond_evalue} \\
        --sensitivity ${params.diamond_sensitivity} \\
        --max-target-seqs ${params.diamond_max_targets} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle \\
        --threads ${task.cpus} \\
        --out ${chunk_id}.nr_cluster_fungi.tsv
    """

    stub:
    """
    touch ${chunk_id}.nr_cluster_fungi.tsv
    """
}
