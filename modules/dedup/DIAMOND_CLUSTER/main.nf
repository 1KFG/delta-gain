// Alternative BFD-protein dereplication path, run alongside MMSEQS_LINCLUST_BFD
// for direct comparison (2026-09-08 methodology question: does DIAMOND's
// cascaded clustering -- the same algorithm behind DIAMOND DeepClust,
// Buchfink et al. 2026 Nat Methods, doi:10.1038/s41592-026-03030-z --
// capture materially different structure than MMseqs2 linclust on the same
// protein set?). Not wired into the search/binning steps -- this branch
// exists purely to produce a second clustering for COMPARE_CLUSTERINGS.
//
// diamond/2.2.6 specifically (not the older 2.2.5 already cached): needed
// for --cluster-steps/--round-coverage cascaded-clustering support. New
// container built 2026-09-08 from quay.io/biocontainers/diamond:2.2.6--he361c42_0.
//
// --approx-id/--mutual-cover chosen to match MMSEQS_LINCLUST's existing
// 95%/90% thresholds (params.mmseqs_min_seq_id/params.mmseqs_coverage) for
// an apples-to-apples comparison, NOT independently tuned to DeepClust's
// own published recipe -- revisit if the comparison suggests it's worth
// tuning further.
//
// Only runs `diamond makedb`/`diamond cluster --out` here -- representative
// extraction is a SEPARATE process (EXTRACT_CLUSTER_REPS, below) in the
// python container. Real bug found 2026-09-08 running this for real on the
// actual BFD pilot data: the quay.io/biocontainers DIAMOND image has no
// python3 at all ("env: can't execute 'python3': No such file or
// directory"), so bin/extract_cluster_reps.py (needed anyway to work around
// a SEPARATE real DIAMOND bug -- `--reps`+`--out` together crashes, see
// that script's header) can't run inside this container. Splitting into two
// processes, each in the container that actually has what it needs, fixes
// both problems at once.

process DIAMOND_CLUSTER {
    tag       { prefix }
    label     'diamond_cluster'
    container params.diamond_container_2_2_6
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(prefix), path(combined_fasta)

    output:
        tuple val(prefix), path("${prefix}_diamond_clusters.tsv"), path(combined_fasta), emit: clusters_and_fasta

    script:
    """
    mkdir -p diamond_tmp
    diamond makedb --in ${combined_fasta} --db DB --threads ${task.cpus}

    diamond cluster --db DB.dmnd --out ${prefix}_diamond_clusters.tsv \\
        --approx-id ${params.diamond_cluster_approx_id} \\
        --mutual-cover ${params.diamond_cluster_mutual_cover} \\
        --threads ${task.cpus} \\
        --tmpdir diamond_tmp
    """

    stub:
    """
    printf "STUB_REP1\\tSTUB_REP1\\n" > ${prefix}_diamond_clusters.tsv
    """
}

process EXTRACT_CLUSTER_REPS {
    tag       { prefix }
    label     'extract_cluster_reps'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(prefix), path(clusters_tsv), path(combined_fasta)

    output:
        tuple val(prefix), path("${prefix}_diamond_nr.fasta"),      emit: representatives
        tuple val(prefix), path(clusters_tsv),                     emit: cluster_membership

    script:
    """
    ${projectDir}/bin/extract_cluster_reps.py \\
        --clusters ${clusters_tsv} \\
        --fasta ${combined_fasta} \\
        --out ${prefix}_diamond_nr.fasta
    """

    stub:
    """
    printf ">STUB_REP1 stub\\nMSTUBSEQ\\n" > ${prefix}_diamond_nr.fasta
    """
}
