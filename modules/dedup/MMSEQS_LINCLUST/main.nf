// Dereplicate the combined UniProtKB-fungi set (Swiss-Prot + TrEMBL,
// ~17.5M sequences) via MMseqs2 linclust -- the linear-time clustering
// algorithm, the only one of MMseqs2's clustering modes that scales to this
// size (the design doc's Step 2 already flagged cascaded `mmseqs cluster` as
// not scaling to hundreds-of-millions-of-sequence inputs; the same reasoning
// applies here, at smaller but still large scale).
//
// Output: a representative-sequence FASTA (one sequence per cluster) plus a
// cluster membership TSV (representative accession <TAB> member accession),
// so cluster-to-original-accession provenance is preserved for joining back
// to the InterPro/GO xrefs table from CONVERT_DAT_TO_FASTA.

process MMSEQS_LINCLUST {
    label 'mmseqs_linclust'

    input:
        path(combined_fasta)

    output:
        path("uniprot_fungi_nr.fasta"),      emit: representatives
        path("uniprot_fungi_clusters.tsv"),  emit: cluster_membership

    script:
    """
    mkdir -p mmseqs_tmp
    mmseqs createdb ${combined_fasta} DB

    mmseqs linclust DB DB_clu mmseqs_tmp \\
        --min-seq-id ${params.mmseqs_min_seq_id} \\
        -c ${params.mmseqs_coverage} \\
        --cov-mode ${params.mmseqs_cov_mode} \\
        --threads ${task.cpus}

    mmseqs createsubdb DB_clu DB DB_rep
    mmseqs convert2fasta DB_rep uniprot_fungi_nr.fasta

    mmseqs createtsv DB DB DB_clu uniprot_fungi_clusters.tsv
    """

    stub:
    """
    printf ">STUB_REP1 stub\\nMSTUBSEQ\\n" > uniprot_fungi_nr.fasta
    printf "STUB_REP1\\tSTUB_REP1\\n" > uniprot_fungi_clusters.tsv
    """
}
