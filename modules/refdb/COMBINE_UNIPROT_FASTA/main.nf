// Combine Swiss-Prot + TrEMBL fungi (FASTA and xrefs) into one comprehensive
// UniProtKB-fungi comparison set, feeding MMSEQS_LINCLUST.

process COMBINE_UNIPROT_FASTA {
    label 'combine_uniprot'
    publishDir "${params.uniprot_staging_dir}", mode: 'copy'

    input:
        path(fastas)   // sprot_fungi.fasta, trembl_fungi.fasta
        path(xrefs)    // sprot_fungi.xrefs.tsv, trembl_fungi.xrefs.tsv

    output:
        path("uniprot_fungi_combined.fasta"),      emit: fasta
        path("uniprot_fungi_combined.xrefs.tsv"),  emit: xrefs

    script:
    """
    cat ${fastas} > uniprot_fungi_combined.fasta

    head -n1 ${xrefs[0]} > uniprot_fungi_combined.xrefs.tsv
    for f in ${xrefs}; do
        tail -n +2 "\$f" >> uniprot_fungi_combined.xrefs.tsv
    done
    """

    stub:
    """
    touch uniprot_fungi_combined.fasta uniprot_fungi_combined.xrefs.tsv
    """
}
