// Build a DIAMOND db from the dereplicated UniProt-fungi representative
// FASTA (MMSEQS_LINCLUST's uniprot_fungi_nr.fasta) so BFD proteins can be
// searched against it directly -- unlike nr_cluster_seq.dmnd (already
// staged, taxonomy-aware, referenced by absolute path), this reference set
// is itself a pipeline product, so it needs its own db build step.

process DIAMOND_MAKEDB {
    label     'diamond_makedb'
    container params.diamond_container
    publishDir "${params.refdb_outdir}", mode: 'copy'

    input:
        path(fasta)

    output:
        path("${fasta.baseName}.dmnd"), emit: dmnd

    script:
    """
    diamond makedb --in ${fasta} --db ${fasta.baseName} --threads ${task.cpus}
    """

    stub:
    """
    touch ${fasta.baseName}.dmnd
    """
}
