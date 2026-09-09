// Concatenate all `matched` BFD proteins (per BUILD_PEP_MANIFEST) into one
// combined FASTA for MMSEQS_LINCLUST, prefixing every ID with its ASMID
// (bin/combine_bfd_proteins.py) so per-protein genome provenance survives
// dereplication -- see that script's header comment.
//
// Reads pep_path values straight out of the manifest rather than taking
// each pep FASTA as a Nextflow `path` input: up to ~11k separate files today
// (more once the rest of the 22k genomes are annotated) is not a channel
// Nextflow should stage/symlink individually -- same absolute-path-direct
// convention as BUILD_PEP_MANIFEST's pep_dir.

process COMBINE_BFD_PROTEINS {
    label 'combine_bfd_proteins'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(manifest_tsv)

    output:
        path("bfd_proteins_combined.fasta"),   emit: fasta
        path("bfd_proteins_provenance.tsv"),   emit: provenance

    script:
    """
    ${projectDir}/bin/combine_bfd_proteins.py \\
        --manifest ${manifest_tsv} \\
        --out-fasta bfd_proteins_combined.fasta \\
        --out-provenance bfd_proteins_provenance.tsv
    """

    stub:
    """
    printf ">STUB_ASMID__stub1 stub\\nMSTUBSEQ\\n" > bfd_proteins_combined.fasta
    printf "new_id\\tasmid\\toriginal_id\\n" > bfd_proteins_provenance.tsv
    printf "STUB_ASMID__stub1\\tSTUB_ASMID\\tstub1\\n" >> bfd_proteins_provenance.tsv
    """
}
