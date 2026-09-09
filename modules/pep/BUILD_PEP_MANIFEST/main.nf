// Build the ASMID -> pep-FASTA-path manifest (bin/build_pep_manifest.py) --
// see that script's header for the exact makeSampleTag sanitization rule
// this reuses from Fungi_BFD/nextflow/modules/common/utils.nf, ported to
// Python rather than reimplemented ad hoc.
//
// pep_dir is passed as a val (a plain string path), NOT staged as a
// Nextflow `path` input -- it's a ~10k-file directory on the same shared
// GPFS filesystem every compute node already sees, same convention as
// params.nr_cluster_dmnd in DIAMOND_NOVELTY_SCREEN (an absolute path
// referenced directly, not copied/symlinked into the task work dir).

process BUILD_PEP_MANIFEST {
    label 'build_pep_manifest'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(samples_csv)
        val(pep_dir)
        val(n_test)

    output:
        path("bfd_pep_manifest.tsv"), emit: manifest

    script:
    """
    ${projectDir}/bin/build_pep_manifest.py \\
        --samples ${samples_csv} \\
        --pep-dir ${pep_dir} \\
        --n-test ${n_test} \\
        --out bfd_pep_manifest.tsv
    """

    stub:
    """
    printf "ASMID\\ttag\\tpep_path\\tstatus\\n" > bfd_pep_manifest.tsv
    printf "GCF_TEST001\\tstub_one\\t${projectDir}/tests/data/pep/stub_one.proteins.fa\\tmatched\\n" >> bfd_pep_manifest.tsv
    """
}
