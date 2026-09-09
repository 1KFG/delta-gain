// Build the genome x cluster incidence matrix + QC/taxonomy join (design
// doc: BUILD_CLUSTER_INCIDENCE_MATRIX). Joins BFD.duckdb's busco_genome/
// asm_stats Parquet tables by ASMID -- see that design doc section for why
// samples.csv alone is not enough (no BUSCO/N50 columns there).
//
// genome_classification_tsv is Stage 1's (genome_classify.nf) output, not
// produced by this workflow -- it is read from a completed run via
// params.genome_classification_tsv (see protein_novelty.nf's workflow{}
// block for how that path is turned into a channel), the same
// already-completed-output-by-path convention protein_novelty.nf already
// uses for params.uniprot_fungi_nr_fasta.

process BUILD_CLUSTER_INCIDENCE_MATRIX {
    label 'build_cluster_incidence_matrix'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(clusters_tsv)
        path(provenance_tsv)
        path(novelty_bins_tsv)
        path(samples_csv)
        path(busco_parquet)
        path(asm_stats_parquet)
        path(genome_classification_tsv)

    output:
        path("incidence_matrix.parquet"), emit: incidence_matrix
        path("cluster_external_status.parquet"), emit: external_status
        path("genome_metadata.parquet"), emit: genome_metadata

    script:
    """
    ${projectDir}/bin/build_cluster_incidence_matrix.py \\
        --clusters-tsv ${clusters_tsv} \\
        --provenance-tsv ${provenance_tsv} \\
        --novelty-bins-tsv ${novelty_bins_tsv} \\
        --samples-csv ${samples_csv} \\
        --busco-parquet ${busco_parquet} \\
        --asm-stats-parquet ${asm_stats_parquet} \\
        --genome-classification-tsv ${genome_classification_tsv} \\
        --min-clade-n ${params.min_clade_n} \\
        --outdir .
    """

    stub:
    """
    touch incidence_matrix.parquet cluster_external_status.parquet genome_metadata.parquet
    """
}
