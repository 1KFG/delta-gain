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
//
// Python environment provisioning: this process needs pyarrow/duckdb/numpy,
// which a bare `python3` on this cluster does not have (resolves to
// miniconda py39 with none of them installed). params.python_container
// points at a generic Python Apptainer/Singularity image (real path
// supplied by the calling profile config, e.g. DeltaGain_Fungi's -c
// config, pulled from something like docker://python:3.11-slim into the
// shared singularity_cache); the container image itself does NOT need
// pyarrow/duckdb/numpy baked in -- pip_cache_dir (SETUP_PYTHON_PIP_CACHE's
// output, this run's own workDir/pip_cache) is prepended to PYTHONPATH so
// the container's python can import them without a custom-built image.
// params.python_container is empty by default (undefined here, like
// params.diamond_container) -- real value belongs in the calling profile
// config per the DeltaGain/DeltaGain_Fungi pipeline/data split.

process BUILD_CLUSTER_INCIDENCE_MATRIX {
    label 'build_cluster_incidence_matrix'
    container params.python_container
    publishDir { "${params.outdir}/accumulation_${version_tag}" }, mode: 'copy'

    input:
        path(clusters_tsv)
        path(provenance_tsv)
        path(novelty_bins_tsv)
        path(samples_csv)
        path(busco_parquet)
        path(asm_stats_parquet)
        path(genome_classification_tsv)
        val(version_tag)
        val(pip_cache_dir)

    output:
        path("incidence_matrix.parquet"), emit: incidence_matrix
        path("cluster_external_status.parquet"), emit: external_status
        path("genome_metadata.parquet"), emit: genome_metadata

    script:
    """
    export PYTHONPATH="${pip_cache_dir}:\${PYTHONPATH:-}"
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
