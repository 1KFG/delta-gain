// Resample the incidence matrix into accumulation curves, pan-proteome
// power-law fit, and per-clade Shapley contribution (design doc:
// ACCUMULATION_CURVE). params.n_permutations defaults to 1000, provisional
// pending the permutation-count convergence check (validation #1) -- see
// design doc's "Choosing production P".
//
// Python environment provisioning: same mechanism as
// BUILD_CLUSTER_INCIDENCE_MATRIX (see that module's header comment for the
// full rationale) -- params.python_container (a generic Python Apptainer/
// Singularity image) + pip_cache_dir (SETUP_PYTHON_PIP_CACHE's output,
// this run's own workDir/pip_cache), prepended to PYTHONPATH at runtime.
// params.python_container is empty by default; real value belongs in the
// calling profile config.

process ACCUMULATION_CURVE {
    label 'accumulation_curve'
    container params.python_container
    publishDir { "${params.outdir}/accumulation_${version_tag}" }, mode: 'copy'

    input:
        path(incidence_matrix)
        path(external_status)
        path(genome_metadata)
        val(version_tag)
        val(pip_cache_dir)

    output:
        path("accumulation_curve_summary.tsv"), emit: summary
        path("pangenome_powerlaw_fit.tsv"), emit: fit
        path("clade_contribution.tsv"), emit: clade_contribution
        path("permutation_marginals.parquet"), emit: raw_marginals

    script:
    """
    export PYTHONPATH="${pip_cache_dir}:\${PYTHONPATH:-}"
    ${projectDir}/bin/accumulation_curve.py \\
        --incidence-matrix ${incidence_matrix} \\
        --external-status ${external_status} \\
        --genome-metadata ${genome_metadata} \\
        --n-permutations ${params.n_permutations} \\
        --seed 0 \\
        --version-tag ${version_tag} \\
        --outdir .
    """

    stub:
    """
    touch accumulation_curve_summary.tsv pangenome_powerlaw_fit.tsv \\
          clade_contribution.tsv permutation_marginals.parquet
    """
}
