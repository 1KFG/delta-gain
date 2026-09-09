// Resample the incidence matrix into accumulation curves, pan-proteome
// power-law fit, and per-clade Shapley contribution (design doc:
// ACCUMULATION_CURVE). params.n_permutations defaults to 1000, provisional
// pending the permutation-count convergence check (validation #1) -- see
// design doc's "Choosing production P".
//
// Python environment provisioning (pyarrow/duckdb/numpy) must be supplied by
// the calling profile config for the 'accumulation_curve' label -- this repo
// declares no container/module/pixi directive here, per the DeltaGain/
// DeltaGain_Fungi pipeline/data split.

process ACCUMULATION_CURVE {
    label 'accumulation_curve'
    publishDir { "${params.outdir}/accumulation_${version_tag}" }, mode: 'copy'

    input:
        path(incidence_matrix)
        path(external_status)
        path(genome_metadata)
        val(version_tag)

    output:
        path("accumulation_curve_summary.tsv"), emit: summary
        path("pangenome_powerlaw_fit.tsv"), emit: fit
        path("clade_contribution.tsv"), emit: clade_contribution
        path("permutation_marginals.parquet"), emit: raw_marginals

    script:
    """
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
