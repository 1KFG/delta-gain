// One-time, per-run pip install of the packages BUILD_CLUSTER_INCIDENCE_MATRIX
// and ACCUMULATION_CURVE need (pyarrow/duckdb/numpy) into workDir/pip_cache --
// scoped to THIS pipeline run's own work/ directory (shared storage, not
// $SCRATCH, which is node-local and would need reinstalling per task), not a
// long-lived cache shared across separate runs (2026-09-09 decision) --
// avoids managing package-version drift across pipeline versions, and
// cleans up naturally whenever this run's own work/ directory is cleaned up.
//
// params.python_container supplies a generic Python Apptainer/Singularity
// image with no need for these packages baked in; both downstream processes
// prepend this process's output directory to PYTHONPATH at runtime. Runs
// once per pipeline execution -- called unconditionally before both
// consumers in protein_novelty.nf's workflow{}, so Nextflow's own dependency
// ordering (not a beforeScript on each consumer, which would race if many
// task instances of that label ran concurrently) guarantees it completes
// before either downstream process starts.

process SETUP_PYTHON_PIP_CACHE {
    label 'setup_python_pip_cache'
    container params.python_container

    output:
        val("${workDir}/pip_cache"), emit: pip_cache_dir

    script:
    """
    mkdir -p ${workDir}/pip_cache
    pip install --no-input --target=${workDir}/pip_cache numpy pyarrow duckdb
    """

    stub:
    """
    mkdir -p ${workDir}/pip_cache
    """
}
