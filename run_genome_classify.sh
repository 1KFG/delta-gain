#!/usr/bin/env bash
# run_genome_classify.sh — sbatch launcher for DeltaGain Stage 1
# (genome_classify.nf): label-free BFD genome classification via mash + skani.
#
# Submit with:  sbatch run_genome_classify.sh   [extra nextflow args...]
#
# The Nextflow head process needs to outlive your shell and must not run on
# the login node — it uses a light allocation and submits the real work as
# separate SLURM jobs via the slurm executor. Head partition time limit must
# cover the WHOLE pipeline (it stays alive to submit/monitor every child job).

#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00
#SBATCH --job-name nf-DeltaGain-genome_classify
#SBATCH -o logs/slurm/nf_%j.out
#SBATCH -e logs/slurm/nf_%j.err

set -euo pipefail

mkdir -p logs/slurm logs/nextflow

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow
module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true

# Container cache — same convention as Fungi_BFD/nextflow. Set one of these
# before submitting if not already exported in your environment.
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}"

# Local development: run from this checkout directly (no GitHub repo for
# DeltaGain yet). Once published, switch to the PIPELINE=<org>/<repo> pattern
# used by the sibling Fungi_BFD launchers.
PIPELINE="${PIPELINE:-$PWD}"

nextflow run "${PIPELINE}" \
    -profile genome_classify \
    -resume \
    "$@"

# Cheap validation before a real submit (run on a login node, seconds):
#   nextflow run . -profile test -stub-run --n_test 2
