#!/usr/bin/env bash
# run_protein_novelty.sh — sbatch launcher for DeltaGain Stage 3
# (protein_novelty.nf): stage the comprehensive UniProtKB-fungi reference set
# (FTP bulk download + dat->fasta/xrefs conversion + MMseqs2 linclust
# dereplication) and run the small REST smoke test alongside it.
#
# Submit with:  sbatch run_protein_novelty.sh   [extra nextflow args...]
#
# This is a genuinely long-running pipeline: uniprot_trembl_fungi.dat.gz is
# ~13 GB (resumable download), and MMSEQS_LINCLUST on the combined ~17.5M-
# sequence set is configured for up to 48h on highmem. The head process must
# survive the whole thing (it stays alive to submit/monitor every child job),
# hence the long time limit below despite the head itself being light on
# cpu/mem.

#SBATCH -p epyc
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00
#SBATCH --job-name nf-DeltaGain-protein_novelty
#SBATCH -o logs/slurm/nf_%j.out
#SBATCH -e logs/slurm/nf_%j.err

set -euo pipefail

mkdir -p logs/slurm logs/nextflow lib/uniprot_fungi

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow
module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true

export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}"

PIPELINE="${PIPELINE:-$PWD}"

nextflow run "${PIPELINE}" -main-script protein_novelty.nf \
    -profile protein_novelty \
    -resume \
    "$@"

# Cheap validation before a real submit (run on a login node, seconds):
#   nextflow run . -main-script protein_novelty.nf -profile test -stub-run
