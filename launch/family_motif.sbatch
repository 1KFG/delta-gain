#!/usr/bin/env bash
#SBATCH -J family_motif
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 1:00:00
#SBATCH -o logs/family_motif.%j.out
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}; DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
R=$FUNGI/results
$FUNGI/.venv-pangenome/bin/python -I $DG/bin/family_motif.py --regions-dir $R/family_regions_wide --clusters-dir $R/family_clusters_wide \
  --tables-dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables --outdir $R/family_deepdive
