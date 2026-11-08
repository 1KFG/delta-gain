#!/usr/bin/env bash
#SBATCH -J catalytic
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 2:00:00
#SBATCH -o logs/catalytic_align.%j.out
set -x
source /etc/profile; module load hmmer/3.4; set -e
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}; DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
R=$FUNGI/results
$FUNGI/.venv-pangenome/bin/python -I $DG/bin/catalytic_align.py --regions-dir $R/family_regions_wide --clusters-dir $R/family_clusters_wide \
  --hmm-dir $R/catalytic --outdir $R/catalytic
