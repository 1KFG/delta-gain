#!/usr/bin/env bash
#SBATCH -J hgt_proxy
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 2:00:00
#SBATCH -o logs/hgt_proxy.%j.out
set -eo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}; DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
OUT=$FUNGI/results/hgt_proxy; TMP=$OUT/_tmp.$$; mkdir -p $TMP; trap "rm -rf $TMP" EXIT
$FUNGI/.venv-pangenome/bin/python -I $DG/bin/hgt_proxy.py --survivors $FUNGI/results/survivor_annot/survivor_annot.tsv \
  --tables-dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables --outdir $OUT --threads $SLURM_CPUS_PER_TASK \
  --memory-limit 90GB --tmp-dir $TMP
