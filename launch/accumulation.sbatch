#!/usr/bin/env bash
# E4 species accumulation: stage A builds (genome, tier-50 cluster) pairs, stage B runs the permutation curves.
#SBATCH -J accumulation
#SBATCH -p highclock
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o logs/accumulation.%j.out
set -o pipefail
source /etc/profile; set -e
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}; DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
OUT=${OUT:-$FUNGI/results/accumulation}; mkdir -p $OUT; TMP=$OUT/_tmp.$$; mkdir -p $TMP; trap "rm -rf $TMP" EXIT
PY=$FUNGI/.venv-pangenome/bin/python
if [[ ! -s $OUT/pairs.parquet ]]; then
  echo "[$(date)] stage pairs"
  $PY -I $DG/bin/accumulation.py --stage pairs --funnel-dir $FUNGI/results/delta_funnel --outdir $OUT --threads $SLURM_CPUS_PER_TASK \
     --memory-limit 150GB --tmp-dir $TMP
fi
echo "[$(date)] stage curves"
$PY -I $DG/bin/accumulation.py --stage curves --tables-dir $TABLES --classes $FUNGI/results/delta_funnel/genome_classes.tsv \
   --outdir $OUT --threads 8 --perms ${PERMS:-20}
echo "[$(date)] done"
