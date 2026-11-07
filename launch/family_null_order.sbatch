#!/usr/bin/env bash
# F3b leave-one-order-out null, one array task per family (keys from results/family_regions/regions_summary.tsv)
# Submit: N=$(($(wc -l < $FUNGI/results/family_regions/regions_summary.tsv)-1)); sbatch --array=1-$N launch/family_null_order.sbatch
# Aggregate: python bin/family_null_order.py --aggregate --outdir $OUT --members-dir . --clusters-dir . --tables-dir . --classes .
# Hardware: stajichlab has abu_dhabi nodes (c02,c21) with no AVX2 -> numpy/DuckDB can die; constraint avoids them.
#SBATCH -J family_null_order
#SBATCH -p stajichlab
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 4:00:00
#SBATCH -o logs/family_null_order.%A_%a.out
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
CLASSES=${CLASSES:-$FUNGI/results/delta_funnel/genome_classes.tsv}
REGIONS=${REGIONS:-$FUNGI/results/family_regions}
CLUSTERS=${CLUSTERS:-$FUNGI/results/family_clusters}
OUT=${OUT:-$FUNGI/results/family_null_order}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
KEY=$(tail -n +2 "$REGIONS/regions_summary.tsv" | cut -f1 | sed -n "${SLURM_ARRAY_TASK_ID}p")
mkdir -p "$OUT"
echo "[$(date)] order-null $KEY draws=${DRAWS:-200}"
$PY -I $DG/bin/family_null_order.py --members-dir "$REGIONS" --clusters-dir "$CLUSTERS" --tables-dir "$TABLES" \
  --classes "$CLASSES" --outdir "$OUT" --families "$KEY" --draws ${DRAWS:-200} -K ${KGEN:-10} --min-base ${MINBASE:-20} --min-new ${MINNEW:-10}
echo "[$(date)] done"
