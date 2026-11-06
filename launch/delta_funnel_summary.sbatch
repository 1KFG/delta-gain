#!/usr/bin/env bash
# Re-run ONLY stage 4 (summary) of the funnel on an existing results/delta_funnel dir.
# Hardware: highclock = Xeon Gold (AVX2). Never abu_dhabi (DuckDB dies with Illegal instruction).
#SBATCH -J delta_funnel_sum
#SBATCH -p highclock
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o logs/delta_funnel_sum.%j.out
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
OUT=${OUT:-$FUNGI/results/delta_funnel}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
TMP=$OUT/_duckdb_tmp.$$; mkdir -p "$TMP"; trap "rm -rf $TMP" EXIT
cd "$OUT"
echo "[$(date)] 4. summarise"
$PY $DG/bin/delta_funnel_summary.py --prov funnel_provenance.tsv \
  --tier-clusters tier_95_clusters.tsv tier_70_clusters.tsv tier_50_clusters.tsv \
  --tier-labels 95 70 50 --outdir summary --threads $SLURM_CPUS_PER_TASK \
  --memory-limit $((SLURM_MEM_PER_NODE*3/4/1024))GB --tmp-dir "$TMP"
echo "[$(date)] done"
