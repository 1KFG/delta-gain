#!/usr/bin/env bash
# F6: Pfam/CAZy/MEROPS labels on funnel survivors -> results/survivor_annot/
#SBATCH -J survivor_annot
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 2:00:00
#SBATCH -o logs/survivor_annot.%j.out
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
OUT=${OUT:-$FUNGI/results/survivor_annot}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
mkdir -p "$OUT"; TMP=$OUT/_tmp.$$; mkdir -p $TMP; trap "rm -rf $TMP" EXIT
$PY -I $DG/bin/survivor_annotate.py --survivors $FUNGI/results/delta_funnel/summary/survivors_clusters.tsv \
  --tables-dir "$TABLES" --outdir "$OUT" --threads $SLURM_CPUS_PER_TASK --memory-limit $((SLURM_MEM_PER_NODE*3/4/1024))GB --tmp-dir $TMP
