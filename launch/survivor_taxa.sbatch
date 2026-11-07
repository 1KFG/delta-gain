#!/usr/bin/env bash
# E6: taxonomic breadth of survivor clusters -> results/survivor_taxa/
#SBATCH -J survivor_taxa
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 2:00:00
#SBATCH -o logs/survivor_taxa.%j.out
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
OUT=${OUT:-$FUNGI/results/survivor_taxa}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
mkdir -p "$OUT"; TMP=$OUT/_tmp.$$; mkdir -p $TMP; trap "rm -rf $TMP" EXIT
$PY -I $DG/bin/survivor_taxa.py --funnel-dir $FUNGI/results/delta_funnel --tables-dir "$TABLES" --outdir "$OUT" \
  --threads $SLURM_CPUS_PER_TASK --memory-limit $((SLURM_MEM_PER_NODE*3/4/1024))GB --tmp-dir $TMP
