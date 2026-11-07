#!/usr/bin/env bash
# E5 triage of survivors -> results/survivor_triage/
#SBATCH -J survivor_triage
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 1:00:00
#SBATCH -o logs/survivor_triage.%j.out
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
OUT=${OUT:-$FUNGI/results/survivor_triage}; mkdir -p $OUT
${PY:-$FUNGI/.venv-pangenome/bin/python} -I $DG/bin/survivor_triage.py --results $FUNGI/results --tables-dir $TABLES --outdir $OUT
