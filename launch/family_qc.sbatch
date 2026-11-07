#!/usr/bin/env bash
# F4 artefact/quality diagnostics for family novelty -> results/family_qc/{family_qc,qc_by_quality}.tsv
# Groups file (family_key<TAB>group) is built by launch/make_qc_groups.py
#SBATCH -J family_qc
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -o logs/family_qc.%j.out
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
R=$FUNGI/results
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
mkdir -p $R/family_qc
$PY -I $DG/launch/make_qc_groups.py $R > $R/family_qc/groups.tsv
$PY -I $DG/bin/family_qc.py --members-dir $R/family_regions_wide --clusters-dir $R/family_clusters_wide \
  --tables-dir "$TABLES" --groups $R/family_qc/groups.tsv --outdir $R/family_qc --threads $SLURM_CPUS_PER_TASK
