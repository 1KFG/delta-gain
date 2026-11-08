#!/usr/bin/env bash
# F2 within-family clustering, one array task per family (keys from results/family_regions/regions_summary.tsv)
# Submit: N=$(($(wc -l < $FUNGI/results/family_regions/regions_summary.tsv)-1)); sbatch --array=1-$N launch/family_cluster.sbatch
# Aggregate afterwards: python bin/family_cluster.py --aggregate --outdir $OUT --regions-dir . --tables-dir . --classes .
# Hardware: stajichlab has abu_dhabi nodes (c02,c21) with no AVX2 -> MMseqs2 dies; constraint avoids them.
#SBATCH -J family_cluster
#SBATCH -p stajichlab
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o logs/family_cluster.%A_%a.out
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
CLASSES=${CLASSES:-$FUNGI/results/delta_funnel/genome_classes.tsv}
REGIONS=${REGIONS:-$FUNGI/results/family_regions}
OUT=${OUT:-$FUNGI/results/family_clusters}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
SIF=${MMSEQS_SIF:-${NXF_APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}/depot.galaxyproject.org-singularity-mmseqs2-17.b804f--hd6d6fdc_1.img}
KEY=$(tail -n +2 "$REGIONS/regions_summary.tsv" | cut -f1 | sed -n "${SLURM_ARRAY_TASK_ID}p")
mkdir -p "$OUT"; TMP=${SCRATCH:-$OUT}/_tmp.$KEY.$$; mkdir -p "$TMP"; trap "rm -rf $TMP" EXIT
echo "[$(date)] cluster $KEY"
$PY -I $DG/bin/family_cluster.py --regions-dir "$REGIONS" --outdir "$OUT" --tables-dir "$TABLES" \
  --classes "$CLASSES" --families "$KEY" --threads $SLURM_CPUS_PER_TASK --tmp-dir "$TMP" \
  ${COV:+--cov $COV} ${COVMODE:+--cov-mode $COVMODE} \
  --mmseqs "apptainer exec --bind /bigdata --bind $TMP $SIF mmseqs"
echo "[$(date)] done"
