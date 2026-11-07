#!/usr/bin/env bash
# Stage 5: sensitive search of survivor (novel-at-50%) cluster representatives against baseline-known tier-50 reps.
# Hit (e<=1e-5, query coverage >= 0.4) = distant homolog of a known cluster; no hit = novel even by sensitive search.
# Test: sbatch --export=ALL,MAXQ=20000,OUT=$FUNGI/results/survivor_search_test -c 32 --mem=200G -t 4:00:00 launch/survivor_search.sbatch
#SBATCH -J survivor_search
#SBATCH -p stajichlab
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 64
#SBATCH --mem=400G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/survivor_search.%j.out
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
FUN=${FUN:-$FUNGI/results/delta_funnel}
OUT=${OUT:-$FUNGI/results/survivor_search}
MING=${MING:-5}; MAXQ=${MAXQ:-0}; SENS=${SENS:-6.0}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
SIF=${MMSEQS_SIF:-/bigdata/stajichlab/shared/singularity_cache/depot.galaxyproject.org-singularity-mmseqs2-17.b804f--hd6d6fdc_1.img}
MMSEQS="apptainer exec --bind /bigdata --bind ${SCRATCH:-/tmp} $SIF mmseqs"
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMP=${SCRATCH:-/tmp}/survsearch.$$; mkdir -p "$OUT" "$TMP"; trap 'rm -rf "$TMP"' EXIT
cd "$OUT"
if [[ ! -s query.faa ]]; then
  echo "[$(date)] split"
  $PY -I $DG/bin/survivor_split.py --reps $FUN/tier_50_reps.faa --survivors $FUN/summary/survivors_clusters.tsv \
     --min-genomes $MING --max-query $MAXQ --out-query query.faa --out-target target.faa
fi
echo "[$(date)] createdb"
$MMSEQS createdb query.faa "$TMP/q" -v 1
$MMSEQS createdb target.faa "$TMP/t" -v 1
echo "[$(date)] search s=$SENS"
$MMSEQS search "$TMP/q" "$TMP/t" "$TMP/res" "$TMP/tmp" -s $SENS -e 1e-5 --max-seqs 20 -c 0.4 --cov-mode 2 \
   --threads $THREADS --split-memory-limit $((${SLURM_MEM_PER_NODE:-200000}*7/10/1024))G -v 1
$MMSEQS convertalis "$TMP/q" "$TMP/t" "$TMP/res" hits.tsv --format-output query,target,pident,evalue,bits,qcov,tcov --threads $THREADS -v 1
echo "[$(date)] summarise"
$PY -I $DG/bin/survivor_search_summary.py --hits hits.tsv --survivors $FUN/summary/survivors_clusters.tsv --query-fasta query.faa --outdir .
echo "[$(date)] done"
