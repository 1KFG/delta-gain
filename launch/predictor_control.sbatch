#!/usr/bin/env bash
# Predictor-mismatch control for the funnel: see bin/predictor_control_prep.py
#SBATCH -J pred_control
#SBATCH -p short
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 32
#SBATCH --mem=160G
#SBATCH -t 2:00:00
#SBATCH -o logs/pred_control.%j.out
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
ASMDIR=${ASMDIR:-/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM}
OUT=${OUT:-$FUNGI/results/predictor_control}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
SIF=${MMSEQS_SIF:-/bigdata/stajichlab/shared/singularity_cache/depot.galaxyproject.org-singularity-mmseqs2-17.b804f--hd6d6fdc_1.img}
MMSEQS="apptainer exec --bind /bigdata --bind ${SCRATCH:-/tmp} $SIF mmseqs"
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMP=${SCRATCH:-/tmp}/predctl.$$; mkdir -p "$OUT" "$TMP"; trap 'rm -rf "$TMP"' EXIT
cd "$OUT"
$PY -I $DG/bin/predictor_control_prep.py --status $FUNGI/prior_annotation_status.tsv \
  --classes $FUNGI/results/delta_funnel/genome_classes.tsv \
  --manifest $FUNGI/results/protein_novelty_full/bfd_pep_manifest.tsv --asm-dir $ASMDIR --outdir "$OUT" \
  --test ${TEST:-300} --ref ${REF:-1000}
for RUN in A B; do
  mkdir -p $OUT/$RUN; cd $OUT/$RUN
  echo "[$(date)] run $RUN build"
  $PY -I $DG/bin/build_funnel_input.py --status $OUT/status_ref.tsv --asm-dir $ASMDIR --new-manifest $OUT/manifest_$RUN.tsv \
      --out-fasta funnel_input.faa --out-prov funnel_provenance.tsv --min-len 30
  cur=funnel_input.faa; k=0; files=(); labels=()
  for ID in 0.95 0.70 0.50; do
    lab=$(printf '%.0f' "$(echo "$ID*100" | bc -l)"); if [[ $k -eq 0 ]]; then cov=0.90; else cov=0.80; fi
    echo "[$(date)] $RUN linclust id=$ID cov=$cov"
    rm -rf "$TMP/t$k"; mkdir -p "$TMP/t$k"
    $MMSEQS createdb "$cur" "$TMP/t$k/DB" -v 1
    $MMSEQS linclust "$TMP/t$k/DB" "$TMP/t$k/clu" "$TMP/t$k/tmp" --min-seq-id "$ID" -c "$cov" --cov-mode 0 --threads "$THREADS" -v 1
    $MMSEQS createtsv "$TMP/t$k/DB" "$TMP/t$k/DB" "$TMP/t$k/clu" tier_${lab}_clusters.tsv --threads "$THREADS" -v 1
    $MMSEQS createsubdb "$TMP/t$k/clu" "$TMP/t$k/DB" "$TMP/t$k/rep" -v 1
    $MMSEQS convert2fasta "$TMP/t$k/rep" tier_${lab}_reps.faa -v 1
    files+=(tier_${lab}_clusters.tsv); labels+=("$lab"); cur=tier_${lab}_reps.faa; k=$((k+1))
  done
  $PY -I $DG/bin/delta_funnel_summary.py --prov funnel_provenance.tsv --tier-clusters "${files[@]}" --tier-labels "${labels[@]}" \
      --outdir summary --threads "$THREADS" --memory-limit "$((${SLURM_MEM_PER_NODE:-65536}*3/4/1024))GB" --tmp-dir "$TMP"
  rm -f tier_*_reps.faa funnel_input.faa
done
echo "[$(date)] done"
