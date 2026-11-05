#!/usr/bin/env bash
# Delta-gain funnel: baseline = NCBI proteins of prior-annotated fungal genomes,
# query = BFD proteins of genomes with NO prior annotation.
# Tier 0 clusters baseline+query together at 95%/90% (linclust); each later tier
# clusters the previous tier's representatives at a lower identity.  A new
# protein is "novel" at a tier when its cluster has no baseline member.
# Everything is MMseqs2 linclust, so the whole funnel is CPU-hours, not DIAMOND-days.
#
# Required env (defaults shown):
#   STATUS      prior_annotation_status.tsv           (DeltaGain_Fungi/prior_annotation_status.tsv)
#   ASMDIR      NCBI_ASM folder                        (1KFG/.../source/NCBI_ASM)
#   PEP_MANIFEST bfd_pep_manifest.tsv                  (results/protein_novelty_full/bfd_pep_manifest.tsv)
#   OUT         output directory                       (results/delta_funnel)
#   TIERS       identity per tier, tier 0 first        ("0.95 0.70 0.50")
#   SIBLINGS    optional strain-sibling accession list
# Submit:  sbatch launch/delta_funnel.sbatch
#SBATCH -J delta_funnel
# Hardware: NEVER run on abu_dhabi nodes (all of `batch`, plus c02/c21 etc. in short/preempt/stajichlab).
# They lack AVX2 -> MMseqs2/DuckDB die with "Illegal instruction". epyc = Milan, intel = Broadwell; both modern.
# The constraint also excludes `highmem` Ivy Bridge nodes if a partition list is ever widened.
#SBATCH -p epyc,intel
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 64
#SBATCH --mem=400G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/delta_funnel.%j.out
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
STATUS=${STATUS:-$FUNGI/prior_annotation_status.tsv}
ASMDIR=${ASMDIR:-/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM}
PEP_MANIFEST=${PEP_MANIFEST:-$FUNGI/results/protein_novelty_full/bfd_pep_manifest.tsv}
OUT=${OUT:-$FUNGI/results/delta_funnel}
TIERS=${TIERS:-"0.95 0.70 0.50"}
COV0=${COV0:-0.90}      # tier-0 coverage (matches the pipeline's 95/90)
COV=${COV:-0.80}        # later tiers, bidirectional (cov-mode 0) like the metagenome-families paper
MINLEN=${MINLEN:-30}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
SIF=${MMSEQS_SIF:-${NXF_APPTAINER_CACHEDIR:-${APPTAINER_CACHEDIR:?set MMSEQS_SIF}}/depot.galaxyproject.org-singularity-mmseqs2-17.b804f--hd6d6fdc_1.img}
MMSEQS="apptainer exec --bind /bigdata --bind ${SCRATCH:-/tmp} $SIF mmseqs"
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMP=${SCRATCH:-/tmp}/delta_funnel.$$
mkdir -p "$OUT" "$TMP" "$FUNGI/logs"
cd "$OUT"
trap 'rm -rf "$TMP"' EXIT

echo "[$(date)] 1. classify genomes"
SIB=(); [[ -n "${SIBLINGS:-}" ]] && SIB=(--sibling-flags "$SIBLINGS")
$PY $DG/bin/classify_genomes_by_prior_annotation.py --pep-manifest "$PEP_MANIFEST" --status "$STATUS" \
    "${SIB[@]}" --out-classes genome_classes.tsv --out-new-manifest new_manifest.tsv

echo "[$(date)] 2. build baseline + query FASTA"
$PY $DG/bin/build_funnel_input.py --status "$STATUS" --asm-dir "$ASMDIR" --new-manifest new_manifest.tsv \
    --out-fasta funnel_input.faa --out-prov funnel_provenance.tsv --min-len "$MINLEN"

cur=funnel_input.faa; k=0; files=(); labels=()
for ID in $TIERS; do
  lab=$(printf '%.0f' "$(echo "$ID*100" | bc -l)")
  if [[ $k -eq 0 ]]; then cov=$COV0; else cov=$COV; fi
  echo "[$(date)] 3.$k linclust id=$ID cov=$cov on $cur"
  rm -rf "$TMP/t$k"; mkdir -p "$TMP/t$k"
  $MMSEQS createdb "$cur" "$TMP/t$k/DB" -v 1
  $MMSEQS linclust "$TMP/t$k/DB" "$TMP/t$k/clu" "$TMP/t$k/tmp" --min-seq-id "$ID" -c "$cov" --cov-mode 0 \
      --threads "$THREADS" -v 1
  $MMSEQS createtsv "$TMP/t$k/DB" "$TMP/t$k/DB" "$TMP/t$k/clu" tier_${lab}_clusters.tsv --threads "$THREADS" -v 1
  $MMSEQS createsubdb "$TMP/t$k/clu" "$TMP/t$k/DB" "$TMP/t$k/rep" -v 1
  $MMSEQS convert2fasta "$TMP/t$k/rep" tier_${lab}_reps.faa -v 1
  files+=(tier_${lab}_clusters.tsv); labels+=("$lab"); cur=tier_${lab}_reps.faa; k=$((k+1))
done

echo "[$(date)] 4. summarise"
$PY $DG/bin/delta_funnel_summary.py --prov funnel_provenance.tsv --tier-clusters "${files[@]}" \
    --tier-labels "${labels[@]}" --outdir summary --threads "$THREADS" --memory-limit "$((${SLURM_MEM_PER_NODE:-65536}*3/4/1024))GB" \
    --tmp-dir "$TMP"
echo "[$(date)] done. survivors for DIAMOND/profile tiers: $OUT/summary/survivors_clusters.tsv"
