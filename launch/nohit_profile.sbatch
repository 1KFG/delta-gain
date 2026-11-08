#!/usr/bin/env bash
# Profile search of 'no_hit' survivor clusters vs baseline-known tier-50 reps.
# 1 expand clusters to members  2 mmseqs align members->rep (-a)  3 result2profile  4 profile-vs-sequence search (-s 7.0)
# Hit (e<=1e-3, query coverage >=0.4) = distant homolog found only by profile sensitivity; no hit = orphan even by profile search.
# Test: sbatch --export=ALL,MAXC=2000,OUT=$FUNGI/results/nohit_profile_test ...
#SBATCH -J nohit_profile
#SBATCH -p stajichlab
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 32
#SBATCH --mem=250G
#SBATCH -t 1-12:00:00
#SBATCH -o logs/nohit_profile.%j.out
set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
FUN=${FUN:-$FUNGI/results/delta_funnel}
BEST=${BEST:-$FUNGI/results/survivor_search/survivor_search_best.tsv}
TARGET=${TARGET:-$FUNGI/results/survivor_search/target.faa}
OUT=${OUT:-$FUNGI/results/nohit_profile}
MAXC=${MAXC:-0}; SENS=${SENS:-7.0}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
SIF=${MMSEQS_SIF:-/bigdata/stajichlab/shared/singularity_cache/depot.galaxyproject.org-singularity-mmseqs2-17.b804f--hd6d6fdc_1.img}
MMSEQS="apptainer exec --bind /bigdata --bind ${SCRATCH:-/tmp} $SIF mmseqs"
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMP=${SCRATCH:-/tmp}/nohitprof.$$; mkdir -p "$OUT" "$TMP"; trap 'rm -rf "$TMP"' EXIT
cd "$OUT"
if [[ ! -s members.faa ]]; then
  echo "[$(date)] prep"
  $PY -I $DG/bin/nohit_profile_prep.py --funnel-dir $FUN --best $BEST --outdir . --max-clusters $MAXC \
      --memory-limit $((${SLURM_MEM_PER_NODE:-100000}*3/4/1024))GB --tmp-dir $TMP
fi
echo "[$(date)] createdb + cluster members (30 % id) + align -a + result2profile"
$MMSEQS createdb members.faa $TMP/mem -v 1
$MMSEQS cluster $TMP/mem $TMP/clu $TMP/ctmp --min-seq-id 0.3 -c 0.5 -s 5.0 --threads $THREADS -v 1
$MMSEQS align $TMP/mem $TMP/mem $TMP/clu $TMP/aln -a --threads $THREADS -v 1
$MMSEQS result2profile $TMP/mem $TMP/mem $TMP/aln $TMP/prof --threads $THREADS -v 1
$MMSEQS createtsv $TMP/mem $TMP/mem $TMP/clu profile_clusters.tsv --threads $THREADS -v 1
echo "[$(date)] profile search vs baseline-known reps"
$MMSEQS createdb $TARGET $TMP/t -v 1
$MMSEQS search $TMP/prof $TMP/t $TMP/res $TMP/tmp -s $SENS -e 1e-3 --max-seqs 20 -c 0.4 --cov-mode 2 --threads $THREADS \
   --split-memory-limit $((${SLURM_MEM_PER_NODE:-200000}*7/10/1024))G -v 1
$MMSEQS convertalis $TMP/prof $TMP/t $TMP/res profile_hits.tsv --format-output query,target,pident,evalue,bits,qcov,tcov --threads $THREADS -v 1
python3 - <<PY
import csv
best={}
for q,t,pid,ev,bits,qc,tc in csv.reader(open("profile_hits.tsv"),delimiter="\t"):
    ev=float(ev)
    if q not in best or ev<best[q]: best[q]=ev
n=nh=0; reps=set(); hit=set(); prot=0; prot_hit=0
for rep,mem in csv.reader(open("profile_clusters.tsv"),delimiter="\t"):
    reps.add(rep); prot+=1
    if rep in best:
        prot_hit+=1
for r in reps:
    if r in best: hit.add(r)
txt=f"profile_clusters\t{len(reps)}\nwith_hit\t{len(hit)}\nno_hit\t{len(reps)-len(hit)}\nproteins\t{prot}\nproteins_in_hit_clusters\t{prot_hit}\n"
open("profile_summary.txt","w").write(txt); print(txt)
PY
echo "[$(date)] done"
