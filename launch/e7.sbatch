#!/usr/bin/env bash
# E7: lineage-specific expansions (CAZy/MEROPS copy number) and signal-peptide / SSP-candidate novelty -> results/e7/
# stage 1 DuckDB tables, stage 2 cysteine counts from BFD FASTAs (multiprocessing), stage 3 within-order statistics.
#SBATCH -J e7
#SBATCH -p stajichlab
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 32
#SBATCH --mem=160G
#SBATCH -t 10:00:00
#SBATCH -o logs/e7.%j.out
set -o pipefail
source /etc/profile; set -e; ulimit -n 65535 2>/dev/null || ulimit -n 16384 2>/dev/null || true
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}; DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
OUT=${OUT:-$FUNGI/results/e7}; mkdir -p $OUT; TMP=$OUT/_tmp.$$; mkdir -p $TMP; trap "rm -rf $TMP" EXIT
PY=$FUNGI/.venv-pangenome/bin/python
CLASSES=$FUNGI/results/delta_funnel/genome_classes.tsv
if [[ ! -s $OUT/fam_counts.parquet ]]; then
  echo "[$(date)] stage 1"
  $PY -I $DG/bin/e7_stage1.py --tables-dir $TABLES --classes $CLASSES --outdir $OUT --threads 24 --memory-limit 120GB --tmp-dir $TMP
fi
if [[ ! -s $OUT/cys.parquet ]]; then
  echo "[$(date)] stage 2"
  $PY -I $DG/bin/e7_cys.py --cand $OUT/cand.parquet --meta $OUT/genome_meta.parquet --manifest $FUNGI/results/protein_novelty_full/bfd_pep_manifest.tsv --outdir $OUT --workers 30
fi
echo "[$(date)] stage 3"
$PY -I $DG/bin/e7_analyze.py --dir $OUT --survivors $FUNGI/results/survivor_annot/survivor_annot.tsv --perms ${PERMS:-3000} --boot ${BOOT:-1000}
echo "[$(date)] done"
