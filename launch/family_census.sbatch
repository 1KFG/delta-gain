#!/usr/bin/env bash
# F0 family census (Pfam / CAZy / MEROPS) -> results/family_census/*.tsv
# Hardware: stajichlab has abu_dhabi nodes (c02,c21) with no AVX2 -> DuckDB dies; constraint avoids them.
#SBATCH -J family_census
#SBATCH -p stajichlab
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o logs/family_census.%j.out
# Test:  sbatch --export=ALL,MAXG=200,OUT=$FUNGI/results/family_census_test -c 8 --mem=64G -t 2:00:00 launch/family_census.sbatch
set -euo pipefail
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
CLASSES=${CLASSES:-$FUNGI/results/delta_funnel/genome_classes.tsv}
OUT=${OUT:-$FUNGI/results/family_census}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
mkdir -p "$OUT"; TMP=$OUT/_duckdb_tmp.$$; mkdir -p "$TMP"; trap "rm -rf $TMP" EXIT
echo "[$(date)] family census  tables=$TABLES classes=$CLASSES max_genomes=${MAXG:-0}"
$PY -I $DG/bin/family_census.py --tables-dir "$TABLES" --classes "$CLASSES" --outdir "$OUT" \
  --threads $SLURM_CPUS_PER_TASK --memory-limit $((SLURM_MEM_PER_NODE*3/4/1024))GB --tmp-dir "$TMP" \
  --max-genomes ${MAXG:-0}
echo "[$(date)] done"
