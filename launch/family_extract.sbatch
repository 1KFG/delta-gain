#!/usr/bin/env bash
# F1 region extraction for pilot families (+ top DUFs from the F0 census) -> results/family_regions/
# Hardware: stajichlab has abu_dhabi nodes (c02,c21) with no AVX2 -> DuckDB dies; constraint avoids them.
#SBATCH -J family_extract
#SBATCH -p stajichlab
#SBATCH --constraint=ryzen|broadwell|cascade|sapphire
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 8:00:00
#SBATCH -o logs/family_extract.%j.out
# Wide run: sbatch --export=ALL,OUT=$FUNGI/results/family_regions_wide,AUTOPFAM=200,AUTOLEVELS="cazy_family merops_family" -c 32 --mem=250G -t 24:00:00 launch/family_extract.sbatch
# Test: sbatch --export=ALL,MAXG=100,OUT=$FUNGI/results/family_regions_test -J fam_extract_test -c 8 --mem=32G -t 2:00:00 launch/family_extract.sbatch
set -euo pipefail
ulimit -n 16384 2>/dev/null || true   # one .faa + one .members.tsv open per family
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}
DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
TABLES=${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}
CLASSES=${CLASSES:-$FUNGI/results/delta_funnel/genome_classes.tsv}
MANIFEST=${MANIFEST:-$FUNGI/results/protein_novelty_full/bfd_pep_manifest.tsv}
CENSUS=${CENSUS:-$FUNGI/results/family_census}
OUT=${OUT:-$FUNGI/results/family_regions}
PY=${PY:-$FUNGI/.venv-pangenome/bin/python}
mkdir -p "$OUT"; TMP=$OUT/_duckdb_tmp.$$; mkdir -p "$TMP"; trap "rm -rf $TMP" EXIT
echo "[$(date)] family extract  max_genomes=${MAXG:-0} auto_duf=${AUTODUF:-10}"
$PY -I $DG/bin/family_extract.py --tables-dir "$TABLES" --classes "$CLASSES" --manifest "$MANIFEST" \
  --outdir "$OUT" --census-dir "$CENSUS" --auto-duf ${AUTODUF:-10} --auto-pfam ${AUTOPFAM:-0} ${FAMILIES:+--families $FAMILIES} ${AUTOLEVELS:+--auto-levels $AUTOLEVELS} ${PFAMREGION:+--pfam-region $PFAMREGION} --threads $SLURM_CPUS_PER_TASK \
  --memory-limit $((SLURM_MEM_PER_NODE*3/4/1024))GB --tmp-dir "$TMP" --max-genomes ${MAXG:-0}
echo "[$(date)] done"
