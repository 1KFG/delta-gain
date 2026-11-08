#!/usr/bin/env bash
# E9 parameter sensitivity of the family-level analysis (docs/E9_RESULTS.md).
# Families: GH18, S08, pfam_PF00109 (backbone), pfam_PF10450 (enriched), pfam_PF01568 (random-control).
# Variants: env_cov80_m0 (default), env_cov60_m0, env_cov90_m0, env_cov80_m1 (target-only coverage), ali_cov80_m0 (Pfam regions by alignment span).
# Every variant clusters with the standard tiers (c90, c70, c50, c30) and runs the F3 baseline-vs-baseline null (200 draws).
source /etc/profile; module load slurm
FUNGI=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi; DG=/bigdata/stajichlab/shared/projects/BFD/DeltaGain
R=$FUNGI/results; E=$R/e9; cd $DG; mkdir -p $E
PART="-p short --constraint=ryzen|broadwell|cascade|sapphire"
KEYS="cazy_family_GH18 merops_family_S08 pfam_PF00109 pfam_PF10450 pfam_PF01568"

# --- regions: envelope (existing wide/ctrl extraction) as symlinks ------------------------------------------------------
mk_env() {
  d=$E/regions_env; mkdir -p $d
  head -n 1 $R/family_regions_wide/regions_summary.tsv > $d/regions_summary.tsv
  for k in $KEYS; do
    src=$R/family_regions_wide; [ -e $src/$k.faa ] || src=$R/family_regions_ctrl; [ -e $src/$k.faa ] || src=$R/family_regions_bgc
    ln -sf $src/$k.faa $d/$k.faa; ln -sf $src/$k.members.tsv $d/$k.members.tsv
    awk -F'\t' -v k=$k '$1==k' $src/regions_summary.tsv >> $d/regions_summary.tsv
  done
}
mk_env
echo "regions_env:"; cut -f1-3 $E/regions_env/regions_summary.tsv

# --- regions: alignment span for the three Pfam families --------------------------------------------------------------
if [ ! -s $E/regions_ali/regions_summary.tsv ]; then
  J=$(sbatch --parsable --wait $PART -c 16 --mem=100G -t 2:00:00 \
      --export=ALL,OUT=$E/regions_ali,AUTODUF=0,PFAMREGION=ali,FAMILIES="pfam:PF00109 pfam:PF10450 pfam:PF01568" launch/family_extract.sbatch)
  echo "extract ali job $J"
fi
cut -f1-3 $E/regions_ali/regions_summary.tsv

run_variant() {  # name regions cov covmode
  name=$1; reg=$2; cov=$3; cm=$4
  N=$(($(wc -l < $reg/regions_summary.tsv)-1))
  CL=$E/clusters_$name; NU=$E/null_$name
  sbatch --parsable --wait $PART --array=1-$N -t 2:00:00 \
     --export=ALL,REGIONS=$reg,OUT=$CL,COV=$cov,COVMODE=$cm launch/family_cluster.sbatch > /dev/null
  sbatch --parsable --wait $PART --array=1-$N -t 2:00:00 \
     --export=ALL,REGIONS=$reg,CLUSTERS=$CL,OUT=$NU,DRAWS=200 launch/family_null.sbatch > /dev/null
  echo "variant $name done $(date)"
}
run_variant env_cov80_m0 $E/regions_env 0.80 0 &
run_variant env_cov60_m0 $E/regions_env 0.60 0 &
run_variant env_cov90_m0 $E/regions_env 0.90 0 &
run_variant env_cov80_m1 $E/regions_env 0.80 1 &
run_variant ali_cov80_m0 $E/regions_ali 0.80 0 &
wait

# --- aggregate nulls ------------------------------------------------------------------------------------------------------
PY=$FUNGI/.venv-pangenome/bin/python; T=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables; C=$R/delta_funnel/genome_classes.tsv
for v in env_cov80_m0 env_cov60_m0 env_cov90_m0 env_cov80_m1 ali_cov80_m0; do
  reg=$E/regions_env; [ $v = ali_cov80_m0 ] && reg=$E/regions_ali
  sbatch --parsable --wait $PART -c 4 --mem=24G -t 0:30:00 -o logs/e9_agg_$v.%j.out \
    --wrap "$PY -I bin/family_null.py --aggregate --outdir $E/null_$v --members-dir $reg --clusters-dir $E/clusters_$v --tables-dir $T --classes $C; $PY -I bin/family_cluster.py --aggregate --outdir $E/clusters_$v --regions-dir $reg --tables-dir $T --classes $C"
done
echo ALLDONE
