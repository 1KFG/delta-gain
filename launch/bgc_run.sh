#!/usr/bin/env bash
# BGC-backbone proxy: family pipeline F1-F3b on secondary-metabolism backbone Pfams.
source /etc/profile; module load slurm
FUNGI=/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi; DG=/bigdata/stajichlab/shared/projects/BFD/DeltaGain
R=$FUNGI/results; cd $DG
FAMS="pfam:PF00109 pfam:PF00501 pfam:PF00668 pfam:PF03936 pfam:PF19086 pfam:PF01397 pfam:PF00975 pfam:PF14765 pfam:PF08659 pfam:PF00698 pfam:PF01040"
SB="sbatch --parsable --wait"
$SB --export=ALL,OUT=$R/family_regions_bgc,AUTODUF=0,FAMILIES="$FAMS" -p short -c 16 --mem=100G -t 2:00:00 launch/family_extract.sbatch
N=$(($(wc -l < $R/family_regions_bgc/regions_summary.tsv)-1)); echo "families $N"
sbatch --parsable --wait -p short --array=1-$N --export=ALL,REGIONS=$R/family_regions_bgc,OUT=$R/family_clusters_bgc -t 2:00:00 launch/family_cluster.sbatch
for k in null null_order; do
  sbatch --parsable --wait -p short --array=1-$N --export=ALL,REGIONS=$R/family_regions_bgc,CLUSTERS=$R/family_clusters_bgc,OUT=$R/family_${k}_bgc -t 2:00:00 launch/family_$k.sbatch &
done; wait
echo ALLDONE
