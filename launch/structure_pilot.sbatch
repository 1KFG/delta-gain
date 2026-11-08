#!/usr/bin/env bash
#SBATCH -J esmfold_pilot
# partition given on the command line (exfab and short_gpu both tried; first to start wins)
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 2:00:00
#SBATCH -o logs/structure_pilot.%j.out
# ESMFold pilot on 400 no-hit profile-cluster representatives (stage inputs with bin/structure_pilot_select.py).
set -o pipefail
source /etc/profile; set -e
FUNGI=${FUNGI:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi}; DG=${DG:-/bigdata/stajichlab/shared/projects/BFD/DeltaGain}
export HF_HOME=$FUNGI/.hf_cache HF_HUB_OFFLINE=1
nvidia-smi -L || true
$FUNGI/.venv-esm/bin/python $DG/bin/structure_pilot_esmfold.py $FUNGI/results/structure_pilot/${FASTA:-pilot.faa} $FUNGI/results/structure_pilot/esmfold_${OUTSUB:-run}
