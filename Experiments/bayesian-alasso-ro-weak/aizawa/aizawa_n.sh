#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH --job-name=brow_aizawa_n
# Bayesian-ARGOS (Integration) -- weak-form design from PySINDy, pipeline unchanged.
source activate "${WEAK_CONDA_ENV:-pysindy-prod}"
export OMP_NUM_THREADS=1
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS
export SYSTEM=aizawa
R CMD BATCH ../weak_n.R aizawa_${STATE_VAR}_${CI_LEVEL}_n_${START}_${END}_${SNR}.Rout
