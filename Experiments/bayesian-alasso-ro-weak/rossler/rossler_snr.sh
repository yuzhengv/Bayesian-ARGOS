#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH --job-name=brow_rossler_snr
# Bayesian-ARGOS (Integration) -- weak-form design from PySINDy, pipeline unchanged.
# Needs R with rstanarm/glmnet and a Python with pysindy >= 1.7.3 on the path
# (the pysindy-prod conda environment of the SINDy benchmarks, updated).
source activate "${WEAK_CONDA_ENV:-pysindy-prod}"
export OMP_NUM_THREADS=1
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS
export SYSTEM=rossler
R CMD BATCH ../weak_snr.R rossler_${STATE_VAR}_${CI_LEVEL}_snr_${START}_${END}_${N_OBS}.Rout
