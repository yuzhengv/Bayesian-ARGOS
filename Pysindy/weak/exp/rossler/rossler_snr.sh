#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH --job-name=psw_rossler_snr
# SINDy (Integration) -- PySINDy weak library + STLSQ at the benchmark threshold.
source activate "${WEAK_CONDA_ENV:-pysindy-prod}"
export OMP_NUM_THREADS=1
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS
export SYSTEM=rossler
R CMD BATCH ../../weak_snr.R rossler_snr_${START}_${END}_${N_OBS}.Rout
