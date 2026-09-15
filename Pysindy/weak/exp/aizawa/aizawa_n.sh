#!/bin/bash
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH --job-name=psw_aizawa_n
# SINDy (Integration) -- PySINDy weak library + STLSQ at the benchmark threshold.
source activate "${WEAK_CONDA_ENV:-pysindy-prod}"
export OMP_NUM_THREADS=1
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS
export SYSTEM=aizawa
R CMD BATCH ../../weak_n.R aizawa_n_${START}_${END}_${SNR}.Rout
