#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 12
##SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH --job-name=ps_aizawa_snr

# Load Application
# module purge
# module load r/4.2.1
# module load python/3.9.9
# Load the specific conda environment
source activate pysindy-prod

# try to control automatic multithreading
ompthreads=1
export OMP_NUM_THREADS=$ompthreads
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

# Run application
R CMD BATCH system_snr.R system_snr_${START}_${END}_${N_OBS}.Rout



