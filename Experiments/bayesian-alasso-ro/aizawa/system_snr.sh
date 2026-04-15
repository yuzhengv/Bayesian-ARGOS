#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 24
##SBATCH --mem=48G
#SBATCH -t 72:00:00
#SBATCH --job-name=bro_aizawa_snr

# Load Application
# module purge
# module load r/4.2.1
# module load python/3.9.9

# try to control automatic multithreading
ompthreads=1
export OMP_NUM_THREADS=$ompthreads
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

# Run application
R CMD BATCH system_snr.R system_${STATE_VAR}_${CI_LEVEL}_snr_${START}_${END}_${N_OBS}.Rout



