#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 80
##SBATCH --mem=160G
#SBATCH -t 72:00:00

# Load Application
# module purge
# module load r/4.2.1
# module load python/3.9.9

# try to control automatic multithreading
ompthreads=1
export OMP_NUM_THREADS=$ompthreads
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

# Run application
R CMD BATCH vdp_snr.R vdp_${STATE_VAR}_${SR}_snr_${START}_${END}_${N_OBS}.Rout



