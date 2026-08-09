#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 25
##SBATCH --mem=100G
#SBATCH -t 24:00:00

# Load Application
# module purge
# module load r/4.2.1
# module load python/3.9.9

# try to control automatic multithreading
ompthreads=1
export OMP_NUM_THREADS=$ompthreads
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

# Run application
R CMD BATCH forced_duffing_n.R forced-duffing_${STATE_VAR}_${CI_LEVEL}_n_${START}_${END}_${SNR}.Rout
