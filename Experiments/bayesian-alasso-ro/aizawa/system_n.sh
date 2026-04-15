#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH -n 24
##SBATCH --mem=32G
#SBATCH -t 48:00:00
#SBATCH --job-name=bro_aizawa_n

# Load Application
# module purge
# module load r/4.2.1
# module load python/3.9.9

# try to control automatic multithreading
ompthreads=1
export OMP_NUM_THREADS=$ompthreads
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

# Run application
R CMD BATCH system_n.R system_${STATE_VAR}_${CI_LEVEL}_n_${START}_${END}_${SNR}.Rout



