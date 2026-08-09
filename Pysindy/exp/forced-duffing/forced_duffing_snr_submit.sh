#!/bin/bash
n_obs=5000
num_init=100 # 100
snr_start=1 # 2 Control the number of observations
snr_end=61 # 5
shell_by=60
snr_seq_end=$((snr_end-shell_by))
snr_loop_seq=$(seq $snr_start $shell_by $snr_seq_end)
snr_by=1
dt=0.01 # 0.01
seed=100 # 100

sg_poly_order=4 # 4
library_degree=5 # 5
library_type='poly' # 'poly','four','poly_four'

ncpus=12
memory=36

for loop_var in $snr_loop_seq; do
    sbatch --mem=${memory}G --export=ALL,N_OBS=${n_obs},NUM_INIT=${num_init},START=${loop_var},END=$(echo "$loop_var+$shell_by" | bc),BY_SNR=${snr_by},DT=${dt},SEED=${seed},POLY_ORDER=${sg_poly_order},LIBRARY_DEGREE=${library_degree},LIBRARY_TYPE=${library_type},CPU_NUM=${ncpus} forced_duffing_snr.sh
done
