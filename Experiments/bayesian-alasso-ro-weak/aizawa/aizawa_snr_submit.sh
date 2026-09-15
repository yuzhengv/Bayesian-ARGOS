#!/bin/bash
# SNR sweep, n = 5000, 100 initial conditions, SNR 1..61 dB + infinity (as Experiments/bayesian-alasso-ro/aizawa).
n_obs=5000
num_init=100
snr_start=1
snr_end=61
shell_by=60
snr_seq_end=$((snr_end-shell_by))
snr_loop_seq=$(seq $snr_start $shell_by $snr_seq_end)
snr_by=1
dt=0.01
seed=100
sg_poly_order=4
library_degree=5
library_type='poly'
state_var_seq=$(seq 1 1 3)
ci_level=0.9
ncpus=24
memory=48
for state_var in $state_var_seq; do
    for loop_var in $snr_loop_seq; do
        sbatch --mem=${memory}G --export=ALL,N_OBS=${n_obs},NUM_INIT=${num_init},START=${loop_var},END=$(echo "$loop_var+$shell_by" | bc),BY_SNR=${snr_by},DT=${dt},SEED=${seed},POLY_ORDER=${sg_poly_order},LIBRARY_DEGREE=${library_degree},LIBRARY_TYPE=${library_type},STATE_VAR=${state_var},CI_LEVEL=${ci_level},CPU_NUM=${ncpus} aizawa_snr.sh
    done
done
