#!/bin/bash
# n sweep, SNR 49 dB, 100 initial conditions, n = 10^2 .. 10^5 in steps of 10^0.1 (as Experiments/bayesian-alasso-ro/rossler).
snr=49
num_init=100
shell_start=2
shell_end=5
shell_by=1
shell_seq_end=$(awk "BEGIN {print $shell_end - $shell_by}")
shell_loop_seq=$(awk "BEGIN {for (i=$shell_start; i<=$shell_seq_end; i+=$shell_by) printf i\" \"}")
by_time=0.1
dt=0.01
seed=100
sg_poly_order=4
library_degree=5
library_type='poly'
state_var_seq=$(seq 1 1 3)
ci_level=0.9
ncpus=24
memory=96
for state_var in $state_var_seq; do
    for loop_var in $shell_loop_seq; do
        sbatch --mem=${memory}G --export=ALL,SNR=${snr},NUM_INIT=${num_init},START=${loop_var},END=$(echo "$loop_var+$shell_by" | bc),BY_TIME=${by_time},DT=${dt},SEED=${seed},POLY_ORDER=${sg_poly_order},LIBRARY_DEGREE=${library_degree},LIBRARY_TYPE=${library_type},STATE_VAR=${state_var},CI_LEVEL=${ci_level},CPU_NUM=${ncpus} rossler_n.sh
    done
done
