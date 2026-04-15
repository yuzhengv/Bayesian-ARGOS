#!/bin/bash
snr=49 # 49
num_init=100 # 100 
shell_start=2 # 2 Control the number of observations
shell_end=5 # 5
shell_by=1
shell_seq_end=$((shell_end-shell_by))
shell_loop_seq=$(seq $shell_start $shell_by $shell_seq_end)
by_time=0.1 # 0.1
dt=0.01 # 0.001 typically for lorenz system
seed=100 # 100

sg_poly_order=4 # 4 
library_degree=5 # 5
library_type='poly' # 'poly','four','poly_four'

state_var_start=1
state_var_end=3
state_var_seq=$(seq $state_var_start 1 $state_var_end)
ci_level=0.9 # ~ default setting is 0.9 

ncpus=24 # 25
memory=96

# Run the system identification experiments for the three governing equations of
# the halvorsen system respectively with differnt number of observations.
for state_var in $state_var_seq; do
    for loop_var in $shell_loop_seq; do
        sbatch --mem=${memory}G --export=ALL,SNR=${snr},NUM_INIT=${num_init},START=${loop_var},END=$(echo "$loop_var+$shell_by" | bc),BY_TIME=${by_time},DT=${dt},SEED=${seed},POLY_ORDER=${sg_poly_order},LIBRARY_DEGREE=${library_degree},LIBRARY_TYPE=${library_type},STATE_VAR=${state_var},CI_LEVEL=${ci_level},CPU_NUM=${ncpus} halvorsen_n.sh
    done
done