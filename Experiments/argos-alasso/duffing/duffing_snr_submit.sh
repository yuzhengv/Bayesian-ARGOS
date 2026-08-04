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

state_var_start=1
state_var_end=2
state_var_seq=$(seq $state_var_start 1 $state_var_end)
alpha_level=0.05 # 0.05
num_samples=2000 # 2000
sr_method='alasso' # 'lasso', 'alasso'
weights_method='ridge' # null, 'ols', 'ridge'
ols_ps=true # true
parallel='multicore' # 'multicore'
mc_ncpus=80
ncpus=20
memory=160


# Run the system identification experiments for the two governing equations of
# the duffing system respectively with differnt levels of noise.
for state_var in $state_var_seq; do
    for loop_var in $snr_loop_seq; do
        sbatch --mem=${memory}G --export=ALL,N_OBS=${n_obs},NUM_INIT=${num_init},START=${loop_var},END=$(echo "$loop_var+$shell_by" | bc),BY_SNR=${snr_by},DT=${dt},SEED=${seed},POLY_ORDER=${sg_poly_order},LIBRARY_DEGREE=${library_degree},LIBRARY_TYPE=${library_type},STATE_VAR=${state_var},ALPHA_LEVEL=${alpha_level},NUM_SAMPLE=${num_samples},SR=${sr_method},SR_RW=${weights_method},OLS=${ols_ps},PAR_CON=${parallel},CPU_NUM=${ncpus},MC_CPU_NUM=${mc_ncpus} duffing_snr.sh
    done
done
