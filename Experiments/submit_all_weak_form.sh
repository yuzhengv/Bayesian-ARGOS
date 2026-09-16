#!/bin/bash
# Submit the two weak-form (Integration) arms of the 2 x 2 comparison
#   Bayesian-ARGOS (SG)  x  Bayesian-ARGOS (Integration)
#   SINDy (SG)           x  SINDy (Integration)
# to SLURM. The SG arms are the stored benchmark results
# (Experiments/bayesian-alasso-ro, Pysindy/exp); only the Integration arms run.
#
# Usage:
#   bash submit_all_weak_form.sh                    # aizawa dadras rossler, both arms, both sweeps
#   bash submit_all_weak_form.sh dadras             # one system
#   ARMS="sindy" bash submit_all_weak_form.sh       # ARMS: "bayesian", "sindy" or "bayesian sindy"
#   SWEEPS="snr" bash submit_all_weak_form.sh       # SWEEPS: "n", "snr" or "n snr"
#   DRY_RUN=1 bash submit_all_weak_form.sh          # print, do not submit
#   WEAK_CONDA_ENV=<name> bash submit_all_weak_form.sh  # override pysindy-weak
set -euo pipefail
base_dir="$(cd "$(dirname "$0")" && pwd)"
root_dir="$(cd "${base_dir}/.." && pwd)"

# Hamilton environment verified by HAMILTON_ENV_CHECK.md. The child submit
# scripts use --export=ALL, so these settings also reach the batch jobs.
export WEAK_CONDA_ENV="${WEAK_CONDA_ENV:-pysindy-weak}"
# Site module / Conda initialization scripts may reference unset variables.
set +u
module load r/4.3.1
source activate "${WEAK_CONDA_ENV}"
set -u
export PYTHONNOUSERSITE=1 R_ENVIRON_USER=/dev/null
export RETICULATE_PYTHON="${CONDA_PREFIX}/bin/python"
export ARGOS_PYTHON="${RETICULATE_PYTHON}"
export R_LIBS="${CONDA_PREFIX}/r-library${R_LIBS:+:${R_LIBS}}"
export ARGOS_ROOT="${root_dir}" OMP_NUM_THREADS=1

systems=("$@"); [ "${#systems[@]}" -eq 0 ] && systems=(aizawa dadras rossler)
ARMS="${ARMS:-bayesian sindy}"; SWEEPS="${SWEEPS:-n snr}"
for arm in $ARMS; do
  for system in "${systems[@]}"; do
    case "$arm" in
      bayesian) exp_dir="${base_dir}/bayesian-alasso-ro-weak/${system}" ;;
      sindy)    exp_dir="${root_dir}/Pysindy/weak/exp/${system}" ;;
      *) echo "unknown arm $arm"; exit 1 ;;
    esac
    for sweep in $SWEEPS; do
      script="${exp_dir}/${system}_${sweep}_submit.sh"
      [ -f "$script" ] || { echo "WARNING: $script missing, skipping"; continue; }
      echo "== ${arm} ${system} ${sweep}: bash ${script}"
      [ "${DRY_RUN:-0}" = "1" ] || (cd "${exp_dir}" && bash "${script}")
    done
  done
done
