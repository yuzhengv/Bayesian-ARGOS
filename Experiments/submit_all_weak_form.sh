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
set -euo pipefail
base_dir="$(cd "$(dirname "$0")" && pwd)"
root_dir="$(cd "${base_dir}/.." && pwd)"
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
