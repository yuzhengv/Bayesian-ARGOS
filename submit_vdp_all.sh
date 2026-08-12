#!/bin/bash
# Submit the Van der Pol (vdp) benchmark experiments to SLURM in one go.
#
# The argos-alasso vdp results already exist in
# Experiments/argos-alasso/vdp/results, so by default only the two missing
# methods (bayesian-alasso-ro and Pysindy) are submitted. Pass method names
# explicitly to override, e.g. to re-run argos-alasso as well.
#
# Usage:
#   bash submit_vdp_all.sh                  # submit bayesian-alasso-ro + pysindy
#   bash submit_vdp_all.sh pysindy          # submit one method only
#   bash submit_vdp_all.sh argos-alasso bayesian-alasso-ro pysindy
#   DRY_RUN=1 bash submit_vdp_all.sh        # print the submissions without running
#
# Each experiment folder's vdp_n_submit.sh and vdp_snr_submit.sh are executed
# from inside that folder, since they pass job scripts to sbatch by relative
# name.
set -euo pipefail

base_dir="$(cd "$(dirname "$0")" && pwd)"

declare -A method_dirs=(
    [argos-alasso]="${base_dir}/Experiments/argos-alasso/vdp"
    [bayesian-alasso-ro]="${base_dir}/Experiments/bayesian-alasso-ro/vdp"
    [pysindy]="${base_dir}/Pysindy/exp/vdp"
)

if [ "$#" -gt 0 ]; then
    methods=("$@")
else
    methods=(bayesian-alasso-ro pysindy)
fi

submitted=0
skipped=0
for method in "${methods[@]}"; do
    exp_dir="${method_dirs[$method]:-}"
    if [ -z "${exp_dir}" ]; then
        echo "WARNING: unknown method '${method}' (expected: ${!method_dirs[*]}), skipping"
        skipped=$((skipped + 1))
        continue
    fi
    if [ ! -d "${exp_dir}" ]; then
        echo "WARNING: ${exp_dir} does not exist, skipping"
        skipped=$((skipped + 1))
        continue
    fi
    for submit_script in vdp_n_submit.sh vdp_snr_submit.sh; do
        if [ ! -f "${exp_dir}/${submit_script}" ]; then
            echo "WARNING: ${exp_dir}/${submit_script} not found, skipping"
            skipped=$((skipped + 1))
            continue
        fi
        echo ">>> ${method}/${submit_script}"
        if [ "${DRY_RUN:-0}" != "0" ]; then
            echo "    (dry run) cd ${exp_dir} && bash ${submit_script}"
        else
            (cd "${exp_dir}" && bash "${submit_script}")
        fi
        submitted=$((submitted + 1))
    done
done

echo "----------------------------------------"
echo "Submit scripts executed: ${submitted}, skipped: ${skipped}"
