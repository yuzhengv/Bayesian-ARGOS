#!/bin/bash
# Submit all Bayesian-ARGOS ablation experiments to SLURM in one go.
#
# Usage:
#   bash submit_all_ablation.sh                 # submit thomas + rossler (default)
#   bash submit_all_ablation.sh lorenz          # submit lorenz only
#   bash submit_all_ablation.sh thomas rossler lorenz
#   DRY_RUN=1 bash submit_all_ablation.sh      # print the submissions without running
#
# Each system folder's *_n_submit.sh and *_snr_submit.sh are executed from
# inside that folder, since they pass job scripts to sbatch by relative name.
set -euo pipefail

base_dir="$(cd "$(dirname "$0")" && pwd)"

variants=(
    bayesian-alasso-or
    bayesian-alasso-oo
    bayesian-alasso-rr
    bayesian-alasso-single-ols
    bayesian-alasso-single-ridge
)

if [ "$#" -gt 0 ]; then
    systems=("$@")
else
    systems=(thomas rossler)
fi

submitted=0
skipped=0
for variant in "${variants[@]}"; do
    for system in "${systems[@]}"; do
        exp_dir="${base_dir}/${variant}/${system}"
        if [ ! -d "${exp_dir}" ]; then
            echo "WARNING: ${exp_dir} does not exist, skipping"
            skipped=$((skipped + 1))
            continue
        fi
        for submit_script in "${system}_n_submit.sh" "${system}_snr_submit.sh"; do
            if [ ! -f "${exp_dir}/${submit_script}" ]; then
                echo "WARNING: ${exp_dir}/${submit_script} not found, skipping"
                skipped=$((skipped + 1))
                continue
            fi
            echo ">>> ${variant}/${system}/${submit_script}"
            if [ "${DRY_RUN:-0}" != "0" ]; then
                echo "    (dry run) cd ${exp_dir} && bash ${submit_script}"
            else
                (cd "${exp_dir}" && bash "${submit_script}")
            fi
            submitted=$((submitted + 1))
        done
    done
done

echo "----------------------------------------"
echo "Submit scripts executed: ${submitted}, skipped: ${skipped}"
