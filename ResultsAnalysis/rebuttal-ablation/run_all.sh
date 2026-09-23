#!/usr/bin/env bash
# Regenerates the numeric ablation summary for the reply to Referee 1 Point 1.
set -euo pipefail
cd "$(dirname "$0")"
Rscript ablation_success_table.R > /dev/null
echo "written: results/ablation_success_rates.csv tables/ablation_summary.md"
