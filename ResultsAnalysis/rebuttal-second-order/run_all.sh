#!/usr/bin/env bash
# Regenerates the numeric summary of the second-order benchmark systems (reply to Referee 2 Point 2).
set -euo pipefail
cd "$(dirname "$0")"
Rscript second_order_success_table.R > /dev/null
echo "written: results/second_order_success_rates.csv tables/second_order_summary.md"
