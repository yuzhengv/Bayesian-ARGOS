#!/usr/bin/env bash
# Regenerates the R-pipeline evidence for the reply to Referee 1 Point 6
# (and Referee 2 minor point 2).  Run from this folder.
#   Rscript exp1_collinearity.R  [trials_snr=20] [trials_n=10]
#   Rscript exp2_heteroscedasticity.R [trials=20]
#   Rscript aizawa_screening_vs_hmc.R      (stored benchmark results)
# then figures (make_figures.R, house style) and tables (the shared Python script, --tables-only)
set -euo pipefail
cd "$(dirname "$0")"
Rscript exp1_collinearity.R 20 10
Rscript exp2_heteroscedasticity.R 20
Rscript aizawa_screening_vs_hmc.R
Rscript aizawa_success_vs_snr_all_methods.R
Rscript exp6_lambda_rule.R 20
Rscript dadras_psis_loo.R                   # Referee 1 Point 7 (20 trials x 5 n)
Rscript make_figures.R                      # figures in the manuscript house style
Rscript make_figures_dadras_loo.R
PYARGOS_FIG="${PYARGOS_FIG:-../../../pyargos_with_sindy_shred/pyargos/results-diagnostics/rebuttal-diagnostics/make_figures.py}"
PY="${PY:-python}"
if [ -f "$PYARGOS_FIG" ]; then "$PY" "$PYARGOS_FIG" --res results --tables-only; fi   # tables/summary_tables.md
