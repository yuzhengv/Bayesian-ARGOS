#!/usr/bin/env bash
# Weak-form check for Referee 1 Points 6-7 (SG derivative vs integral formulation).
#   Rscript weak_form_check.R [trials=20] [cases=dadras,aizawa,rossler]   ~35 min on 16 cores
#   Rscript make_figures_weak_form.R                                        seconds
set -euo pipefail
cd "$(dirname "$0")"
ARGOS_NCORES="${ARGOS_NCORES:-16}" Rscript weak_form_check.R 20
Rscript make_figures_weak_form.R

# per-equation exact recovery and spurious-intercept counts, Dadras n sweep (Supplementary Note 4, 2026-11-08)
Rscript dadras_intercepts_by_equation.R

# conditioning at the paper's convention and Aizawa failure classes (Supplementary Note 4, 2026-11-08)
Rscript diag_conditioning.R
Rscript aizawa_failure_patterns_by_arm.R
