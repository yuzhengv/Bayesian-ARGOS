#!/usr/bin/env bash
# Weak-form check for Referee 1 Points 6-7 (SG derivative vs integral formulation).
#   Rscript weak_form_check.R [trials=20] [cases=dadras,aizawa,rossler]   ~35 min on 16 cores
#   Rscript make_figures_weak_form.R                                        seconds
set -euo pipefail
cd "$(dirname "$0")"
ARGOS_NCORES="${ARGOS_NCORES:-16}" Rscript weak_form_check.R 20
Rscript make_figures_weak_form.R
