#!/usr/bin/env bash
# Regenerates the spurious-terms evidence for the reply to Referee 1 Point 3:
# the term-frequency and summary tables, the revised Fig. 6, the seven
# Supplementary term-selection figures and the Supplementary spurious-terms-by-
# system figure. The figure scripts copy their outputs to imgs/, si-plots/ and
# reply-to-reviewers/Figures/ of the manuscript repository when those folders
# exist (the submodule sits at <manuscript>/code-repos/Bayesian-ARGOS).
set -euo pipefail
cd "$(dirname "$0")"
Rscript spurious_terms_table.R > /dev/null
Rscript make_fig6.R > /dev/null
Rscript make_si_term_selection.R > /dev/null
Rscript make_si_spurious_by_system.R > /dev/null
echo "written: results/term_frequencies.csv results/spurious_summary.csv figures/shrinkage_spurious_terms.pdf figures/term_selection_<system>.pdf figures/spurious_terms_by_system.pdf"
