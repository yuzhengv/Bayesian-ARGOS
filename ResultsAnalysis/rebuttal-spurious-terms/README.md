# Rebuttal spurious terms (Referee 1, Point 3)

**Status (2026-10-14): complete.** Aggregate evidence on how many spurious
terms each method selects and how strongly, computed from the stored 100-trial
benchmark results of the seven systems, plus the redesigned term-selection
figures: the revised Fig. 6 (aggregate panels and one circular panel with the
three methods side by side), the seven Supplementary term-selection figures
(every system, every method, both sweeps) and a Supplementary figure of the
spurious terms per equation for every system. Used for the reply to Referee 1
Point 3; the figure files are copied to `imgs/`, `si-plots/` and
`reply-to-reviewers/Figures/` of the manuscript repository.

## Referee's question

> "Circular diagrams are hard to read: they either need an instruction on how
> to read them, or a different visualization. Fig 5 makes comparison
> unnecessarily hard - panel a compares Bayesian-ARGOS vs SINDy on Lorenz under
> varying n; panel b compares Bayesian-ARGOS to ARGOS on Thomas under varying
> SNR. This feels like a strong cherrypicking from all available evidence in
> parameter sweeps. It is less crucial to know *which* spurious terms appear in
> hard-identification cases, so maybe a different visualization could be
> designed to convey *how many and how strongly* do spurious terms appear."

## What is computed

For every system (Lorenz, Thomas, Rössler, Dadras, Aizawa, Sprott, Halvorsen),
method (Bayesian-ARGOS = `Experiments/bayesian-alasso-ro`, ARGOS =
`Experiments/argos-alasso`, SINDy = `Pysindy/exp`), sweep (n from 10^2 to 10^5
in steps of 0.1 in log10 n at SNR = 49 dB; SNR from 1 to 61 dB in 1 dB steps
plus the noiseless case at n = 5000) and grid point, from the final identified
models of the 100 trials:

- **Selection frequency** of a term in an equation: the fraction of the 100
  trials in which the term is in the final model of that equation
  (`results/term_frequencies.csv`; true terms never selected are listed with
  frequency 0). Term labels are the LaTeX labels the circular plots use, so
  that the R candidate library of ARGOS and Bayesian-ARGOS and the PySINDy
  library of SINDy meet in one naming; the intercept is `$C$`.
- **Spurious terms per equation**: per trial, the number of selected terms that
  are not true terms, summed over the three equations and divided by 3; then
  the mean over the 100 trials.
- **True terms missed per equation**: per trial, the number of true terms not
  selected, summed over the three equations and divided by 3; then the mean
  over trials.
- **Mean selection frequency of the spurious terms** ("how strongly"): the
  mean, over the spurious terms selected in at least one trial at that
  condition (the three equations pooled), of their selection frequency; 0 when
  no spurious term is selected in any trial. A value near 1 means the same
  spurious terms return in every trial (systematic); a value near 0.02 means
  many different terms appear once or twice (sporadic). The table also carries
  the largest selection frequency among the spurious terms and the number of
  distinct spurious terms, for reference.
- **Success**: the fraction of trials with all three equations recovered
  exactly (no spurious and no missed term). This reproduces the success-rate
  notebooks' criterion and is checked against the stored success-rate tables
  of `rebuttal-ablation/` (Ridge-OLS baseline, three systems) and
  `rebuttal-diagnostics/aizawa_success_vs_snr_all_methods.csv` (three
  methods): 465 grid points, maximum absolute difference 0
  (`results/success_rate_check.csv`).

Means over the seven systems per method and grid point, with the minimum and
maximum over the systems, are in `results/spurious_summary_aggregate.csv`
(block a of Fig. 6).

## Departures from the manuscript pipeline

None. The scripts only re-read the stored results through the circular-plot
notebooks' own loaders (`enhance_total_results_summary_df` of
`MethodsEvaluation/evaluation-functions-for-plots/` and
`enhance_total_results_summary_df_py` of `Pysindy/utils/utils-functions-for-plots/`,
unchanged, each sourced into its own environment because the four modules
define the same function names). The true terms are those of the success-rate
notebooks in both naming schemes (`common.R`). No package code is modified; the
circular-plot function is copied into `circular_plot_house.R` and modified
there only.

## Figure design (for the captions)

- **Circular plots** (`circular_plot_house.R`, a copy of
  `plot_freq_of_identified_terms_circular_equations_as_sectors_refined`): one
  sector per equation (labelled with the derivative), one ring per sweep value,
  one bar per selected term whose height is its selection frequency (a
  frequency of 1 reaches 35 % of the ring thickness, the rest is label room);
  dark fill = true term, light tint of the same hue = spurious term; hue per
  method (Bayesian-ARGOS blue, ARGOS red, SINDy green). Rings are labelled in a
  column at the top; the ring order puts the hardest condition outside
  (smallest n or lowest SNR in the outermost ring, 10^5 or the highest SNR or
  ∞ innermost), because the crowded cells need the longest arcs. Labels are
  drawn for terms whose frequency reaches 0.1; in a cell where more labels are
  requested than fit at 5 pt, the labelled bars are spread over the cell and
  the most frequent terms (true terms first) keep their labels; the number
  dropped per cell is in `results/fig6_circular_label_stats.csv` and
  `results/si_circular_label_stats.csv` (only SINDy cells at small n or low
  SNR are affected, e.g. Lorenz at 10^2.8: 14 to 18 of 38 to 48 candidate
  labels). Text: sans-serif (Helvetica, Symbol for ∞), term labels 5 pt, ring
  labels 6.5 pt, headings 7 pt at the stated final width; no radial axis, no
  grid beyond the ring boundaries.
- **Sweep panels** (`panel_helpers.R`): the spurious-term count is drawn on a
  pseudo-logarithmic axis (asinh(y)/ln 10: linear near 0, logarithmic above
  about 1, labelled 0, 1, 2, 5, 10, 20, 50) because SINDy reaches 20 to 50
  spurious terms per equation at low SNR or small n while ARGOS and
  Bayesian-ARGOS stay between 0 and 2; the other two quantities use linear
  axes from 0. SNR panels follow the house axis (10 to 60 dB, ∞ beyond the
  break; the axis starts at 6 dB, so the 1 to 5 dB points are not drawn).
- **Fig. 6**: block a, 2 x 3 panels, rows = sweep (SNR = 49 dB; n = 5000),
  columns = spurious terms per equation, true terms missed per equation, mean
  selection frequency of the spurious terms; curve = mean over the seven
  systems, ribbon = minimum to maximum over the systems; tags a to f. Block b,
  Lorenz, n sweep, rings 10^2.8, 10^3, 10^3.5, 10^4, 10^5, the three methods
  side by side, tags g to i. Composite 183 x 183 mm.
- **Supplementary term-selection figures**: 3 rows (Bayesian-ARGOS, ARGOS,
  SINDy) x 2 columns (n sweep, rings as in Fig. 6b; SNR sweep, rings as the
  submitted SI figures used them, see `RINGS_SNR` in `common.R`), tags a to f
  row-major, 76 mm plots, 160 x 242 mm (three stacked square plots give an
  aspect near 1.5, so 88 mm plots would exceed the 247 mm page depth; the size
  is `PLOT_MM` in the script).
- **Supplementary spurious terms by system**: 7 rows (systems) x 2 columns
  (sweeps) of spurious terms per equation for the three methods, one shared y
  axis, 183 mm wide, aspect 1.3.

## Scripts

| Script | What it does | Runtime |
|---|---|---|
| `common.R` | paths, the four loader modules in separate environments, systems with their true terms in both naming schemes, methods, ring values, the per-trial loader | (sourced) |
| `spurious_terms_table.R` | reads the 42 system x method x sweep result sets; writes `results/term_frequencies.csv`, `results/spurious_summary.csv`, `results/success_rate_check.csv` | 24 s |
| `circular_plot_house.R` | the house copy of the circular plot and `draw_circular_pdf()` | (sourced) |
| `panel_helpers.R` | the sweep panels (house style, pseudo-log axis, SNR axis) | (sourced) |
| `tex_compose.R` | writes and compiles the standalone TikZ composites (helvet; tags 8 pt bold, labels 7 pt) | (sourced) |
| `make_fig6.R` | `figures/fig6_block_a.pdf` (18 x 10.8 in canvas), `figures/circ_lorenz_n_{ba,argos,sindy}.pdf` (12 x 12 in), `figures/shrinkage_spurious_terms.tex/.pdf`; copies the PDF to `imgs/` and `reply-to-reviewers/Figures/`; `results/spurious_summary_aggregate.csv`, `results/fig6_circular_label_stats.csv` | 3 s |
| `make_si_term_selection.R` | 42 circular plots `figures/circ_<system>_<sweep>_<method>.pdf` and `figures/term_selection_<system>.tex/.pdf`; copies to `si-plots/` and `reply-to-reviewers/Figures/`; `results/si_circular_label_stats.csv` | 18 s |
| `make_si_spurious_by_system.R` | `figures/spurious_terms_by_system.pdf` (18 x 23.4 in canvas); copies to `si-plots/` and `reply-to-reviewers/Figures/` | 2 s |
| `run_all.sh` | runs the four scripts in order | 50 s |

```bash
cd ResultsAnalysis/rebuttal-spurious-terms && ./run_all.sh
```

Requirements: R >= 4.3 with dplyr, tidyr, stringr, ggplot2 >= 3.5, scales,
patchwork, circlize (the notebooks' modules load it at source time), latex2exp,
ggh4x, plus `../plot-style/house_style.R`; pdflatex with standalone, tikz and
helvet; pdfinfo.

## Outputs

`results/term_frequencies.csv`, `results/spurious_summary.csv`,
`results/spurious_summary_aggregate.csv`, `results/success_rate_check.csv`,
`results/fig6_circular_label_stats.csv`, `results/si_circular_label_stats.csv`;
`figures/shrinkage_spurious_terms.pdf` (revised Fig. 6, also in `imgs/`),
`figures/term_selection_<system>.pdf` (seven Supplementary figures, also in
`si-plots/`), `figures/spurious_terms_by_system.pdf` (Supplementary figure,
also in `si-plots/`), all three also in `reply-to-reviewers/Figures/`; the
generated `.tex` composites and the panel PDFs beside them.
