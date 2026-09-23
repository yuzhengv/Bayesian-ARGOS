# Rebuttal ablation summary (Referee 1, Point 1a/1b)

**Status (2026-09-21): complete.** Numeric summary of the screening-stage
ablation, exported from the stored 100-trial results; used for the reply to
Referee 1 Point 1 (the ablation studies the referee and the editor asked for)
and cross-referenced from Referee 2 Point 1. Internal reference only: the
reply quotes these numbers in prose, the SI shows the figures.

## Referee's question

> "the two-pass frequentist stage is explained in detail including ad hoc
> weight estimates, but is never really motivated - is it really necessary and
> does it improve screening?" (1b); "it is completely unclear which steps are
> truly necessary or novel" (1a).

## What the ablation is

Six variants of the screening stage, each followed by the unchanged Bayesian
stage (`rstanarm::stan_glm`, 90 % posterior-interval rule), run with the
benchmark sweeps (n = 10^2 … 10^5 in 31 steps at 49 dB; SNR = 1 … 61 dB in 1 dB
steps plus the noiseless case at n = 5000; 100 initial conditions per grid
point) on the Lorenz, Rössler and Thomas systems. Arms and folders:

| Arm | Pass-1 pilot | Pass-2 pilot | Folder |
|---|---|---|---|
| Ridge–OLS (baseline, the manuscript's pipeline) | ridge | OLS | `Experiments/bayesian-alasso-ro` |
| OLS–Ridge | OLS | ridge | `Experiments/bayesian-alasso-or` |
| OLS–OLS | OLS | OLS | `Experiments/bayesian-alasso-oo` |
| Ridge–Ridge | ridge | ridge | `Experiments/bayesian-alasso-rr` |
| Single OLS | — (no first pass, no library shrinkage) | OLS | `Experiments/bayesian-alasso-single-ols` |
| Single Ridge | — | ridge | `Experiments/bayesian-alasso-single-ridge` |

Arm code: `R/bayesian_alasso_abl.R`; submission: `Experiments/submit_all_ablation.sh`;
figures: `ResultsAnalysis/success-rate-plots/imgs/{lorenz,rossler,thomas}_ablation_success_rate_plot.pdf`
(notebooks `*_ablation_results_compare.ipynb`). No arm without the Bayesian
stage was run.

## Departures from the manuscript pipeline

None. The script only re-reads the stored results with the notebooks' own
evaluation function (`generate_total_success_rate_table` in
`MethodsEvaluation/results_processing_for_analysis.R`, unchanged); success
means all three equations recovered exactly, as in the manuscript.

## Scripts

| Script | What it does | Runtime |
|---|---|---|
| `ablation_success_table.R` | reads the 36 arm × system × sweep result sets; writes `results/ablation_success_rates.csv` (long: system, sweep, arm, x, success; x = log10 n or SNR in dB with Inf for the noiseless case) and `tables/ablation_summary.md` (80 % crossing, plateau, noiseless point, success at reference points, difference to the baseline per arm) | ≈ 1 min |
| `run_all.sh` | runs the script from any working directory | |

```bash
cd ResultsAnalysis/rebuttal-ablation && ./run_all.sh
```

Requirements: R ≥ 4.3 with dplyr, tidyr, stringr, plus the packages the
processing script loads.

## Outputs

`results/ablation_success_rates.csv`, `tables/ablation_summary.md`.
