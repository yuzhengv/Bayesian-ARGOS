# Weak-form (Integration) arms of the 2 × 2 derivative-vs-integral comparison

> **Status:** READY TO SUBMIT on Hamilton (2026-09-15). Drivers and SLURM
> scripts written and smoke-tested locally (Dadras SNR sweep, Rössler n sweep,
> two trials each, both arms; evaluation loader verified). No Hamilton run yet.
> Requires `pysindy >= 1.7.3` in the `pysindy-prod` conda environment on
> Hamilton (tested with 2.1.0 locally); check with
> `python -c "import pysindy; print(pysindy.__version__)"` before submitting.

## The comparison

| | Savitzky–Golay derivative (SG) | weak form (Integration) |
|---|---|---|
| **Bayesian-ARGOS** | `Experiments/bayesian-alasso-ro/<system>/results` (stored benchmark, Figs 2–3) | **this folder**: `Experiments/bayesian-alasso-ro-weak/<system>/results` |
| **SINDy** | `Pysindy/exp/<system>/results` (stored benchmark, STLSQ 0.005) | `Pysindy/weak/exp/<system>/results` |

Response construction (differentiate vs integrate) and selection stage
(STLSQ vs adaptive-lasso screening + Stan) are separate factors, so the
result answers both "why not integrate?" (Referee 1, Points 6–7 context) and
"is the Bayesian-ARGOS gain a property of the selection stage or of the
derivative?" (novelty, both referees). Only the two Integration arms are new;
the SG arms are the stored 100-trial results.

## What the Integration arms do

`weak_form_common.R` (shared by both arms) replaces the SG step of
`build_design_matrix` by PySINDy's `WeakPDELibrary` (Messenger & Bortz 2021):
row k of the regression is b_k = −∫ φ_k′ x dt against G_k = ∫ φ_k Θ(x) dt,
with test functions (1 − s²)^p on K subdomains of half-width H_xt. Settings
(environment variables, defaults in brackets): `WEAK_H` [0.25 time units =
support of 50 samples at dt 0.01], `WEAK_P` [8], `WEAK_KFRAC` [0.25: K =
⌈0.25 n⌉ subdomains, so the number of rows scales with n]. Subdomain centres
are drawn at random by PySINDy; NumPy is seeded per (condition, trial) with
the same seed in both arms, so both arms use the same subdomains.

* **Bayesian-ARGOS (Integration)** (`weak_n.R`, `weak_snr.R`): the weak design is
  rescaled row-wise by ∫ φ_k dt (constant column = 1, dropped; coefficients keep
  their scale), columns renamed to the pipeline's monomial names and ordered by
  degree, and passed to `bayesian_alasso_ro()` **unchanged** (two adaptive-lasso
  passes, `stan_glm`, 90 % posterior-interval rule). Results have the same
  structure and file names as `Experiments/bayesian-alasso-ro/<system>/results`.
* **SINDy (Integration)** (`Pysindy/weak/weak_n.R`, `weak_snr.R`): `ps.SINDy(
  optimizer = STLSQ(threshold = 0.005), feature_library = WeakPDELibrary(...))`
  fitted on the raw noisy samples (the weak library computes its own response).
  Same threshold as the SG baseline (`STLSQ_THRESHOLD` overrides). Results have
  the same structure and file names as `Pysindy/exp/<system>/results`.

Data, initial conditions (seed 100, same `runif` sequences as the original
drivers), grids (n = 10² … 10⁵ at 49 dB; SNR 1 … 61 dB and ∞ at n = 5000),
100 trials, library degree 5: all as in the stored benchmarks.

## Departures from the manuscript pipeline

* Response and library are integrated against test functions instead of being
  differentiated and smoothed; nothing else changes.
* Rows of the weak design overlap (K = n/4 centres, support 50 samples), so
  their residuals are correlated; the Gaussian likelihood ignores this, as it
  ignores the lag-1 autocorrelation ≈ 0.7 of the SG residual.
* A degree-5 library on raw noisy samples carries a small even-power bias at
  low SNR (E[x²] = x² + σ²); this is the standard weak-form design and is left
  as is.
* The STLSQ threshold (0.005) was chosen for the SG response; it is kept fixed
  for comparability (the coefficient scale is unchanged by the weak form).
* Noise realisations are seeded (the stored SG runs were not), so the two arms
  do not share noise draws with the stored results.

## Running on Hamilton

```bash
cd Experiments
bash submit_all_weak_form.sh                 # all systems, both arms, both sweeps
bash submit_all_weak_form.sh dadras          # one system
ARMS=sindy SWEEPS=snr bash submit_all_weak_form.sh
DRY_RUN=1 bash submit_all_weak_form.sh       # show what would be submitted
```

Each system folder holds `<system>_{n,snr}.sh` (SLURM job, 24 cores for the
Bayesian arm, 12 for SINDy, `source activate "${WEAK_CONDA_ENV:-pysindy-prod}"`, so a rebuilt environment can be selected with `WEAK_CONDA_ENV=<name>` at submission) and
`<system>_{n,snr}_submit.sh` (one job per equation and per n-decade / SNR
block, as the original scripts). The drivers locate the repository through
`ARGOS_ROOT`, else `/nobackup/qtzk83/Projects/Bayesian-ARGOS`, else relative to
the job folder; `ARGOS_PYTHON` overrides the Python reticulate uses.

Local runtimes (DGX Spark, 8 cores): Bayesian arm 1.7 s per Stan fit at
n = 5000 (K = 1250 rows), SINDy arm 0.1 s. On Hamilton the Bayesian
SNR sweep is 62 × 100 × 3 fits per system; the n sweep 31 × 100 × 3 with
K up to 25 000 rows at n = 10⁵.

## Analysis

`ResultsAnalysis/rebuttal-weak-form/success_2x2.R` builds the success-rate
tables for the four arms (both sweeps, three systems) with the paper's own
evaluation functions (`MethodsEvaluation/`, `Pysindy/utils/`), writing
`results/success_2x2_{n,snr}.csv`. The evaluation code dispatches on a fixed
list of method names, so the weak arms are loaded under the label of their SG
counterpart and relabelled; no evaluation code is modified.

| script | role |
|---|---|
| `HAMILTON_ENV_CHECK.md` | prompt for checking the Hamilton environment and running a one-minute smoke test before submitting |
| `weak_form_common.R` | root/Python setup, system definitions, weak design, one-fit wrappers for both arms |
| `weak_n.R`, `weak_snr.R` | Bayesian-ARGOS (Integration) drivers (generic, `SYSTEM` from the job script) |
| `<system>/<system>_{n,snr}.sh`, `_submit.sh` | SLURM job and submission scripts |
| `../submit_all_weak_form.sh` | submits both arms for the chosen systems |
| `../../Pysindy/weak/weak_{n,snr}.R`, `exp/<system>/*.sh` | SINDy (Integration) drivers and jobs |
