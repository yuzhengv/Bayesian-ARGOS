# Weak-form (Integration) arms of the 2 × 2 derivative-vs-integral comparison

> **Status:** RUN on Hamilton 2026-09-16 (commit 418874e, environment
> `pysindy-weak`, pysindy 2.1.0; all 62 result files present). Evaluated with
> `ResultsAnalysis/rebuttal-weak-form/success_2x2.R` (tables
> `results/success_2x2_{n,snr}.csv`, figure `figures/fig_success_2x2.pdf`).
> Decisions 2026-09-16 (Yuzheng): no re-run. **The evidence for the reply and the
> SI is the Bayesian-ARGOS pair only** (SG step vs weak-form design, everything
> else unchanged): `ResultsAnalysis/rebuttal-weak-form/figures/fig_sg_vs_weak.pdf`.
> The two SINDy arms are kept as run for internal reference
> (`fig_success_2x2_internal.pdf`, tables in `success_2x2_{n,snr}.csv`) and are
> not cited: within that pair the response changes but the fixed threshold
> was set for the SG response, so it is not a clean test of integration, and
> Weak SINDy as published (own threshold selection, generalised least squares)
> is a different method that was not benchmarked. Framing: integration does
> not improve *this pipeline*; making it work would need per-problem settings,
> which the pipeline is designed to avoid. Nothing cited in the reply drafts
> yet; proposed wording in `notes/rebuttal-memos/reviewer1/points6-7-weak-form-check.md`.

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

## Results (2026-09-16)

Success rate = all three equations exactly recovered, 100 trials.

| | Aizawa n=10⁴ | Aizawa 61 dB | Aizawa ∞ | Dadras n=10⁴ | Dadras n=10⁵ | Dadras 49 dB | Dadras ∞ | Rössler n=10³ | Rössler 49 dB | Rössler ∞ |
|---|---|---|---|---|---|---|---|---|---|---|
| Bayesian-ARGOS (SG) | 0.64 | 0.71 | 0.00 | 0.98 | 0.70 | 0.99 | 0.98 | 0.90 | 1.00 | 0.26 |
| Bayesian-ARGOS (Integration) | 0.30 | 0.38 | 0.25 | 0.68 | 0.88 | 0.77 | 1.00 | 0.84 | 0.98 | 1.00 |
| SINDy (SG) | 0.25 | 0.21 | 0.13 | 0.99 | 1.00 | 1.00 | 1.00 | 0.29 | 0.97 | 0.86 |
| SINDy (Integration) | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | 0.69 | 0.60 | 0.00 | 0.54 | 0.55 |

Rows 1–2 are the evidence used (`fig_sg_vs_weak.pdf`); rows 3–4 are internal.
The Integration arm of Bayesian-ARGOS wins exactly where the SG derivative
fails (Dadras at n = 10⁵, Rössler and Aizawa at SNR = ∞) and loses at finite
SNR and moderate n. Two effects, both properties of the weak form rather than
of a particular setting, contribute to the losses and were diagnosed with `ResultsAnalysis/rebuttal-weak-form/diag_overlap.R` and
`diag_conditioning.R`:

1. **Window overlap inflates the effective sample size (Bayesian arm).** With
   K = n/4 subdomains of 50 samples every sample sits in ≈ 12 windows; the
   Bayesian stage treats the rows as independent, so intervals are too narrow
   and a spurious intercept survives at finite SNR (Dadras x₃-dot: 16–30 % of
   trials at 30–61 dB, 0 % at ∞). With K = n/25 (50 % overlap) or n/50 (no
   overlap) the Dadras x₃-dot recovery is 20/20. For Aizawa x₃-dot the fewer
   rows starve the screening instead (0.40 → 0.20 → 0.10 for K = n/4, n/25,
   n/50): at n = 5000 the record holds only ≈ 100 independent windows for a
   55-column library.
2. **The weak transform worsens the conditioning of the polynomial library
   (both arms).** Aizawa, 49 dB, degree-≤4 library: Belsley κ 38 (SG) vs
   104–109 (weak); median VIF of the true x₃-dot terms 87 (SG) vs 93–242
   (weak). Integrating monomials against a wide test function is a low-pass
   filter that makes the chain terms (x₁²x₃, x₁²x₃², …) more alike, so the
   chain-substitution failure of Point 6 becomes more frequent.
3. **SINDy (Integration) keeps the baseline's STLSQ settings, so within the
   SINDy pair the weak form is the only difference.** With the baseline's
   threshold 0.005 (and PySINDy's default ridge penalty) the weak features give
   dense models on Aizawa (50+ terms, 0 % exact at every SNR) and sparse ones on
   Dadras and Rössler. A coarser fixed threshold does not change the Aizawa
   verdict: on 10 trials at 30, 49 and 61 dB, thresholds 0.01, 0.05 and 0.1 give
   0/10 exact (the smallest true coefficient, 0.1 on x₁³x₃, is not separable
   from the weak form's spurious coefficients of 0.01–0.05), whereas Dadras
   reaches 10/10 from 0.01 upwards. PySINDy's optional column normalisation
   (off in the baseline, as in the original SINDy) leaves the models dense
   (53–56 terms with or without it on five Aizawa trials). Weak SINDy's own
   threshold selection (a scan over a λ grid minimising a residual-plus-sparsity
   loss, Messenger & Bortz 2021) is not adopted: it is not implemented in
   PySINDy, it still requires a user-specified λ grid, and it carries no
   statistical justification comparable to the credible-interval rule.

The Dadras and Rössler n-sweep gaps at small n (10² – 10³) have the same
origin as item 1: with a 50-sample window a record of 100–1000 samples holds
2–20 independent windows.
