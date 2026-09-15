# Weak-form check (Referee 1, Points 6 and 7)

> **Status:** DONE 2026-09-15 (20 trials per case; Dadras + Aizawa in one run of
> 23 min on 16 cores, Rössler re-run separately after a missing-coefficient fix).
> Numbers below are from `results/weak_form_summary.csv`. Not yet cited in the
> reply drafts; proposed wording is in
> `notes/rebuttal-memos/reviewer1/points6-7-weak-form-check.md` (manuscript repo).

## The question

The replies to Points 6 and 7 trace the remaining benchmark failures to the
Savitzky–Golay (SG) derivative estimate: in the noiseless limit its
deterministic truncation error is the whole residual, the residual variance
collapses, and the credible-interval rule declares spurious terms significant
(Aizawa, Rössler, Sprott at SNR = ∞); for the fast Dadras system the same error
carries a systematic component (+0.005 on the constant term) that becomes
significant once n ≳ 6 × 10⁴ shrinks the interval below it. A referee can
reasonably ask: *why differentiate at all?* Weak-form / integral formulations
(Schaeffer & McCalla 2017; Reinbold, Gurevich & Grigoriev 2020; Messenger &
Bortz 2021) avoid pointwise differentiation and are the standard remedy.

This folder answers by substitution. The SG step of `build_design_matrix` is
replaced by a weak-form design and **everything downstream is the pipeline's
own code, unchanged**: the two adaptive-lasso passes (`alasso` from
`R/argos_files.R`, reproduced by `screen()` in
`../rebuttal-diagnostics/common.R`), `rstanarm::stan_glm` with the default
priors and the 90 % posterior-interval rule (`bayes_stage()`). Data come from
`DataGeneration/ode_auto.py` through reticulate, with the initial conditions
and noise seeds of `dadras_psis_loo.R` (Point 7) and
`exp2_heteroscedasticity.R` (Point 6), so the `sg` arm reproduces those runs.

## Weak form used

For ẋ = Θ(x)β and a test function φ_k supported on [t_k − h, t_k + h] with
φ_k = 0 at the ends, integration by parts gives

    ∫ φ_k ẋ dt = −∫ φ_k′ x dt = ( ∫ φ_k Θ(x) dt ) β ,

so row k of the regression is b_k = −∫ φ_k′ x dt against G_k = ∫ φ_k Θ(x) dt.
φ(t) = (1 − ((t − t_k)/h)²)^p on 2m + 1 samples (h = m Δt), normalised so that
∫ φ dt = 1 (the intercept column stays 1 and coefficients keep their scale);
trapezoid quadrature on the samples; query points every m samples (50 %
overlap). This is the WSINDy construction of Messenger & Bortz (2021) with a
fixed test-function family rather than their spectral selection of (m, p).

Arms per trial (same data, same SG parameters):

| arm | states the library is evaluated on | response |
|---|---|---|
| `sg` | SG-smoothed | SG derivative (the pipeline as submitted) |
| `weak/raw/m={1,2}w/p={4,8}` | raw noisy samples (the WSINDy design) | −∫ φ′ x dt on the raw samples |
| `weak/smoothed/m={1,2}w/p={4,8}` | SG-smoothed (same smoothing as `sg`) | −∫ φ′ x_smoothed dt (isolates the derivative step) |

w is the SG window length (in samples) the pipeline's own cross-validated grid
search selects for the response coordinate, so `m = 1w` gives the weak form a
support twice as wide as the SG stencil and `m = 2w` four times.

## Cases

| case | equation | SNR | n | trials | seeds |
|---|---|---|---|---|---|
| Dadras (Point 7) | x₁-dot | 49 dB | 10⁴, 10⁵, 10^5.5 | 20 (benchmark pool rows 1–20) | as `dadras_psis_loo.R` |
| Aizawa (Point 6) | x₃-dot | 60 dB, ∞ | 5000 | 20 (pool rows 1–20) | as `exp2_heteroscedasticity.R` |
| Rössler (Point 6, secondary) | x₃-dot | 61 dB, ∞ | 5000 | 20 | as `exp2_heteroscedasticity.R` |

n = 10^5.5 lies beyond the benchmark grid (10⁴–10⁵) and is included only to
show where the SG failure rate is heading.

## Departures from the manuscript pipeline

* The weak-form arms replace the SG derivative (and, in the `raw` arms, the SG
  smoothing) by the construction above. Nothing else changes.
* The `sg` arm is the pipeline; it re-runs the Point 6 / Point 7 cases on the
  same seeds rather than reading their result files, so every arm sees
  identical data.
* Stan and screening are run per trial with one core per forked worker
  (`ARGOS_NCORES`, default min(16, cores − 4)).

## Scripts

| script | what it does | runtime (16 cores, DGX Spark) |
|---|---|---|
| `weak_form_check.R [trials=20] [cases]` | all arms for all cases; writes `results/weak_form_trials.csv` (one row per case × trial × arm), `results/weak_form_summary.csv`, `results/weak_form_settings.csv` | see status line |
| `make_figures_weak_form.R [primary arm]` | `figures/fig_weak_form_check.pdf` (house style) | seconds |
| `run_all.sh` | both, in order | |

```bash
cd ResultsAnalysis/rebuttal-weak-form
ARGOS_NCORES=16 Rscript weak_form_check.R 20
Rscript make_figures_weak_form.R
```

Per row the CSV records the screened and kept supports, exact recovery, the
OLS intercept on the true support with its s.e. and t value (the Point 7 bias),
the residual sd and lag-1 autocorrelation on the true support, the
deterministic error of the response computed on the *clean* signal (SG
truncation error for `sg`; quadrature error for the weak arms), the number of
weak-form rows K, and the Stan time.

## Headline findings

Primary weak-form setting quoted below: `weak/raw/m=2w/p=8` (raw samples,
support 4w + 1 samples, test-function degree 8 — the arm with the smallest
quadrature error). The other seven arms give the same conclusions; ranges are
given where they differ.

**Point 7 — Dadras x₁-dot, 49 dB (the intercept bias is a property of the SG
derivative and disappears under the weak form).**

| | n = 10⁴ | n = 10⁵ | n = 10^5.5 |
|---|---|---|---|
| SG: OLS intercept on true support, median (t) | +0.0050 (1.2) | +0.0043 (3.1) | +0.0044 (5.8) |
| SG: positive in | 19/20 | 20/20 | 20/20 |
| SG: spurious intercept kept / exact recovery | 0/20, 20/20 | 6/20, 14/20 | 20/20, 0/20 |
| weak: OLS intercept, median (t) | +0.0004 (0.3) | −0.00003 (−0.07) | +0.00004 (0.2) |
| weak: positive in | 13/20 | 9/20 | 13/20 |
| weak: spurious intercept kept / exact recovery | 0/20, 20/20 | 0/20, 20/20 | 0/20, 20/20 |
| residual sd, true support: SG / weak | 0.35 / 0.032 | 0.42 / 0.032 | 0.43 / 0.032 |
| deterministic response error sd (clean signal): SG / weak | 0.29 / 6e-7 | 0.37 / 6e-7 | 0.38 / 6e-7 |
| Stan time per fit: SG / weak | 3 s / 1 s | 26 s / 4 s | 96 s / 8 s |

All eight weak-form arms recover the equation in 20/20 trials at every n; their
intercept medians are within ±0.0001 of zero with no consistent sign. The SG
arm reproduces `dadras_psis_loo.R` (same seeds): 6/20 at 10⁵, median intercept
+0.0043. The SG residual is dominated by the truncation error (sd 0.37 against
noise-only sd ≈ 0.01); the weak-form residual is 13× smaller and its
deterministic part is 10⁻⁶.

**Point 6 — Aizawa x₃-dot, n = 5000 (the weak form does not remove the
noiseless-limit failure; it changes which deterministic error collapses the
intervals).**

| | SG, 60 dB | weak, 60 dB | SG, ∞ | weak, ∞ |
|---|---|---|---|---|
| exact recovery | 14/20 | 13/20 (10–15 across arms) | 0/20 | 6/20 (0–6 across arms) |
| terms kept after Stan, median | 7 | 7 | 22.5 | 8 (8–14.5 across arms) |
| spurious terms kept, median (max) | 0 (24) | 0 (16) | 16 (35) | 1 (15); 1–8 across arms |
| residual sd, true support | 1.5e-2 | 1.9e-3 | 1.0e-5 | 1.9e-7 |
| OLS intercept t-value, true support | 1.7e3 | 3.5e3 | 2.2e6 | 3.6e7 |
| deterministic response error sd | 1.1e-5 | 1.9e-7 | 1.1e-5 | 1.9e-7 |

At 60 dB the weak form gains nothing (the plateau is the library
collinearity, reply to Point 6). At ∞ the residual of the correct model is the
quadrature error of the weak form exactly as it was the truncation error of the
SG derivative; it is 50× smaller, so the credible intervals collapse further
(intercept t-value 4 × 10⁷), and although fewer spurious terms survive (median 1
against 16) exact recovery is restored in at most 6/20 trials. The mechanism is
estimator-agnostic: with no measurement noise, whatever deterministic error the
discretisation leaves becomes "significant" under an interval rule.

**Rössler x₃-dot, n = 5000, 61 dB and ∞ (secondary, a control).** Every arm
recovers the equation in 20/20 trials at both SNRs; this equation does not fail
in the noiseless limit under the R pipeline, so it is uninformative for the
question and is kept only as a check that the weak-form design behaves
normally on a second system (residual sd on the true support: SG 5.2e-2 →
1.9e-3 from 61 dB to ∞; weak 8.5e-3 → 6.8e-7).

**Caveats.** The weak form is applied with a fixed test-function family and
50 % overlap between query windows (K = n/m rows; residuals of neighbouring
rows are correlated, ρ₁ from −0.4 to +0.1 depending on p, against +0.68–0.78
for the SG residual). No tuning of (m, p) beyond the 2 × 2 grid was done. The
result is a demonstration that the regression stage accepts an integral design
unchanged, not a benchmark of weak-form identification.
