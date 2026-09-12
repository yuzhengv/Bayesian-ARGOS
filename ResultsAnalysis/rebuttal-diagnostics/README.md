# Rebuttal diagnostics (R pipeline)

Evidence for the reply to Referee 1 Point 6 and Referee 2 minor point 2 of
COMMSPHYS-26-0738-T, computed with the R implementation that produced the
manuscript's benchmarks (`R/argos_files.R`, `R/bayesian_alasso_ro.R`). A
Python counterpart with the same output schema lives in
`pyargos_with_sindy_shred/pyargos/results-diagnostics/rebuttal-diagnostics/`;
its `make_figures.py --res results --tables-only` writes `tables/summary_tables.md`,
while the figures are drawn here by `make_figures.R`.

| Script | What it does |
|---|---|
| `common.R` | systems in the `ode_auto.py` format; data generation through the benchmark code path (`DataGeneration/ode_auto.py`, scipy odeint via reticulate, NumPy seeded before the noisy call); `build_design_matrix`; the two `alasso` passes; the Bayesian stage of `bayesian_alasso_ro` (`rstanarm::stan_glm`, 90 % posterior-interval rule) as `bayes_stage()`; VIF, Belsley diagnostics, residual statistics |
| `exp1_collinearity.R` | per-term VIF and condition number on the full / degree-trimmed / pass-2-selected / ground-truth designs over the Fig. 3b grids; same monomials on one trajectory, a uniform box, 50 short transients; clean trajectory vs n; Belsley near-dependency groups |
| `exp2_heteroscedasticity.R` | Breusch–Pagan, White-type, decile variance ratio, lag-1 autocorrelation, Durbin–Watson vs SNR for the ground-truth support and for the support the screening selects, both fitted with the pipeline's Stan stage (residuals about the posterior mean); final model after the posterior-interval rule and per-trial success; residual decomposition into derivative-estimation (incl. SG truncation), smoothing and coefficient error |
| `aizawa_screening_vs_hmc.R` | re-reads the stored 100-trial Aizawa results: pass-2 support size, survival of spurious candidates through the 90 % credible-interval rule, success rate; exports per-trial supports (`aizawa_eq3_supports.csv`) |
| `exp6_lambda_rule.R` | tests whether tuning the adaptive lasso by `lambda.1se` instead of `lambda.min` (pass 2, or both passes) removes the large-n decline; the Bayesian stage is the pipeline's own (`stan_glm`, posterior-interval rule), so the reported success is the pipeline's |
| `aizawa_success_vs_snr_all_methods.R` | success vs SNR for Bayesian-ARGOS, ARGOS and SINDy from the stored 100-trial results (the paper's own evaluation functions), plus the ARGOS x3-dot failure decomposition (missing vs spurious terms, bootstrap-active chain neighbours) |
| `dadras_psis_loo.R` | Referee 1 Point 7: PSIS-LOO Pareto k-hat for Dadras x1-dot over the first 20 initial conditions of the benchmark pool at n = 1e4, 1e4.4, 1e4.6, 1e4.8, 1e5 with the pipeline's own screening and Stan stage (k-hat from the pointwise Gaussian log-likelihood with `loo::psis` + `relative_eff`, chunked); where the flagged points sit in time and on the attractor (leverage, residual, derivative-estimation error); refit without them; intercept bias vs its shrinking interval (`results/dadras_loo_*.csv`; `results/dadras_stored_intercept_vs_n.csv` from the stored 100-trial results; `results/dadras_intercept_bias_source.csv`) |
| `make_figures_dadras_loo.R` | the Point 7 figure (`figures/fig_dadras_psis_loo.pdf`): k-hat vs time at two n on one scale, flagged points on the attractor, intercept vs 90 % half-width vs n |
| `make_figures.R` | the four figures in the manuscript's house style (theme, palette, fonts and broken ∞-axis of the success-rate notebooks), following Referee 1's figure requirements: only the true terms highlighted in the VIF panel, ∞ symbol and 10-dB ticks on every SNR axis, sans-serif with two axis borders and no grid, compared subplots on one vertical scale (standardised residuals) and labelled by system and condition |
| `run_all.sh` | runs everything, then `make_figures.R` and the Python `make_figures.py --tables-only` |

Requirements: R ≥ 4.3 with tidyverse, glmnet, signal, Metrics, lmtest,
reticulate, rstan, rstanarm and loo, plus a Python with numpy and scipy for
`ode_auto.py` (set `ARGOS_PYTHON`; default is the `env_dgx_spark` conda
interpreter). `ARGOS_NCORES` overrides the number of forked workers (default
min(16, cores − 4)); Stan runs with one core per worker.

Runtime on a 20-core machine: exp1 ≈ 5 min; exp2 and exp6 are dominated by
the `stan_glm` fits (roughly 30 min and 1–2 h); dadras_psis_loo ≈ 17 min.

Outputs: `results/*.csv|json`, `figures/*.pdf`, `tables/summary_tables.md`.
