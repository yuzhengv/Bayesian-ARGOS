# Rebuttal diagnostics (R pipeline)

> **Status:** DONE 2026-10-24 (exp1 re-run with 20 trials on both sweeps, the
> $n$ sweep having had 10; figures for the revised manuscript). The
> experiments (exp1, exp2, exp6, the Aizawa and Rossler stored-results reads
> and the Dadras PSIS-LOO run) are complete; the figure scripts of the plotting
> standard (`reply-to-reviewers/Notes/plotting-standard.md`, blocks M1, M4, M5
> and S-D) draw from `results/` with `../plot-style/house_style.R` and write
> the manuscript files `imgs/bayesian_argos_schematic.pdf` (Fig. 1),
> `imgs/diagnostics_aizawa_chains.pdf` (Fig. 4),
> `imgs/diagnostics_derivative_error.pdf` (Fig. 5) and
> `si-plots/residual_statistics.pdf` (SI), with copies under
> `reply-to-reviewers/Figures/`. The older `fig_*.pdf` files in `figures/`
> are the reply drafts and are kept unchanged. The sampling comparison
> (one trajectory against transients and box samples) stays internal and
> appears in no manuscript figure.

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
| `exp1_collinearity.R` | per-term VIF and condition number on the full / degree-trimmed / pass-2-selected / ground-truth designs over the Fig. 3b grids, 20 trials on each sweep; same monomials on one trajectory, a uniform box, 50 short transients; clean trajectory vs n; Belsley near-dependency groups |
| `table_s1_nested_vif_49db.R` | Supplementary Table S1 at the benchmark convention: VIFs of the seven true terms and condition number on the four nested designs at n = 5000, SNR = 49 dB for the same 20 trajectories as exp1 (exp1's `one_condition()` and `designs_for()` reused by evaluating only the function definitions of `exp1_collinearity.R`); writes `results/table_s1_{vif,kappa,meta}_49db.csv` and `tables/table_s1_49db.md` (LaTeX rows, support sizes); 13 s on 16 cores |
| `exp2_heteroscedasticity.R` | Breusch–Pagan, White-type, decile variance ratio, lag-1 autocorrelation, Durbin–Watson vs SNR for the ground-truth support and for the support the screening selects, both fitted with the pipeline's Stan stage (residuals about the posterior mean); final model after the posterior-interval rule and per-trial success; residual decomposition into derivative-estimation (incl. SG truncation), smoothing and coefficient error; SNR grids 10, 20, 30, 40, 45, 50, 55, 57, 60 dB and ∞ for Aizawa and the same plus 58 and 61 dB for Rössler (the Rössler points 10, 45, 57 and 60 dB added 2026-11-07 with `Rscript exp2_heteroscedasticity.R 20 missing`, which runs only the conditions absent from `results/exp2_stats.csv` and appends them, stored rows and `exp2_examples.csv` untouched) |
| `aizawa_screening_vs_hmc.R` | re-reads the stored 100-trial Aizawa results: pass-2 support size, survival of spurious candidates through the 90 % credible-interval rule, success rate; exports per-trial supports (`aizawa_eq3_supports.csv`) |
| `exp6_lambda_rule.R` | tests whether tuning the adaptive lasso by `lambda.1se` instead of `lambda.min` (pass 2, or both passes) removes the large-n decline; the Bayesian stage is the pipeline's own (`stan_glm`, posterior-interval rule), so the reported success is the pipeline's |
| `aizawa_success_vs_snr_all_methods.R` | success vs SNR for Bayesian-ARGOS, ARGOS and SINDy from the stored 100-trial results (the paper's own evaluation functions), plus the ARGOS x3-dot failure decomposition (missing vs spurious terms, bootstrap-active chain neighbours) |
| `dadras_psis_loo.R` | Referee 1 Point 7: PSIS-LOO Pareto k-hat for Dadras x1-dot over the first 20 initial conditions of the benchmark pool at n = 1e4, 1e4.4, 1e4.6, 1e4.8, 1e5 with the pipeline's own screening and Stan stage (k-hat from the pointwise Gaussian log-likelihood with `loo::psis` + `relative_eff`, chunked); where the flagged points sit in time and on the attractor (leverage, residual, derivative-estimation error); refit without them; intercept bias vs its shrinking interval (`results/dadras_loo_*.csv`; `results/dadras_stored_intercept_vs_n.csv` from the stored 100-trial results) |
| `dadras_residual_statistics.R` | residual statistics of the Dadras x1-dot fit on the true support over the first 20 initial conditions of the benchmark pool (the 49 dB trajectories of `dadras_psis_loo.R`, same noise seeds) at n = 1e4 and 1e5, SNR = 49 dB and ∞, with the pipeline's Stan stage as in exp2: Breusch-Pagan R^2, decile variance ratio, lag-1 autocorrelation, residual sd against the noise sd on x1, the decomposition with the SG truncation error, the 99.9 % quantile and maximum of |d|, the OLS intercept of y on the true terms and of the truncation error on the true terms at the clean states, and the OLS residual statistics alongside (identical to three digits); `results/dadras_residual_statistics.csv`, `tables/dadras_residual_statistics.md` (medians and ranges). Added 2026-11-07 (Yuzheng: "Why only 3 trials here?"); supersedes the interactive three-trial spot check of 2026-09-12 kept as `results/dadras_residual_structure.csv` and `results/dadras_intercept_bias_source.csv` (no script); 280 s on 16 cores |
| `make_figures_dadras_loo.R` | the Point 7 figure (`figures/fig_dadras_psis_loo.pdf`): k-hat vs time at two n on one scale, flagged points on the attractor, intercept vs 90 % half-width vs n |
| `make_figures.R` | the four figures in the manuscript's house style (theme, palette, fonts and broken ∞-axis of the success-rate notebooks), following Referee 1's figure requirements: only the true terms highlighted in the VIF panel, ∞ symbol and 10-dB ticks on every SNR axis, sans-serif with two axis borders and no grid, compared subplots on one vertical scale (standardised residuals) and labelled by system and condition |
| `aizawa_attractor_geometry.R` | why the near-dependent chains of the Aizawa x3-dot library differ by a factor of x3: for the 20 initial conditions of the pool (clean, n = 5000) the range, sign and mean frequency of each state, the median |corr(theta, x_k theta)| over the degree-1-to-4 monomials for k = 1, 2, 3, and |corr| of the exchanged pairs and the x3-power chain on the trajectory, on uniform samples of the initial-condition box, and on uniform samples with the attractor's x3 range (`results/aizawa_attractor_geometry.csv`); x3 is positive on about 72 % of the trajectory while x1 and x2 are symmetric about zero, its spread is not narrower than theirs, and uniform samples with the attractor's x3 range reproduce the trajectory correlations, so the chains come from the sign and range of x3, not from its slow variation; about 20 s |
| `rossler_screening_vs_hmc.R` | re-reads the stored 100-trial Rossler x1-dot results over both sweeps (every SNR level including the noiseless case): pass-2 support size, truth admitted, success, survival of spurious candidates through the credible-interval rule, and the per-term FINAL-model selection frequency of every term (`results/rossler_eq1_screening_vs_hmc.csv`, `results/rossler_eq1_term_frequencies.csv`); the same per-term frequencies for Dadras x1-dot (stored results, `results/dadras_eq1_term_frequencies.csv`) and Aizawa x3-dot (from `aizawa_eq3_supports.csv`, `results/aizawa_eq3_term_frequencies.csv`); about 10 s |
| `fig1_insets.R [snr=57] [n=5000]` | the three diagnostic insets of Fig. 1e from ONE successfully identified Aizawa trial (pool rows 1, 2, ... with the exp2 noise seed of the row, all three equations through screening and the Stan stage, until every equation is recovered exactly; pool row 1 succeeds): standardised residuals of x3-dot against the fitted value, Pareto k-hat against time (`loo::psis` as in `dadras_psis_loo.R`), VIF of every term of the design entering the second pass; `results/fig1_insets_{trial,residuals,khat,vif}.csv`, `figures/fig1_inset_{residuals,khat,vif}.pdf` (34 mm final width, aspect 0.95; the three panels aligned with `cowplot::align_plots`, equal plot areas on one level, 2026-11-01; the seven true terms of the VIF inset in the colours of the lasso-path plots of Fig. 1b, `AIZAWA_TRUTH_COL_FIG1` in `fig_helpers.R`, 2026-11-02); about 20 s; `Rscript fig1_insets.R plot-only` redraws the insets from the stored CSVs in a second |
| `compose_fig1_schematic.sh` | overlays the three insets on the submitted drawing `schematic-plots/algorithm-description.pdf` (standalone TikZ, generated `compose_fig1_schematic.tex` kept beside the script: white rounded rectangles inside the three pink frames of panel e, new titles "Residual analysis", "Influential observations", "Multicollinearity" in Helvetica 14 pt on the page, insets scaled to 186 pt) and writes `imgs/bayesian_argos_schematic.pdf` (copy in `figures/`); about 2 s |
| `fig_helpers.R` | helpers shared by the figure scripts below: plotmath term labels, the fixed colours and shapes of the seven true terms of Aizawa x3-dot, the vertical keys drawn as panels (standard 4.7.2), the term-selection "skyline" panels (decision D6), the standardised-residual panel, `snr_axis_log()` |
| `make_fig4_aizawa_chains.R` | revised Fig. 4 (mechanism I): a, b final-model term-selection frequency (one bar per term, one facet per condition: n = 10^3.7, 10^4.4, 10^5 at 49 dB; SNR = 50, 55, 61 dB and infinity at n = 5000), c, d median VIF of every term of the pass-2 design on one log scale, e, f condition number of the four nested designs (median, interquartile range); one key per row at the left (the key of a lists the intercept first, as the bars do; key spacing widened 2026-11-01; the intercept has no VIF, see `vif()` in `common.R`); `figures/diagnostics_aizawa_chains.pdf` -> `imgs/` (183 mm, aspect 1.1); seconds |
| `make_fig5_derivative_error.R` | revised Fig. 5 (mechanism II): a, b term-selection frequency for Rossler x1-dot (61 dB, infinity) and Dadras x1-dot (n = 10^4, 10^5), in the order of the main text, c, d Rossler x1-dot standardised residuals at 61 dB and infinity, e, f Dadras block-maximum k-hat (blocks of 100 observations) against time for trial 6 at n = 10^4 and 10^5 on one scale, g flagged observations of the 20 trials at n = 10^5 on the attractor, h the admitted intercept against n with the 90 % half-width; keys at the left (spacing widened 2026-11-01, the row-4 key in two headed groups, g and h); `figures/diagnostics_derivative_error.pdf` -> `imgs/` (183 mm, aspect 1.15); about 30 s (one clean Dadras trajectory through `ode_auto.py`) |
| `make_si_residual_statistics.R` | SI figure S-D: six panels in the house style, Breusch-Pagan R^2 (a, added 2026-11-07, log axis), decile variance ratio (b), lag-1 autocorrelation (c) and the variance decomposition (d, with the truncation error of the SG derivative as a dashed series) against SNR, both systems on the Aizawa grid 10 to 60 dB and ∞ (the Rössler rows at 58 and 61 dB are not drawn), Aizawa x3-dot residuals at 60 dB (e) and infinity (f), from `results/exp2_stats.csv` and `exp2_examples.csv` (the statistics panels of `make_figures.R`'s `fig_heteroscedasticity.pdf`); `figures/residual_statistics.pdf` -> `si-plots/` (183 mm, aspect 0.7); seconds |
| `run_all.sh` | runs everything, then `make_figures.R` and the Python `make_figures.py --tables-only` |

Requirements: R ≥ 4.3 with tidyverse, glmnet, signal, Metrics, lmtest,
reticulate, rstan, rstanarm and loo, plus a Python with numpy and scipy for
`ode_auto.py` (set `ARGOS_PYTHON`; default is the `env_dgx_spark` conda
interpreter). `ARGOS_NCORES` overrides the number of forked workers (default
min(16, cores − 4)); Stan runs with one core per worker.

Runtime on a 20-core machine: exp1 ≈ 10 min (20 trials on both sweeps); exp2 and exp6 are dominated by
the `stan_glm` fits (roughly 30 min and 1–2 h); dadras_psis_loo ≈ 17 min.

Reproducibility of exp1: the trajectories and noise are seeded, so the full and
true-term designs reproduce exactly between runs. The two `alasso` passes call
`cv.glmnet` without a fold seed (`R/argos_files.R`, package code, unchanged), so
between the 2026-09-09 and 2026-10-24 runs the degree trim differed in about
3 % of the conditions and the screened support in about 20 %. The medians drawn
in the figures and every column of the then Supplementary Table S1 (n = 5000, 50 dB)
were unchanged. Since 2026-11-06 Table S1 is drawn from `table_s1_nested_vif_49db.R`
at 49 dB, the SNR of the benchmark's n sweep, instead of the 50 dB point of exp1's grid.

Figure scripts for the revised manuscript (2026-10-14), in order:

```bash
cd ResultsAnalysis/rebuttal-diagnostics
Rscript rossler_screening_vs_hmc.R        # term frequencies for Figs. 4a and 5a
Rscript fig1_insets.R                     # Fig. 1e insets (one Aizawa trial, 57 dB, n = 5000)
bash compose_fig1_schematic.sh            # imgs/bayesian_argos_schematic.pdf
Rscript make_fig4_aizawa_chains.R         # imgs/diagnostics_aizawa_chains.pdf
Rscript make_fig5_derivative_error.R      # imgs/diagnostics_derivative_error.pdf
Rscript make_si_residual_statistics.R     # si-plots/residual_statistics.pdf
```

They source `../plot-style/house_style.R` (theme, palette, axes, `save_figure()`
on the k = 2.5 canvas, `cairo_pdf` so that the fonts are embedded) and
`fig_helpers.R`. The manuscript directory is found four levels up or given
by `ARGOS_MANUSCRIPT_DIR`; a submitted figure file is never overwritten.

Outputs: `results/*.csv|json`, `figures/*.pdf`, `tables/summary_tables.md`.
