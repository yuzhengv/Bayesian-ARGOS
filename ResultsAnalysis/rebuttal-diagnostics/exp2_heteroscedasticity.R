# Experiment 2 (R) -- quantify the residual structure behind Fig. 3d / Fig. 4c
# with the manuscript's own pipeline.
#
# For every trial and SNR (n = 5000) the ground-truth support and the support
# selected by the two alasso passes are each fitted with the package's Bayesian
# stage (rstanarm::stan_glm, default priors, 90 % posterior-interval keep rule,
# as in bayesian_alasso_ro); residuals are taken about the posterior-mean fit.
# Reported per fit: Breusch-Pagan (auxiliary regressor = fitted value), a
# White-type test (fitted + fitted^2), the max/min residual variance over
# deciles of the fitted value, the lag-1 autocorrelation and Durbin-Watson in
# time order.  The final model after the keep rule is recorded too, so the
# per-trial success at n = 5000 is available for comparison with the stored
# benchmark results.  The residual e of the true
# support is decomposed into
#   d       = xdot_SG - f(x_clean)          derivative-estimation error
#   d_trunc = SG-derivative of the clean signal - f(x_clean)  (deterministic part of d)
#   g       = f(x_clean) - Theta(x_SG) b*    state-smoothing error
#   h       = Theta(x_SG) (b* - b_hat)       coefficient-estimation error
#
# Systems: Aizawa x3-dot (Fig. 3d), Rossler x1-dot (Fig. 4c).  Stan runs with
# cores = 1 inside each mclapply worker (4 chains sequentially).
# Outputs: results/exp2_stats.csv, results/exp2_examples.csv (long format).
# Usage: Rscript exp2_heteroscedasticity.R [trials=20]

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
args <- commandArgs(trailingOnly = TRUE)
TRIALS <- if (length(args) >= 1) as.integer(args[1]) else 20L
OUT <- file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "results")
dir.create(OUT, showWarnings = FALSE)

N_FIXED <- 5000
CASES <- list(aizawa = list(eq = 3, snrs = c(10, 20, 30, 40, 45, 50, 55, 57, 60, Inf)),
              rossler = list(eq = 1, snrs = c(20, 30, 40, 50, 55, 58, 61, Inf)))
EXAMPLE_TRIAL <- 0
EXAMPLE_SNRS <- list(aizawa = c(50, 60, Inf), rossler = c(61, Inf))

one <- function(system, trial, snr) {
  eq <- CASES[[system]]$eq
  tr <- SYSTEMS[[system]]$truth[[as.character(eq)]]
  ic <- initial_conditions(system)[trial + 1, ]
  s <- simulate(system, N_FIXED, snr, ic, 20000 + 100 * trial + if (is.infinite(snr)) 99 else snr)
  dm <- library_dm(s$noisy)
  nm <- colnames(dm$sorted_theta)
  y <- dm$xdot_filtered[, eq]
  X <- theta_columns(dm$sorted_theta, nm, tr$terms); colnames(X) <- tr$terms
  bs <- bayes_stage(y, X, tr$intercept)
  st <- residual_stats_from(bs$residuals, bs$fitted)
  e <- st$residuals
  f_clean <- rhs_matrix(system, s$clean)[, eq]
  theta_bstar <- as.numeric(X %*% tr$beta) + if (tr$intercept) tr$beta0 else 0
  d <- y - f_clean; g <- f_clean - theta_bstar; h <- theta_bstar - st$fitted
  sgp <- sg_params(s$noisy[, eq])
  d_trunc <- sg_derivative(s$clean[, eq], sgp["order"], sgp["window"]) - f_clean
  d_noise <- d - d_trunc

  scr <- screen(dm, eq)
  sel <- scr$selected_terms
  bs_sel <- NULL
  if (length(sel) >= 1) {
    Xs <- theta_columns(dm$sorted_theta, nm, sel); colnames(Xs) <- sel
    bs_sel <- bayes_stage(y, Xs, scr$selected_intercept)
  }
  st_sel <- if (!is.null(bs_sel)) residual_stats_from(bs_sel$residuals, bs_sel$fitted) else NULL
  stat_names <- c("sigma", "bp_lm", "bp_p", "bp_r2", "white_lm", "white_p", "var_ratio", "rho1", "durbin_watson")
  row <- data.frame(system = system, eq = eq, trial = trial, snr = snr, n = N_FIXED,
                    sg_order = unname(sgp["order"]), sg_window = unname(sgp["window"]))
  for (k in stat_names) row[[k]] <- st[[k]]
  for (k in stat_names) row[[paste0("sel_", k)]] <- if (is.null(st_sel)) NA else st_sel[[k]]
  row$selected <- paste(to_py_name(sel), collapse = ";")
  row$selected_intercept <- scr$selected_intercept
  row$screen_exact_truth <- setequal(sel, tr$terms) && scr$selected_intercept == tr$intercept
  row$truth_in_selected <- all(tr$terms %in% sel)
  # final model after the posterior-interval rule (the pipeline's output)
  fin_terms <- if (is.null(bs_sel)) character(0) else bs_sel$kept_terms
  fin_int <- if (is.null(bs_sel)) scr$selected_intercept else bs_sel$kept_intercept
  row$final <- paste(to_py_name(fin_terms), collapse = ";")
  row$final_intercept <- fin_int
  row$success <- setequal(fin_terms, tr$terms) && fin_int == tr$intercept
  row$truth_in_final <- all(tr$terms %in% fin_terms) && fin_int == tr$intercept
  row$n_spurious_final <- length(setdiff(fin_terms, tr$terms))
  row$truth_support_all_kept <- setequal(bs$kept_terms, tr$terms) && bs$kept_intercept == tr$intercept
  row$var_e <- var(e); row$var_d <- var(d); row$var_g <- var(g); row$var_h <- var(h)
  row$var_dtrunc <- var(d_trunc); row$var_dnoise <- var(d_noise)
  row$corr_e_d <- cor(e, d); row$corr_e_dtrunc <- cor(e, d_trunc); row$corr_e_g <- cor(e, g)
  row$r2_e_on_d <- cor(e, d)^2
  row$rho1_d <- cor(d[-length(d)], d[-1]); row$rho1_dtrunc <- cor(d_trunc[-length(d_trunc)], d_trunc[-1])
  row$noise_sd_state <- if (is.infinite(snr)) 0 else 10^(-snr / 20) * sqrt(mean((s$clean[, eq] - mean(s$clean[, eq]))^2))

  example <- NULL
  if (trial == EXAMPLE_TRIAL && any(sapply(EXAMPLE_SNRS[[system]], function(v) identical(v, snr)))) {
    arrs <- list(residuals = e, fitted = st$fitted, d = d, d_trunc = d_trunc, g = g, y = y,
                 sel_residuals = if (is.null(st_sel)) rep(NA, length(e)) else st_sel$residuals,
                 sel_fitted = if (is.null(st_sel)) rep(NA, length(e)) else st_sel$fitted)
    example <- do.call(rbind, lapply(names(arrs), function(k)
      data.frame(system = system, snr = snr_label(snr), key = k, index = seq_along(arrs[[k]]) - 1, value = arrs[[k]])))
  }
  list(row = row, example = example)
}

tasks <- list()
for (system in names(CASES)) for (trial in seq_len(TRIALS) - 1) for (snr in CASES[[system]]$snrs)
  tasks[[length(tasks) + 1]] <- list(system, trial, snr)
cat(length(tasks), "conditions on", NCORES, "cores\n")
t0 <- Sys.time()
res <- mclapply(tasks, function(t) one(t[[1]], t[[2]], t[[3]]), mc.cores = NCORES)
ok <- !sapply(res, inherits, "try-error")
write.csv(do.call(rbind, lapply(res[ok], `[[`, "row")), file.path(OUT, "exp2_stats.csv"), row.names = FALSE)
ex <- do.call(rbind, Filter(Negate(is.null), lapply(res[ok], `[[`, "example")))
write.csv(ex, file.path(OUT, "exp2_examples.csv"), row.names = FALSE)
cat("done in", round(as.numeric(Sys.time() - t0, units = "secs")), "s;", sum(!ok), "failures\n")
