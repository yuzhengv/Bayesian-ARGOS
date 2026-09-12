# Dadras x1-dot PSIS-LOO check (Referee 1 Point 7) with the manuscript's own
# R pipeline.
#
# Fig. 4b of the submitted manuscript shows Pareto k-hat for one Dadras trial:
# no k > 0.7 at n = 1e4, a cluster of k > 0.7 near the end of the record at
# n = 1e5.  The referee asks what k-hat measures, whether the flagged points
# are a property of the trajectory "that far away" in time or appear
# sporadically, and whether many short trajectories would avoid them.  The
# manuscript also attributes the spurious intercept at large n to those points
# without a test.  This script answers, for several trials and n:
#   (1) where the flagged points sit in time (clustered / sporadic) and in phase
#       space (state, |x2 x3|, leverage, residual, derivative-estimation error);
#   (2) whether trials that keep the spurious intercept have more flagged
#       points than trials that do not, and whether refitting without the
#       flagged points removes the intercept (the causal claim);
#   (3) the competing explanation: the intercept posterior narrows as
#       1/sqrt(n) around a small fixed bias (mean derivative-estimation error),
#       so its 90 % interval eventually excludes zero regardless of any
#       influential point.  Recorded: intercept posterior mean / interval, the
#       OLS intercept and its s.e. on the true support, and mean(d).
#
# Pipeline per trial: benchmark data path (ode_auto.py, SNR 49 dB, dt 0.01,
# the same seeded pool of initial conditions as dadras_n.R), build_design_matrix,
# the two alasso passes (screen), then rstanarm::stan_glm on the selected
# support with the 90 % posterior-interval keep rule (bayes_stage).  k-hat is
# computed from the Gaussian pointwise log-likelihood of the posterior draws
# with loo::psis (r_eff from loo::relative_eff), in chunks of observations so
# that n = 1e5 fits in memory; this is what rstanarm::loo() does internally.
#
# Outputs (results/):
#   dadras_loo_trials.csv    one row per (trial, n): support, intercept kept,
#                            intercept posterior, counts of k > 0.7 / 0.5,
#                            time span of flagged points, refit without them
#   dadras_loo_flagged.csv   one row per observation with k > 0.5: time, k,
#                            state, leverage, residual, derivative error
#   dadras_loo_khat_series.csv  block-maximum k-hat (blocks of 100 obs) for
#                            every fit, plus the full k-hat vector for the
#                            first FULL_SERIES trials at each n
# Usage: Rscript dadras_psis_loo.R [trials="1,2,...,20"] [log10n="4,4.4,4.6,4.8,5"]
#   Default: the first 20 initial conditions of the benchmark pool (the same
#   convention as exp2 for Point 6) at five points of the stored n grid.  Trial
#   j uses the j-th initial condition of the pool; the noise realisation is
#   re-drawn (the benchmark did not seed it), so the outcome is classified from
#   this run, not from the stored one.

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
suppressPackageStartupMessages(library(loo))
args <- commandArgs(trailingOnly = TRUE)
TRIALS <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:20
LOG10N <- if (length(args) >= 2) as.numeric(strsplit(args[2], ",")[[1]]) else c(4, 4.4, 4.6, 4.8, 5)
OUT <- file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "results")
dir.create(OUT, showWarnings = FALSE)

SYSTEM <- "dadras"; EQ <- 1; SNR <- 49
K_FLAG <- 0.7; K_RECORD <- 0.5
BLOCK <- 100L          # block length for the block-maximum k-hat series
FULL_SERIES <- 1L      # keep the full k-hat vector for this many trials per n (trial 1, used by the figure)
CHUNK <- 4000L         # observations per log-likelihood chunk (4000 draws x 4000 obs x 8 B = 128 MB)

# The benchmark's pool: set.seed(100); 100 draws of runif(3, -4, 4), one per trial.
dadras_ic_pool <- function(num = 100, seed = 100) {
  set.seed(seed)
  t(sapply(seq_len(num), function(i) runif(3, min = -4, max = 4)))
}

# Pointwise Pareto k-hat from an rstanarm Gaussian fit, chunked over observations.
pareto_k_chunked <- function(fit, y, X, intercept) {
  draws <- as.matrix(fit)                       # S x (p [+1] + 1)
  cn <- gsub("`", "", colnames(draws))
  sigma <- draws[, cn == "sigma"]
  beta_names <- c(if (intercept) "(Intercept)", colnames(X))
  beta <- draws[, match(beta_names, cn), drop = FALSE]
  Xd <- if (intercept) cbind(1, as.matrix(X)) else as.matrix(X)
  S <- nrow(draws)
  n_chain <- length(fit$stanfit@sim$samples)
  chain_id <- rep(seq_len(n_chain), each = S / n_chain)
  k <- numeric(length(y))
  for (start in seq(1, length(y), by = CHUNK)) {
    idx <- start:min(start + CHUNK - 1, length(y))
    mu <- beta %*% t(Xd[idx, , drop = FALSE])   # S x n_chunk
    ll <- dnorm(matrix(y[idx], S, length(idx), byrow = TRUE), mu, sigma, log = TRUE)
    r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id, cores = 1)
    ps <- suppressWarnings(loo::psis(-ll, r_eff = r_eff, cores = 1))
    k[idx] <- loo::pareto_k_values(ps)
  }
  k
}

leverage <- function(X, intercept) {
  Xd <- if (intercept) cbind(1, as.matrix(X)) else as.matrix(X)
  q <- qr(Xd)
  Q <- qr.Q(q)
  rowSums(Q^2)
}

cluster_count <- function(idx, gap = 1000L) {
  if (length(idx) == 0) return(0L)
  idx <- sort(idx)
  1L + sum(diff(idx) > gap)
}

one <- function(trial, log10n, keep_full) {
  n <- round(10^log10n)
  tr <- SYSTEMS[[SYSTEM]]$truth[[as.character(EQ)]]
  ic <- dadras_ic_pool()[trial, ]
  s <- simulate(SYSTEM, n, SNR, ic, 70000 + 1000 * trial + round(10 * log10n))
  dm <- library_dm(s$noisy)
  nm <- colnames(dm$sorted_theta)
  y <- dm$xdot_filtered[, EQ]
  tt <- (seq_along(y) - 1) * DT
  f_clean <- rhs_matrix(SYSTEM, s$clean)[, EQ]
  d <- y - f_clean                                   # derivative-estimation error

  # --- pipeline: screening, then the Bayesian stage on the selected support
  scr <- screen(dm, EQ)
  sel <- scr$selected_terms
  if (length(sel) == 0) sel <- tr$terms              # never happens for Dadras x1-dot; guard only
  Xs <- theta_columns(dm$sorted_theta, nm, sel); colnames(Xs) <- sel
  t0 <- Sys.time()
  bs <- bayes_stage(y, Xs, scr$selected_intercept)
  t_stan <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time()
  k <- pareto_k_chunked(bs$fit, y, Xs, scr$selected_intercept)
  t_loo <- as.numeric(Sys.time() - t0, units = "secs")
  h <- leverage(Xs, scr$selected_intercept)
  e <- bs$residuals
  flagged <- which(k > K_FLAG)
  recorded <- which(k > K_RECORD)

  # --- the causal claim: refit on the same support without the flagged points
  refit_int_kept <- NA; refit_int_mean <- NA
  if (length(flagged) > 0 && scr$selected_intercept) {
    bs2 <- bayes_stage(y[-flagged], Xs[-flagged, , drop = FALSE], TRUE)
    refit_int_kept <- bs2$kept_intercept
    refit_int_mean <- unname(bs2$coef["(Intercept)"])
  }

  # --- the competing explanation: intercept bias vs its shrinking uncertainty
  Xt <- theta_columns(dm$sorted_theta, nm, tr$terms); colnames(Xt) <- tr$terms
  ols <- lm(y ~ Xt)
  ols_int <- summary(ols)$coefficients["(Intercept)", ]
  # bias decomposition of the true-support residual, e = d + g + h (as in exp2):
  #   d = xdot_SG - f(x_clean)   derivative-estimation error
  #   g = f(x_clean) - Theta(x_SG) b*   state-smoothing error (library built on smoothed states)
  g <- f_clean - as.numeric(Xt %*% tr$beta)
  # the intercept the OLS fit would need if the library were evaluated on the clean states
  Xc <- poly_features(s$clean)$theta[, tr$terms, drop = FALSE]
  ols_clean_int <- summary(lm(y ~ Xc))$coefficients["(Intercept)", ]

  int_mean <- if (scr$selected_intercept) unname(bs$coef["(Intercept)"]) else NA
  int_lo <- if (scr$selected_intercept) unname(bs$ci[1, "(Intercept)"]) else NA
  int_hi <- if (scr$selected_intercept) unname(bs$ci[2, "(Intercept)"]) else NA

  row <- data.frame(
    trial = trial, log10n = log10n, n = n, snr = SNR,
    selected = paste(to_py_name(sel), collapse = ";"), selected_intercept = scr$selected_intercept,
    final = paste(to_py_name(bs$kept_terms), collapse = ";"), final_intercept = bs$kept_intercept,
    success = setequal(bs$kept_terms, tr$terms) && !bs$kept_intercept,
    int_mean = int_mean, int_lo = int_lo, int_hi = int_hi,
    ols_int_true_support = unname(ols_int["Estimate"]), ols_int_se = unname(ols_int["Std. Error"]),
    mean_d = mean(d), sd_d = sd(d), mean_g = mean(g), sd_g = sd(g),
    ols_int_clean_library = unname(ols_clean_int["Estimate"]), ols_int_clean_library_se = unname(ols_clean_int["Std. Error"]),
    sigma = bs$sigma,
    n_k_gt_0.7 = length(flagged), n_k_gt_0.5 = length(recorded), k_max = max(k),
    frac_k_gt_0.7 = length(flagged) / n,
    flagged_t_min = if (length(flagged)) tt[min(flagged)] else NA,
    flagged_t_max = if (length(flagged)) tt[max(flagged)] else NA,
    flagged_clusters = cluster_count(flagged), recorded_clusters = cluster_count(recorded),
    refit_without_flagged_int_kept = refit_int_kept, refit_without_flagged_int_mean = refit_int_mean,
    t_stan_s = t_stan, t_loo_s = t_loo
  )

  q_x2x3 <- ecdf(abs(s$clean[, 2] * s$clean[, 3]))
  q_r <- ecdf(sqrt(rowSums(s$clean^2)))
  q_h <- ecdf(h)
  q_absd <- ecdf(abs(d))
  fl <- if (length(recorded)) data.frame(
    trial = trial, log10n = log10n, n = n, index = recorded - 1, t = tt[recorded], k = k[recorded],
    x1 = s$clean[recorded, 1], x2 = s$clean[recorded, 2], x3 = s$clean[recorded, 3],
    abs_x2x3 = abs(s$clean[recorded, 2] * s$clean[recorded, 3]), abs_x2x3_quantile = q_x2x3(abs(s$clean[recorded, 2] * s$clean[recorded, 3])),
    radius_quantile = q_r(sqrt(rowSums(s$clean[recorded, , drop = FALSE]^2))),
    leverage = h[recorded], leverage_quantile = q_h(h[recorded]), leverage_x_n_over_p = h[recorded] * n / ncol(Xs),
    residual = e[recorded], std_residual = e[recorded] / bs$sigma,
    deriv_error = d[recorded], deriv_error_quantile = q_absd(abs(d[recorded])),
    y = y[recorded], fitted = bs$fitted[recorded]
  ) else NULL

  blk <- ceiling(seq_along(k) / BLOCK)
  series <- data.frame(trial = trial, log10n = log10n, kind = "block_max",
                       index = (tapply(seq_along(k), blk, min) - 1), t = (tapply(seq_along(k), blk, min) - 1) * DT,
                       k = tapply(k, blk, max))
  if (keep_full && log10n %in% c(4, 5)) series <- rbind(series, data.frame(trial = trial, log10n = log10n, kind = "full",   # the two panels of the figure
                                                    index = seq_along(k) - 1, t = tt, k = k))
  list(row = row, flagged = fl, series = series)
}

tasks <- list()
for (log10n in LOG10N) for (i in seq_along(TRIALS))
  tasks[[length(tasks) + 1]] <- list(TRIALS[i], log10n, i <= FULL_SERIES)
# largest n first so the long fits start immediately
tasks <- tasks[order(-sapply(tasks, `[[`, 2))]
cat(length(tasks), "fits on", NCORES, "cores\n")
t0 <- Sys.time()
res <- mclapply(tasks, function(t) try(one(t[[1]], t[[2]], t[[3]])), mc.cores = NCORES)
ok <- !sapply(res, inherits, "try-error")
if (any(!ok)) for (r in res[!ok]) cat("FAILED:", conditionMessage(attr(r, "condition")), "\n")
rows <- do.call(rbind, lapply(res[ok], `[[`, "row"))
write.csv(rows[order(rows$log10n, rows$trial), ], file.path(OUT, "dadras_loo_trials.csv"), row.names = FALSE)
fl <- do.call(rbind, Filter(Negate(is.null), lapply(res[ok], `[[`, "flagged")))
write.csv(fl, file.path(OUT, "dadras_loo_flagged.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(res[ok], `[[`, "series")), file.path(OUT, "dadras_loo_khat_series.csv"), row.names = FALSE)
cat("done in", round(as.numeric(Sys.time() - t0, units = "secs")), "s;", sum(!ok), "failures\n")
print(rows[order(rows$log10n, rows$trial), c("trial", "log10n", "final_intercept", "int_mean", "int_lo", "int_hi",
                                              "ols_int_true_support", "ols_int_se", "mean_d", "mean_g",
                                              "n_k_gt_0.7", "flagged_clusters", "refit_without_flagged_int_kept", "t_stan_s")])
