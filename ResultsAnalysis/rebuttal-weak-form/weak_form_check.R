# Weak-form (integral) check for the derivative-error failures of Referee 1
# Points 6 and 7, with the manuscript's own pipeline.
#
# Both replies trace the remaining failures to the Savitzky-Golay (SG)
# derivative estimate: its deterministic truncation error becomes
# "significant" under the credible-interval rule once the residual variance
# vanishes (Aizawa / Rossler at SNR = infinity, Point 6) or once n is large
# enough to shrink the interval below a fixed +0.005 intercept bias (Dadras
# x1-dot at n >= 6e4, Point 7).  A referee may ask why the pipeline
# differentiates at all when weak-form / integral formulations (Schaeffer &
# McCalla 2017; Reinbold et al. 2020; Messenger & Bortz 2021) avoid pointwise
# differentiation.  This script answers by substitution: the SG step of
# build_design_matrix is replaced by a weak-form design and everything
# downstream (two adaptive-lasso passes, stan_glm, 90 % posterior-interval
# rule) is the pipeline's own code, unchanged.
#
# Weak form.  For xdot = Theta(x) beta and a test function phi_k compactly
# supported on [t_k - h, t_k + h] with phi_k = 0 at the ends,
#     int phi_k xdot dt = - int phi_k' x dt = int phi_k Theta(x) dt  beta,
# so row k of the regression is  b_k = -int phi_k' x dt  against
# G_k = int phi_k Theta(x) dt.  phi(t) = (1 - ((t - t_k)/h)^2)^p on a support
# of 2m + 1 samples (h = m dt), normalised so that int phi dt = 1 (the constant
# column of the library stays 1 and the coefficients keep their meaning);
# trapezoid quadrature on the samples; query points every m samples (50 %
# overlap).  Reference: Messenger & Bortz, Multiscale Model. Simul. 19 (2021).
#
# Arms per trial (all share the same data and the same SG parameters):
#   sg              the pipeline as submitted (SG smoothing + SG derivative)
#   weak/raw        weak form on the raw noisy states (the WSINDy design)
#   weak/smoothed   weak form on the SG-smoothed states (isolates the
#                   derivative step: same smoothing, no differentiation)
# each weak arm at half-width m = w and 2w (w = the SG window length the
# pipeline's own grid search selects for the response coordinate) and test-
# function degree p = 4 and 8.
#
# Cases:
#   dadras  x1-dot, 49 dB, n = 1e4, 1e5, 10^5.5 (Point 7; the last n is beyond
#           the benchmark grid and shows where the SG failure rate is heading)
#   aizawa  x3-dot, n = 5000, SNR 60 dB and infinity (Point 6)
#   rossler x3-dot, n = 5000, SNR 61 dB and infinity (Point 6, secondary)
# Initial conditions and noise seeds are those of dadras_psis_loo.R and
# exp2_heteroscedasticity.R, so the sg arm reproduces those runs.
#
# Per (case, trial, arm) the output records the screened and kept supports,
# exact recovery, the OLS intercept on the true support with its s.e. (the
# bias of Point 7), the residual sd and lag-1 autocorrelation on the true
# support, the deterministic error of the response on the *clean* signal
# (SG truncation error / weak-form quadrature error), and the coefficient
# error of the kept model.
#
# Usage: Rscript weak_form_check.R [trials=20] [cases="dadras,aizawa,rossler"]
# Outputs: results/weak_form_trials.csv, results/weak_form_settings.csv

args <- commandArgs(trailingOnly = TRUE)
TRIALS <- if (length(args) >= 1) as.integer(args[1]) else 20L
CASES_RUN <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("dadras", "aizawa", "rossler")
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
OUT <- file.path(HERE, "results"); dir.create(OUT, showWarnings = FALSE)
setwd(file.path(HERE, "..", "rebuttal-diagnostics"))
source("common.R")
setwd(HERE)
# common.R carries beta only for the equations its own experiments fit; add
# Rossler x3-dot (0.2 + x1 x3 - 5.7 x3) locally, without touching common.R.
SYSTEMS$rossler$truth[["3"]]$beta <- c(-5.7, 1); SYSTEMS$rossler$truth[["3"]]$beta0 <- 0.2

CASES <- list(
  dadras  = list(eq = 1, snrs = c(49), log10n = c(4, 5, 5.5)),
  aizawa  = list(eq = 3, snrs = c(60, Inf), log10n = c(log10(5000))),
  rossler = list(eq = 3, snrs = c(61, Inf), log10n = c(log10(5000)))
)
WEAK_GRID <- expand.grid(states = c("raw", "smoothed"), m_mult = c(1, 2), p = c(4, 8),
                         stringsAsFactors = FALSE)
STRIDE_FRAC <- 0.5     # query points every m samples (50 % overlap)
CI_LEVEL <- 0.90

# The benchmark's Dadras pool: set.seed(100); 100 draws of runif(3, -4, 4).
dadras_ic_pool <- function(num = 100, seed = 100) {
  set.seed(seed)
  t(sapply(seq_len(num), function(i) runif(3, min = -4, max = 4)))
}
case_ic_seed <- function(system, trial, snr, log10n) {
  if (system == "dadras") {
    list(ic = dadras_ic_pool()[trial, ], seed = 70000 + 1000 * trial + round(10 * log10n))   # dadras_psis_loo.R
  } else {
    t0 <- trial - 1                                                                            # exp2: trials 0..19
    list(ic = initial_conditions(system)[t0 + 1, ],
         seed = 20000 + 100 * t0 + if (is.infinite(snr)) 99 else snr)
  }
}

# ---------------------------------------------------------------- weak form
# Test function and its derivative on the 2m+1 support samples, with the
# trapezoid weights folded in and normalised so that sum(phi * w) = 1.
test_function <- function(m, p, dt = DT) {
  h <- m * dt; u <- (-m:m) * dt / h
  phi <- (1 - u^2)^p
  dphi <- -2 * p * u * (1 - u^2)^(p - 1) / h
  w <- rep(dt, 2 * m + 1); w[c(1, 2 * m + 1)] <- dt / 2
  c0 <- sum(phi * w)
  list(phi_w = phi * w / c0, dphi_w = dphi * w / c0)
}

# Weak-form design in the layout build_design_matrix returns: sorted_theta
# (columns in the pipeline's order `nm`), monomial_orders, xdot_filtered
# (here b_k = -int phi_k' x dt for every coordinate).  `x_lib` are the states
# the library is evaluated on, `x_resp` the states the response integrates.
weak_design <- function(x_lib, x_resp, m, p, stride, nm) {
  n <- nrow(x_lib)
  tf <- test_function(m, p)
  centres <- seq(m + 1, n - m, by = stride)
  pf <- poly_features(x_lib)
  ord <- match(nm, pf$names)
  stopifnot(!anyNA(ord), length(nm) == length(pf$names))
  theta <- pf$theta[, ord, drop = FALSE]
  G <- matrix(0, length(centres), ncol(theta)); B <- matrix(0, length(centres), ncol(x_resp))
  for (j in seq_len(2 * m + 1)) {
    idx <- centres + j - m - 1
    G <- G + tf$phi_w[j] * theta[idx, , drop = FALSE]
    B <- B - tf$dphi_w[j] * x_resp[idx, , drop = FALSE]
  }
  colnames(G) <- nm
  list(sorted_theta = as.data.frame(G), monomial_orders = pf$degree[ord], xdot_filtered = B,
       centres = centres, K = length(centres))
}

sg_smooth <- function(x_col, order, window) sgolayfilt(x_col, p = order, n = window, m = 0, ts = DT)

rho1 <- function(e) cor(e[-length(e)], e[-1])

# One arm: screening + Stan on the selected support, OLS on the true support,
# deterministic response error on the clean signal.
run_arm <- function(arm, dm, eq, tr, y_clean_resp, theta_clean_truth, snr) {
  nm <- colnames(dm$sorted_theta)
  y <- dm$xdot_filtered[, eq]
  t0 <- Sys.time()
  scr <- screen(dm, eq)
  sel <- scr$selected_terms
  kept <- character(0); kept_int <- FALSE; int_mean <- NA; int_lo <- NA; int_hi <- NA; coef_err <- NA; sigma_kept <- NA
  if (length(sel) >= 1 || scr$selected_intercept) {
    Xs <- theta_columns(dm$sorted_theta, nm, sel); colnames(Xs) <- sel
    if (length(sel) == 0) Xs <- matrix(numeric(0), nrow = length(y), ncol = 0)
    bs <- bayes_stage(y, Xs, scr$selected_intercept, ci_level = CI_LEVEL)
    kept <- bs$kept_terms; kept_int <- bs$kept_intercept; sigma_kept <- bs$sigma
    if (scr$selected_intercept) {
      int_mean <- unname(bs$coef["(Intercept)"]); int_lo <- unname(bs$ci[1, "(Intercept)"]); int_hi <- unname(bs$ci[2, "(Intercept)"])
    }
    # coefficient error of the kept model on the true terms (posterior means; missing term = 0)
    cf <- bs$coef; names(cf) <- gsub("`", "", names(cf))
    est <- sapply(tr$terms, function(t) if (t %in% names(cf) && t %in% kept) unname(cf[t]) else 0)
    coef_err <- max(abs(est - tr$beta) / abs(tr$beta))
  }
  t_stan <- as.numeric(Sys.time() - t0, units = "secs")
  # OLS on the true support: intercept bias and residual statistics
  Xt <- theta_columns(dm$sorted_theta, nm, tr$terms); colnames(Xt) <- tr$terms
  ols <- lm(y ~ Xt)
  s <- summary(ols)$coefficients
  e <- residuals(ols)
  # deterministic error of the response on the clean signal: y_clean_resp - Theta(x_clean) beta*
  d_det <- y_clean_resp - as.numeric(theta_clean_truth %*% tr$beta) - if (tr$intercept) tr$beta0 else 0
  success <- setequal(kept, tr$terms) && (kept_int == tr$intercept)
  data.frame(
    arm = arm, K = length(y),
    selected = paste(to_py_name(sel), collapse = ";"), selected_intercept = scr$selected_intercept, n_selected = length(sel),
    kept = paste(to_py_name(kept), collapse = ";"), kept_intercept = kept_int, n_kept = length(kept),
    n_true_kept = length(intersect(kept, tr$terms)), n_spurious_selected = length(setdiff(sel, tr$terms)),
    n_spurious_kept = length(setdiff(kept, tr$terms)), success = success, max_rel_coef_err = coef_err,
    int_mean = int_mean, int_lo = int_lo, int_hi = int_hi,
    ols_int_true_support = unname(s["(Intercept)", "Estimate"]), ols_int_se = unname(s["(Intercept)", "Std. Error"]),
    ols_int_t = unname(s["(Intercept)", "t value"]),
    sigma_true_support = sd(e), rho1_true_support = rho1(e), sigma_kept = sigma_kept,
    det_error_mean = mean(d_det), det_error_sd = sd(d_det), det_error_max = max(abs(d_det)),
    t_stan_s = t_stan, stringsAsFactors = FALSE
  )
}

one <- function(system, trial, snr, log10n) {
  eq <- CASES[[system]]$eq
  tr <- SYSTEMS[[system]]$truth[[as.character(eq)]]
  n <- round(10^log10n)
  cs <- case_ic_seed(system, trial, snr, log10n)
  s <- simulate(system, n, snr, cs$ic, cs$seed)
  dm <- library_dm(s$noisy)
  nm <- colnames(dm$sorted_theta)
  sgp <- t(sapply(1:3, function(i) sg_params(s$noisy[, i])))     # (order, window) per coordinate
  w <- as.integer(sgp[eq, "window"])
  x_smoothed <- sapply(1:3, function(i) sg_smooth(s$noisy[, i], sgp[i, "order"], sgp[i, "window"]))
  theta_clean <- poly_features(s$clean)$theta
  f_clean <- rhs_matrix(system, s$clean)[, eq]
  base <- data.frame(system = system, eq = eq, trial = trial, snr = snr, log10n = log10n, n = n,
                     sg_order = sgp[eq, "order"], sg_window = w, stringsAsFactors = FALSE)
  rows <- list()
  # --- sg arm: response error on the clean signal = SG truncation error
  y_clean_sg <- sg_derivative(s$clean[, eq], sgp[eq, "order"], sgp[eq, "window"])
  theta_clean_sg <- sapply(1:3, function(i) sg_smooth(s$clean[, i], sgp[i, "order"], sgp[i, "window"]))
  theta_clean_sg <- poly_features(theta_clean_sg)$theta[, tr$terms, drop = FALSE]
  r <- run_arm("sg", dm, eq, tr, y_clean_sg, theta_clean_sg, snr)
  rows[[1]] <- cbind(base, data.frame(states = "smoothed", m = NA, p = NA, stride = NA), r)
  # --- weak arms
  for (g in seq_len(nrow(WEAK_GRID))) {
    m <- as.integer(WEAK_GRID$m_mult[g] * w); p <- WEAK_GRID$p[g]; stride <- max(1L, as.integer(round(STRIDE_FRAC * m)))
    x_use <- if (WEAK_GRID$states[g] == "raw") s$noisy else x_smoothed
    dmw <- weak_design(x_use, x_use, m, p, stride, nm)
    # deterministic response error on the clean signal (quadrature error): weak form of the clean states
    dmc <- weak_design(s$clean, s$clean, m, p, stride, nm)
    th_c <- theta_columns(dmc$sorted_theta, nm, tr$terms)
    arm <- sprintf("weak/%s/m=%dw/p=%d", WEAK_GRID$states[g], WEAK_GRID$m_mult[g], p)
    r <- run_arm(arm, dmw, eq, tr, dmc$xdot_filtered[, eq], th_c, snr)
    rows[[length(rows) + 1]] <- cbind(base, data.frame(states = WEAK_GRID$states[g], m = m, p = p, stride = stride), r)
  }
  do.call(rbind, rows)
}

tasks <- list()
for (system in CASES_RUN) for (snr in CASES[[system]]$snrs) for (log10n in CASES[[system]]$log10n)
  for (trial in seq_len(TRIALS)) tasks[[length(tasks) + 1]] <- list(system, trial, snr, log10n)
tasks <- tasks[order(-sapply(tasks, `[[`, 4))]          # largest n first
cat(length(tasks), "tasks x", 1 + nrow(WEAK_GRID), "arms on", NCORES, "cores\n")
t0 <- Sys.time()
res <- mclapply(tasks, function(t) try(one(t[[1]], t[[2]], t[[3]], t[[4]])), mc.cores = NCORES)
ok <- !sapply(res, inherits, "try-error")
if (any(!ok)) for (i in which(!ok)) cat("FAILED:", paste(unlist(tasks[[i]]), collapse = " "), conditionMessage(attr(res[[i]], "condition")), "\n")
rows <- do.call(rbind, res[ok])
rows <- rows[order(rows$system, rows$snr, rows$log10n, rows$trial, rows$arm), ]
# one file per system so that a case can be re-run alone; the combined file and
# the summary are rebuilt from every per-system file present
for (sys in unique(rows$system)) write.csv(rows[rows$system == sys, ], file.path(OUT, sprintf("weak_form_trials_%s.csv", sys)), row.names = FALSE)
rows <- do.call(rbind, lapply(list.files(OUT, "^weak_form_trials_.*\\.csv$", full.names = TRUE), read.csv))
write.csv(rows, file.path(OUT, "weak_form_trials.csv"), row.names = FALSE)
write.csv(cbind(WEAK_GRID, stride_frac = STRIDE_FRAC, ci_level = CI_LEVEL), file.path(OUT, "weak_form_settings.csv"), row.names = FALSE)
cat("done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min;", sum(!ok), "failures\n")

# ---------------------------------------------------------------- summary
summ <- rows %>% group_by(system, eq, snr, log10n, arm) %>% summarise(
  trials = n(), success = mean(success), intercept_kept = mean(kept_intercept),
  n_selected = median(n_selected), n_kept = median(n_kept), spurious_kept = median(n_spurious_kept),
  ols_int_med = median(ols_int_true_support), ols_int_se_med = median(ols_int_se), ols_int_t_med = median(ols_int_t),
  sigma_med = median(sigma_true_support), rho1_med = median(rho1_true_support),
  det_sd_med = median(det_error_sd), det_max_med = median(det_error_max), K_med = median(K),
  t_stan_med = median(t_stan_s), .groups = "drop")
write.csv(summ, file.path(OUT, "weak_form_summary.csv"), row.names = FALSE)
print(as.data.frame(summ[, c("system", "snr", "log10n", "arm", "success", "intercept_kept", "n_kept", "spurious_kept", "ols_int_med", "ols_int_se_med", "sigma_med", "det_sd_med")]), digits = 3)
