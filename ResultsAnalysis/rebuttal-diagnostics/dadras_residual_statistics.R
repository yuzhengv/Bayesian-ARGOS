# Residual statistics of the Dadras x1-dot fit on the true support over 20
# trials (Supplementary Note 3, "Residual statistics").
#
# Replaces the three-trial spot check stored as results/dadras_residual_structure.csv
# and results/dadras_intercept_bias_source.csv (2026-09-12, interactive, OLS on the
# true support, no script), which the SI quoted as "in three trials" (Yuzheng,
# 2026-11-07: "Why only 3 trials here?").  Same conventions as exp2 for Aizawa
# and Rossler: the pipeline's Stan stage (bayes_stage) on the true support,
# residuals about the posterior mean, the statistics of residual_stats_from
# (Breusch-Pagan R^2, decile variance ratio, lag-1 autocorrelation), and the
# decomposition e = d + g + h with the truncation error d_trunc of the SG
# derivative computed from the noise-free signal.  The OLS fit on the same
# support is recorded alongside (ols_* columns) to show that the two residuals
# give the same statistics, as the stored spot check assumed.
#
# Trials: the first 20 initial conditions of the benchmark pool with the noise
# seeds of dadras_psis_loo.R, so the 49 dB trajectories are those of the
# PSIS-LOO run; n = 10^4 and 10^5; SNR = 49 dB (the n sweep's level) and Inf.
# Also recorded per trial: the 99.9 % quantile and the maximum of |d| (the
# derivative-estimation error over the record), the OLS intercept of y on the
# true terms and of d_trunc on the true terms evaluated at the clean states (the
# systematic component of the truncation error that a constant absorbs), and the
# measurement-noise sd on x1 implied by the SNR.
#
# Usage: Rscript dadras_residual_statistics.R [trials="1,2,...,20"] [log10n="4,5"]
# Outputs: results/dadras_residual_statistics.csv (one row per trial, n, SNR),
#          tables/dadras_residual_statistics.md (medians and ranges per n, SNR)

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
args <- commandArgs(trailingOnly = TRUE)
TRIALS <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:20
LOG10N <- if (length(args) >= 2) as.numeric(strsplit(args[2], ",")[[1]]) else c(4, 5)
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
OUT <- file.path(HERE, "results"); TAB <- file.path(HERE, "tables")
dir.create(OUT, showWarnings = FALSE); dir.create(TAB, showWarnings = FALSE)

SYSTEM <- "dadras"; EQ <- 1; SNRS <- c(49, Inf)
TRUTH <- SYSTEMS[[SYSTEM]]$truth[[as.character(EQ)]]
# the benchmark's pool as drawn in dadras_psis_loo.R (set.seed(100); runif(3, -4, 4) per trial)
ICS <- { set.seed(100); t(sapply(seq_len(100), function(i) runif(3, min = -4, max = 4))) }
stopifnot(isTRUE(all.equal(ICS, unname(as.matrix(initial_conditions(SYSTEM))), check.attributes = FALSE)))

one <- function(trial, log10n, snr) {
  t0 <- Sys.time()
  n <- round(10^log10n)
  s <- simulate(SYSTEM, n, snr, ICS[trial, ], 70000 + 1000 * trial + round(10 * log10n))
  dm <- library_dm(s$noisy)
  nm <- colnames(dm$sorted_theta)
  y <- dm$xdot_filtered[, EQ]
  f_clean <- rhs_matrix(SYSTEM, s$clean)[, EQ]
  Xt <- theta_columns(dm$sorted_theta, nm, TRUTH$terms); colnames(Xt) <- TRUTH$terms
  theta_bstar <- as.numeric(Xt %*% TRUTH$beta) + if (TRUTH$intercept) TRUTH$beta0 else 0
  sgp <- sg_params(s$noisy[, EQ])
  d <- y - f_clean
  d_trunc <- sg_derivative(s$clean[, EQ], sgp["order"], sgp["window"]) - f_clean
  d_noise <- d - d_trunc
  g <- f_clean - theta_bstar

  bs <- bayes_stage(y, Xt, TRUTH$intercept)
  st <- residual_stats_from(bs$residuals, bs$fitted)
  e <- st$residuals
  h <- theta_bstar - st$fitted
  ols <- if (TRUTH$intercept) lm(y ~ Xt) else lm(y ~ 0 + Xt)
  st_ols <- residual_stats_from(as.numeric(residuals(ols)), as.numeric(fitted(ols)))
  Xc <- poly_features(s$clean)$theta[, TRUTH$terms, drop = FALSE]

  row <- data.frame(system = SYSTEM, eq = EQ, trial = trial, log10n = log10n, n = n, snr = snr,
                    sg_order = unname(sgp["order"]), sg_window = unname(sgp["window"]),
                    noise_sd_state = if (is.infinite(snr)) 0 else 10^(-snr / 20) * sqrt(mean((s$clean[, EQ] - mean(s$clean[, EQ]))^2)))
  for (k in c("sigma", "bp_r2", "var_ratio", "rho1", "durbin_watson")) { row[[k]] <- st[[k]]; row[[paste0("ols_", k)]] <- st_ols[[k]] }
  row$truth_support_all_kept <- setequal(bs$kept_terms, TRUTH$terms) && bs$kept_intercept == TRUTH$intercept
  row$var_e <- var(e); row$var_d <- var(d); row$var_dtrunc <- var(d_trunc); row$var_dnoise <- var(d_noise); row$var_g <- var(g); row$var_h <- var(h)
  row$corr_e_d <- cor(e, d); row$corr_e_dtrunc <- cor(e, d_trunc); row$corr_d_dtrunc <- cor(d, d_trunc)
  row$rho1_d <- cor(d[-length(d)], d[-1]); row$rho1_dtrunc <- cor(d_trunc[-length(d_trunc)], d_trunc[-1])
  row$q999_absd <- unname(quantile(abs(d), 0.999)); row$max_absd <- max(abs(d))
  row$q999_absdtrunc <- unname(quantile(abs(d_trunc), 0.999)); row$max_absdtrunc <- max(abs(d_trunc))
  row$int_y_on_theta <- unname(coef(lm(y ~ Xt))["(Intercept)"])
  row$int_dtrunc_on_clean <- unname(coef(lm(d_trunc ~ Xc))["(Intercept)"])
  row$seconds <- as.numeric(Sys.time() - t0, units = "secs")
  row
}

tasks <- list()
for (trial in TRIALS) for (l in LOG10N) for (snr in SNRS) tasks[[length(tasks) + 1]] <- list(trial, l, snr)
tasks <- tasks[order(-sapply(tasks, function(t) t[[2]]))]
cat(length(tasks), "conditions on", NCORES, "cores\n")
t0 <- Sys.time()
res <- mclapply(tasks, function(t) one(t[[1]], t[[2]], t[[3]]), mc.cores = NCORES)
ok <- !sapply(res, inherits, "try-error")
if (any(!ok)) print(res[!ok])
stats <- do.call(rbind, res[ok])
stats <- stats[order(stats$log10n, stats$snr, stats$trial), ]
write.csv(stats, file.path(OUT, "dadras_residual_statistics.csv"), row.names = FALSE)

fmt <- function(x) formatC(signif(x, 3), digits = 3, format = "g")
cell <- function(v) sprintf("%s [%s, %s]", fmt(median(v)), fmt(min(v)), fmt(max(v)))
lines <- c(sprintf("# Dadras x1-dot, true support, Stan stage: %d trials, %s", length(unique(stats$trial)), format(Sys.Date())),
           "median [min, max] over the trials", "",
           "| n | SNR | decile variance ratio | rho1 | Breusch-Pagan R2 | residual sd | noise sd on x1 | Var(d_trunc)/Var(e) | corr(e, d_trunc) | q99.9 |d| | max |d| | int(y on theta) | int(d_trunc on clean) | OLS decile ratio | OLS rho1 | OLS sd |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (l in sort(unique(stats$log10n))) for (snr in SNRS) {
  r <- stats[stats$log10n == l & stats$snr == snr, ]
  lines <- c(lines, sprintf("| 1e%g | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |", l, if (is.infinite(snr)) "Inf" else snr,
                            cell(r$var_ratio), cell(r$rho1), cell(r$bp_r2), cell(r$sigma), cell(r$noise_sd_state),
                            cell(r$var_dtrunc / r$var_e), cell(r$corr_e_dtrunc), cell(r$q999_absd), cell(r$max_absd),
                            cell(r$int_y_on_theta), cell(r$int_dtrunc_on_clean), cell(r$ols_var_ratio), cell(r$ols_rho1), cell(r$ols_sigma)))
}
lines <- c(lines, "", sprintf("truth kept by the credible-interval rule in %d of %d fits", sum(stats$truth_support_all_kept), nrow(stats)),
           sprintf("runtime: %.0f s on %d cores (%d failures)", as.numeric(Sys.time() - t0, units = "secs"), NCORES, sum(!ok)))
writeLines(lines, file.path(TAB, "dadras_residual_statistics.md")); cat(lines, sep = "\n")
