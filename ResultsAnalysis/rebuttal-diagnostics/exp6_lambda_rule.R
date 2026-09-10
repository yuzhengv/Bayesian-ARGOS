# Experiment 6 (R) -- is the large-n decline a consequence of tuning the
# adaptive lasso by the prediction-optimal lambda.min?
#
# The pipeline's `alasso` picks cv.glmnet's lambda.min in both passes.  Here
# the same function is patched in memory to use lambda.1se (the largest lambda
# within one standard error of the CV minimum) in pass 2 only, or in both
# passes, and the three variants are compared on identical Aizawa x3-dot data
# over the n grid (49 dB) and at n = 5000 for SNR 49, 60 and infinity.
#
# The Bayesian stage is the package's own (rstanarm::stan_glm on the selected
# support, 90 % posterior-interval keep rule as in bayesian_alasso_ro), so the
# reported success is the pipeline's success; the baseline variant should
# reproduce the stored 100-trial rates (0.76 at 10^3.7, 0.45 at 10^5).
#
# The package code is NOT modified.  Output: results/exp6_lambda_rule.csv
# Usage: Rscript exp6_lambda_rule.R [trials=20]

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
args <- commandArgs(trailingOnly = TRUE)
TRIALS <- if (length(args) >= 1) as.integer(args[1]) else 20L
OUT <- file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "results")
dir.create(OUT, showWarnings = FALSE)

EQ <- 3; SYSTEM <- "aizawa"
TR <- SYSTEMS[[SYSTEM]]$truth[[as.character(EQ)]]
TRUTH <- TR$terms
ICS <- initial_conditions(SYSTEM)
N_GRID <- log_grid(3.4, 5.0, 0.2)
SNR_FIXED <- 49

# --- in-memory variant of alasso that uses lambda.1se -----------------------
src <- deparse(alasso)
src <- gsub("alasso_init\\$lambda\\.min", "alasso_init$lambda.1se", src)
src <- gsub("alasso_new\\$lambda\\.min", "alasso_new$lambda.1se", src)
alasso_1se <- eval(parse(text = src))
environment(alasso_1se) <- environment(alasso)

screen_variant <- function(dm, eq, rule2, rule1 = "min") {
  f1 <- if (rule1 == "1se") alasso_1se else alasso
  f2 <- if (rule2 == "1se") alasso_1se else alasso
  theta <- dm$sorted_theta; target <- dm$xdot_filtered[, eq]
  data <- cbind.data.frame(target = target, theta)
  init <- f1(data, weights_method = "ridge", ols_ps = TRUE); init[is.na(init)] <- 0
  nz_max <- max(which(init != 0))
  cut <- sum(dm$monomial_orders <= dm$monomial_orders[nz_max])
  post <- if (is.na(cut) || cut == length(dm$monomial_orders)) data else cbind.data.frame(target = target, data[-1][, 1:cut])
  final <- f2(post, weights_method = "ols", ols_ps = TRUE); final[is.na(final)] <- 0
  nz <- final != 0
  list(selected_intercept = unname(nz[1]), selected_terms = colnames(post)[-1][nz[-1]], target = target)
}

# the pipeline's credible-interval rule (stan_glm, 90 % posterior interval)
ci_prune <- function(y, X, intercept, level = 0.90) {
  if (ncol(X) == 0) return(list(intercept = intercept, terms = character(0)))
  bs <- bayes_stage(y, X, intercept, ci_level = level)
  list(intercept = bs$kept_intercept, terms = bs$kept_terms)
}

one <- function(trial, k, n, snr) {
  s <- simulate(SYSTEM, n, snr, ICS[trial + 1, ], 40000 + 100 * trial + round(10 * k) + if (is.infinite(snr)) 99 else snr)
  dm <- library_dm(s$noisy); nm <- colnames(dm$sorted_theta)
  rows <- list()
  for (variant in c("min/min", "min/1se", "1se/1se")) {
    r <- strsplit(variant, "/")[[1]]
    scr <- screen_variant(dm, EQ, rule2 = r[2], rule1 = r[1])
    sel <- scr$selected_terms
    X <- theta_columns(dm$sorted_theta, nm, sel); colnames(X) <- sel
    fin <- ci_prune(scr$target, X, scr$selected_intercept)
    rows[[variant]] <- data.frame(
      trial = trial, log10n = k, n = n, snr = snr, variant = variant,
      p_selected = length(sel), truth_in_selected = all(TRUTH %in% sel),
      n_spurious_selected = length(setdiff(sel, TRUTH)),
      p_final = length(fin$terms) + fin$intercept,
      success = setequal(fin$terms, TRUTH) && fin$intercept == TR$intercept,
      truth_in_final = all(TRUTH %in% fin$terms) && fin$intercept == TR$intercept,
      n_spurious_final = length(setdiff(fin$terms, TRUTH)),
      selected = paste(to_py_name(sel), collapse = ";"))
  }
  do.call(rbind, rows)
}

tasks <- list()
for (trial in seq_len(TRIALS) - 1) {
  for (i in seq_len(nrow(N_GRID))) tasks[[length(tasks) + 1]] <- list(trial, N_GRID$log10n[i], N_GRID$n[i], SNR_FIXED)
  for (snr in c(60, Inf)) tasks[[length(tasks) + 1]] <- list(trial, log10(5000), 5000, snr)
}
tasks <- tasks[order(-sapply(tasks, function(t) t[[3]]))]
cat(length(tasks), "conditions x 3 variants on", NCORES, "cores\n")
t0 <- Sys.time()
res <- mclapply(tasks, function(t) one(t[[1]], t[[2]], t[[3]], t[[4]]), mc.cores = NCORES)
ok <- !sapply(res, inherits, "try-error")
out <- do.call(rbind, res[ok])
write.csv(out, file.path(OUT, "exp6_lambda_rule.csv"), row.names = FALSE)
cat("done in", round(as.numeric(Sys.time() - t0, units = "secs")), "s;", sum(!ok), "failures\n")

summ <- aggregate(cbind(p_selected, truth_in_selected, n_spurious_selected, p_final, success, truth_in_final, n_spurious_final) ~ variant + snr + log10n,
                  data = out, FUN = mean)
summ <- summ[order(summ$snr, summ$log10n, summ$variant), ]
print(summ, digits = 2, row.names = FALSE)
