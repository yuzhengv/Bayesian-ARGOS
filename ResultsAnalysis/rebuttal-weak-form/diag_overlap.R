# Diagnostic: does the overlap of the weak-form windows (K = K_FRAC * n
# subdomains of 50 samples) inflate the effective sample size seen by the
# Bayesian stage and produce spurious intercepts at finite SNR?
# Bayesian-ARGOS (Integration) on Dadras x3-dot and Aizawa x3-dot, n = 5000,
# 49 dB, 20 trials, K_FRAC in {0.25 (Hamilton run), 0.04 (50 % overlap),
# 0.02 (no overlap)}.  Usage: ARGOS_ROOT=<repo> ARGOS_PYTHON=<python> Rscript diag_overlap.R
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "..", "..", "Experiments", "bayesian-alasso-ro-weak", "weak_form_common.R"))
TRUTH <- list(dadras = list(int = FALSE, terms = c("x1x2", "x3")), aizawa = list(int = TRUE, terms = c("x3","x1^2","x2^2","x1^2x3","x2^2x3","x1^3x3","x3^3")))
KF <- c(0.25, 0.04, 0.02); N <- 5000; SNR <- 49; TRIALS <- 20
rows <- list()
for (system in c("dadras", "aizawa")) {
  init <- SYSTEMS[[system]]$init(100, 100)
  res <- mclapply(seq_len(TRIALS), function(j) {
    out <- list()
    for (kf in KF) { K_FRAC <<- kf
      r <- run_bayesian_argos_weak(system, N, init[[j]], 0.01, SNR, 5, "poly", 3, 0.9, 1, np_seed_for(100, 1, j))
      im <- r$argos_bi_output$identified_model[, 1]; sel <- names(im)[im != 0]; tr <- TRUTH[[system]]
      out[[length(out) + 1]] <- data.frame(system = system, trial = j, K_frac = kf, K = r$weak$K,
        exact = setequal(sel, c(if (tr$int) "(Intercept)", tr$terms)), intercept = "(Intercept)" %in% sel,
        n_spurious = length(setdiff(sel, c("(Intercept)", tr$terms))), n_missing = length(setdiff(tr$terms, sel))) }
    do.call(rbind, out) }, mc.cores = 8)
  rows[[system]] <- do.call(rbind, res)
}
d <- do.call(rbind, rows)
s <- aggregate(cbind(exact, intercept, n_spurious, n_missing) ~ system + K_frac + K, d, mean)
print(s[order(s$system, -s$K_frac), ], row.names = FALSE, digits = 2)
write.csv(d, file.path(here, "results", "diag_overlap.csv"), row.names = FALSE)
