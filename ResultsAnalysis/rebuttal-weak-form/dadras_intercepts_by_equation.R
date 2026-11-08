# Per-equation exact recovery and spurious-intercept counts for the Dadras system
# in the n sweep at 49 dB (n = 10^2 to 10^5, 100 trials per n), for the two
# Bayesian-ARGOS arms of the 2 x 2 comparison: the Savitzky-Golay derivative
# (stored benchmark, Experiments/bayesian-alasso-ro/dadras/results) and the
# weak-form design (Experiments/bayesian-alasso-ro-weak/dadras/results).
# Written 2026-11-08 for Supplementary Note 4 (Yuzheng: "Is removes
# appropriate here?"): the system-level success rate of the weak form at
# n = 10^5 is 88 %, and this table shows which equation the remaining
# failures come from and whether they are intercepts.
# A term counts as identified when its coefficient is nonzero (the stored
# identified_model keeps zero rows; MethodsEvaluation's model_correctness_checker
# drops them before matching, and so does this script).
# Output: results/dadras_intercepts_by_equation.csv (one row per arm x
# equation x n: exact, intercept, n_trials) and the table on stdout.
# Usage: Rscript dadras_intercepts_by_equation.R      (seconds)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
root <- normalizePath(file.path(here, "..", ".."))
truth <- list(c("x1", "x2", "x2x3"), c("x2", "x1x3", "x3"), c("x1x2", "x3"))
arms <- c(sg = "Experiments/bayesian-alasso-ro/dadras/results", weak = "Experiments/bayesian-alasso-ro-weak/dadras/results")
rows <- list()
for (arm in names(arms)) for (eq in 1:3) {
  nout <- list()
  for (rng in c("e2_e3", "e3_e4", "e4_e5")) {   # the three files of the n sweep, 10^2 to 10^5
    e <- new.env(); load(file.path(root, arms[[arm]], sprintf("dadras_%d_%s_se100_snr49.RData", eq, rng)), envir = e)
    for (nkey in names(e$n_output)) if (is.null(nout[[nkey]])) nout[[nkey]] <- e$n_output[[nkey]]   # the shared end point is stored once
  }
  for (nkey in names(nout)) {
    terms <- lapply(nout[[nkey]], function(t) { m <- t$id_result$argos_bi_output$identified_model; rownames(m)[m[, 1] != 0] })   # nonzero coefficients only, as model_correctness_checker() does
    rows[[length(rows) + 1]] <- data.frame(arm = arm, equation = eq, n = sub("n=", "", nkey),
      exact = sum(sapply(terms, function(tt) !is.null(tt) && setequal(tt, truth[[eq]]))),
      intercept = sum(sapply(terms, function(tt) !is.null(tt) && any(grepl("Intercept", tt)))),
      n_trials = length(terms))
  }
}
out <- do.call(rbind, rows)
dir.create(file.path(here, "results"), showWarnings = FALSE)
write.csv(out, file.path(here, "results", "dadras_intercepts_by_equation.csv"), row.names = FALSE)
for (arm in names(arms)) { cat("\n==", arm, "arm: exact / intercept out of 100\n")
  w <- out[out$arm == arm, ]; m <- t(sapply(split(w, w$equation), function(d) sprintf("%d/%d", d$exact, d$intercept)))
  colnames(m) <- unique(w$n); rownames(m) <- paste0("eq", rownames(m)); print(m, quote = FALSE) }
