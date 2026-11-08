# How the Aizawa x3-dot fits fail under the Savitzky-Golay derivative (stored
# benchmark, Experiments/bayesian-alasso-ro/aizawa/results) and under the
# weak-form design (Experiments/bayesian-alasso-ro-weak/aizawa/results), at
# n = 10^4 (49 dB, the lower-range file's copy of the shared end point, as the
# benchmark table uses) and at 61 dB (n = 5000), 100 trials each.
# A term counts as identified when its coefficient is nonzero (the stored
# identified_model keeps zero rows; MethodsEvaluation's model_correctness_checker
# drops them before matching).  Failures are classed as
#   substitution   at least one true term missing and at least one spurious term,
#   missing only   true terms missing, no spurious term,
#   spurious only  every true term kept, spurious terms added.
# Written 2026-11-08 for Supplementary Note 4 (the loss of the weak form at
# finite SNR: added failures are mostly spurious-only, the pattern of credible
# intervals that are too narrow, not chain substitutions).
# Output: results/aizawa_failure_patterns_by_arm.csv (one row per arm x
# condition: exact, substitution, missing_only, spurious_only, median terms)
# and the most frequent missing and spurious terms on stdout.
# Usage: Rscript aizawa_failure_patterns_by_arm.R      (about a minute)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
root <- normalizePath(file.path(here, "..", ".."))
truth <- c("(Intercept)", "x3", "x1^2", "x2^2", "x1^2x3", "x2^2x3", "x1^3x3", "x3^3")
arms <- c(sg = "Experiments/bayesian-alasso-ro/aizawa/results", weak = "Experiments/bayesian-alasso-ro-weak/aizawa/results")
conds <- list(sg = list(c("aizawa_3_e3.5_e4.0_se100_snr49.RData", "n=1e+4"), c("aizawa_3_snr1_snr61_se100_n5000.RData", "snr=61")),
              weak = list(c("aizawa_3_e3.5_e4_se100_snr49.RData", "n=1e+4"), c("aizawa_3_snr1_snr61_se100_n5000.RData", "snr=61")))
rows <- list()
for (arm in names(arms)) for (cf in conds[[arm]]) {
  e <- new.env(); nm <- load(file.path(root, arms[[arm]], cf[1]), envir = e); x <- get(nm[1], envir = e)[[cf[2]]]
  terms <- lapply(x, function(t) { m <- t$id_result$argos_bi_output$identified_model; rownames(m)[m[, 1] != 0] })
  miss <- sapply(terms, function(tt) length(setdiff(truth, tt))); spur <- sapply(terms, function(tt) length(setdiff(tt, truth)))
  ok <- miss == 0 & spur == 0
  rows[[length(rows) + 1]] <- data.frame(arm = arm, condition = cf[2], n_trials = length(terms), exact = sum(ok),
    substitution = sum(miss > 0 & spur > 0), missing_only = sum(miss > 0 & spur == 0), spurious_only = sum(miss == 0 & spur > 0),
    median_terms = median(sapply(terms, length)))
  cat(sprintf("\n%-5s %-7s exact %3d  substitution %2d  missing only %2d  spurious only %2d  (median %g terms)\n", arm, cf[2], sum(ok), sum(miss > 0 & spur > 0), sum(miss > 0 & spur == 0), sum(miss == 0 & spur > 0), median(sapply(terms, length))))
  ms <- sort(table(unlist(lapply(terms[!ok], function(tt) setdiff(truth, tt)))), decreasing = TRUE); if (length(ms)) cat("  missing true terms:", paste(sprintf("%s %d", names(ms), ms), collapse = ", "), "\n")
  sp <- sort(table(unlist(lapply(terms[!ok], function(tt) setdiff(tt, truth)))), decreasing = TRUE); if (length(sp)) cat("  most frequent spurious:", paste(sprintf("%s %d", names(sp)[1:min(6, length(sp))], sp[1:min(6, length(sp))]), collapse = ", "), "\n")
}
out <- do.call(rbind, rows); write.csv(out, file.path(here, "results", "aizawa_failure_patterns_by_arm.csv"), row.names = FALSE); print(out, row.names = FALSE)
