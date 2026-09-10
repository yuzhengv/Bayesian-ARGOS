# Screening support vs final HMC model for the Aizawa x3-dot equation, read
# from the stored benchmark results (Experiments/bayesian-alasso-ro/aizawa).
# For every trial: size of the adaptive-lasso pass-2 support that enters HMC
# (`identified_model` rows), size of the final model after the 90% credible
# interval rule (`trimmed_model` rows), whether the final model is exactly the
# ground truth, and how many spurious pass-2 candidates survive the HMC rule.
# Outputs: aizawa_eq3_screening_vs_hmc.csv (one row per (sweep, value)) and
# aizawa_eq3_supports.csv (one row per trial: the pass-2 and final supports as
# "|"-separated term names, consumed by the Python exp4_screening_substitution.py).
root <- normalizePath(file.path(getwd(), "..", ".."))
res_dir <- file.path(root, "Experiments", "bayesian-alasso-ro", "aizawa", "results")
truth <- c("(Intercept)", "x3", "x1^2", "x2^2", "x1^2x3", "x2^2x3", "x3^3", "x1^3x3")

summarise_level <- function(trials, sweep, value) {
  rows <- lapply(trials, function(t) {
    out <- t$id_result$argos_bi_output
    if (is.null(out) || is.null(out$identified_model)) return(NULL)
    pass2 <- rownames(out$identified_model)
    final <- rownames(out$trimmed_model)
    spurious2 <- setdiff(pass2, truth)
    data.frame(p_pass2 = length(pass2), p_final = length(final),
               truth_in_pass2 = all(truth %in% pass2),
               exact = setequal(final, truth),
               n_spurious_pass2 = length(spurious2),
               n_spurious_final = length(setdiff(final, truth)),
               n_missing_final = length(setdiff(truth, final)),
               spurious_survive = if (length(spurious2)) mean(spurious2 %in% final) else NA)
  })
  d <- do.call(rbind, rows)
  data.frame(sweep = sweep, value = value, trials = nrow(d),
             p_pass2_median = median(d$p_pass2), p_final_median = median(d$p_final),
             truth_in_pass2 = mean(d$truth_in_pass2), success = mean(d$exact),
             n_spurious_pass2_mean = mean(d$n_spurious_pass2),
             spurious_survival = mean(d$spurious_survive, na.rm = TRUE),
             frac_any_spurious_final = mean(d$n_spurious_final > 0),
             frac_any_missing_final = mean(d$n_missing_final > 0))
}

supports <- list()
export_supports <- function(trials, sweep, value) {
  for (j in seq_along(trials)) {
    o <- trials[[j]]$id_result$argos_bi_output
    if (is.null(o) || is.null(o$identified_model)) next
    supports[[length(supports) + 1]] <<- data.frame(
      sweep = sweep, value = value, trial = j,
      pass2 = paste(rownames(o$identified_model), collapse = "|"),
      final = paste(rownames(o$trimmed_model), collapse = "|"))
  }
}

out <- list()
for (f in list.files(res_dir, pattern = "^aizawa_3_e.*_snr49\\.RData$", full.names = TRUE)) {
  load(f)
  for (nm in names(n_output)) {
    k <- round(as.numeric(sub("n=1e\\+", "", nm)), 1)  # names are "n=1e+2.5"
    if (any(sapply(out, function(o) o$sweep == "n" && o$value == k))) next  # chunk files overlap
    out[[length(out) + 1]] <- summarise_level(n_output[[nm]], "n", k)
    export_supports(n_output[[nm]], "n", k)
  }
}
load(file.path(res_dir, "aizawa_3_snr1_snr61_se100_n5000.RData"))
for (nm in names(snr_output)) {
  v <- sub("snr=", "", nm)
  vv <- if (v == "Inf") Inf else as.numeric(v)
  out[[length(out) + 1]] <- summarise_level(snr_output[[nm]], "snr", vv)
  export_supports(snr_output[[nm]], "snr", vv)
}
write.csv(do.call(rbind, supports), "aizawa_eq3_supports.csv", row.names = FALSE)
res <- do.call(rbind, out)
res <- res[order(res$sweep, res$value), ]
write.csv(res, "aizawa_eq3_screening_vs_hmc.csv", row.names = FALSE)
print(res[res$sweep == "n", ], digits = 3, row.names = FALSE)
print(res[res$sweep == "snr" & res$value %in% c(10, 20, 27, 30, 40, 45, 49, 50, 55, 57, 58, 59, 60, 61, Inf), ], digits = 3, row.names = FALSE)
