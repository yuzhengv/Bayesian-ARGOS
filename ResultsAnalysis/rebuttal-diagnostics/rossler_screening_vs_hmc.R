# Screening support vs final HMC model for the Rossler x1-dot equation, read
# from the stored 100-trial benchmark results (Experiments/bayesian-alasso-ro/
# rossler), the counterpart of aizawa_screening_vs_hmc.R for the second
# mechanism of the reply (Referee 1 Points 6 and 7: the spurious intercept the
# credible-interval rule admits in the noiseless limit).  Per level of the SNR
# sweep (1 ... 61 dB and Inf, n = 5000) and of the n sweep (10^2 ... 10^5,
# 49 dB): median pass-2 support size, fraction of trials with the truth in the
# pass-2 support, success rate (final model exactly the truth), mean share of
# the spurious pass-2 candidates that survive the 90 % credible-interval rule,
# and the per-term selection frequency of the FINAL model (fraction of the
# 100 trials) for every term that appears.
#
# The same per-term final-model frequencies are exported for Dadras x1-dot
# (stored results, both sweeps) and for Aizawa x3-dot (from
# aizawa_eq3_supports.csv written by aizawa_screening_vs_hmc.R, column `final`).
# Term names are the row names of the stored models ("(Intercept)", "x1^2x3").
#
# Outputs (results/):
#   rossler_eq1_screening_vs_hmc.csv   one row per (sweep, value)
#   rossler_eq1_term_frequencies.csv   long: sweep, value, trials, term, frequency, is_truth
#   dadras_eq1_term_frequencies.csv    same schema
#   aizawa_eq3_term_frequencies.csv    same schema
# Usage: Rscript rossler_screening_vs_hmc.R        (runs in about 10 s)

here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
root <- normalizePath(file.path(here, "..", ".."))
OUT <- file.path(here, "results"); dir.create(OUT, showWarnings = FALSE)
exp_dir <- function(system) file.path(root, "Experiments", "bayesian-alasso-ro", system, "results")

TRUTH <- list(
  rossler = c("x2", "x3"),                          # x1-dot = -x2 - x3, no intercept
  dadras  = c("x1", "x2", "x2x3"),                  # x1-dot = x2 - 3 x1 + 2.7 x2 x3, no intercept
  aizawa  = c("(Intercept)", "x3", "x1^2", "x2^2", "x1^2x3", "x2^2x3", "x3^3", "x1^3x3")
)

# ------------------------------------------------------------- stored results
# Returns, per level, the list of (pass2, final) supports of the trials.
read_levels <- function(system, eq) {
  res_dir <- exp_dir(system)
  levels <- list()
  add <- function(trials, sweep, value) {
    key <- paste(sweep, value)
    if (!is.null(levels[[key]])) return(invisible())      # chunk files overlap at 10^3 and 10^4
    sup <- lapply(trials, function(t) {
      o <- t$id_result$argos_bi_output
      if (is.null(o) || is.null(o$identified_model)) return(NULL)
      list(pass2 = rownames(o$identified_model), final = rownames(o$trimmed_model))
    })
    levels[[key]] <<- list(sweep = sweep, value = value, supports = Filter(Negate(is.null), sup))
  }
  for (f in sort(list.files(res_dir, pattern = sprintf("^%s_%d_e.*_snr49\\.RData$", system, eq), full.names = TRUE))) {
    e <- new.env(); load(f, envir = e)
    for (nm in names(e$n_output)) add(e$n_output[[nm]], "n", round(as.numeric(sub("n=1e\\+", "", nm)), 1))
  }
  e <- new.env(); load(file.path(res_dir, sprintf("%s_%d_snr1_snr61_se100_n5000.RData", system, eq)), envir = e)
  for (nm in names(e$snr_output)) {
    v <- sub("snr=", "", nm)
    add(e$snr_output[[nm]], "snr", if (v == "Inf") Inf else as.numeric(v))
  }
  levels
}

summarise_level <- function(lv, truth) {
  d <- do.call(rbind, lapply(lv$supports, function(s) {
    spurious2 <- setdiff(s$pass2, truth)
    data.frame(p_pass2 = length(s$pass2), p_final = length(s$final),
               truth_in_pass2 = all(truth %in% s$pass2), exact = setequal(s$final, truth),
               n_spurious_pass2 = length(spurious2), n_spurious_final = length(setdiff(s$final, truth)),
               n_missing_final = length(setdiff(truth, s$final)),
               spurious_survive = if (length(spurious2)) mean(spurious2 %in% s$final) else NA)
  }))
  data.frame(sweep = lv$sweep, value = lv$value, trials = nrow(d),
             p_pass2_median = median(d$p_pass2), p_final_median = median(d$p_final),
             truth_in_pass2 = mean(d$truth_in_pass2), success = mean(d$exact),
             n_spurious_pass2_mean = mean(d$n_spurious_pass2),
             spurious_survival = mean(d$spurious_survive, na.rm = TRUE),
             frac_any_spurious_final = mean(d$n_spurious_final > 0),
             frac_any_missing_final = mean(d$n_missing_final > 0),
             frac_intercept_final = mean(sapply(lv$supports, function(s) "(Intercept)" %in% s$final)))
}

term_frequencies <- function(lv, truth) {
  finals <- lapply(lv$supports, `[[`, "final")
  terms <- sort(unique(unlist(finals)))
  if (length(terms) == 0) return(NULL)
  data.frame(sweep = lv$sweep, value = lv$value, trials = length(finals), term = terms,
             frequency = sapply(terms, function(t) mean(sapply(finals, function(f) t %in% f))),
             is_truth = terms %in% truth, row.names = NULL)
}

order_rows <- function(d) if ("term" %in% names(d)) d[order(d$sweep, d$value, d$term), ] else d[order(d$sweep, d$value), ]

# --------------------------------------------------------------- Rossler x1-dot
ros <- read_levels("rossler", 1)
res <- order_rows(do.call(rbind, lapply(ros, summarise_level, truth = TRUTH$rossler)))
write.csv(res, file.path(OUT, "rossler_eq1_screening_vs_hmc.csv"), row.names = FALSE)
fr <- order_rows(do.call(rbind, lapply(ros, term_frequencies, truth = TRUTH$rossler)))
write.csv(fr, file.path(OUT, "rossler_eq1_term_frequencies.csv"), row.names = FALSE)

# ---------------------------------------------------------------- Dadras x1-dot
dad <- read_levels("dadras", 1)
fd <- order_rows(do.call(rbind, lapply(dad, term_frequencies, truth = TRUTH$dadras)))
write.csv(fd, file.path(OUT, "dadras_eq1_term_frequencies.csv"), row.names = FALSE)

# ---------------------------------------------------------------- Aizawa x3-dot
sup <- read.csv(file.path(here, "aizawa_eq3_supports.csv"), stringsAsFactors = FALSE)
aiz <- lapply(split(sup, list(sup$sweep, sup$value), drop = TRUE), function(d)
  list(sweep = d$sweep[1], value = d$value[1], supports = lapply(strsplit(d$final, "|", fixed = TRUE), function(f) list(final = f))))
fa <- order_rows(do.call(rbind, lapply(aiz, term_frequencies, truth = TRUTH$aizawa)))
write.csv(fa, file.path(OUT, "aizawa_eq3_term_frequencies.csv"), row.names = FALSE)

# ----------------------------------------------------------------- report
cat("Rossler x1-dot, SNR sweep (n = 5000), 55 to 61 dB and Inf:\n")
print(res[res$sweep == "snr" & res$value %in% c(55:61, Inf), ], digits = 3, row.names = FALSE)
cat("\nRossler x1-dot, n sweep (49 dB), 10^4 to 10^5:\n")
print(res[res$sweep == "n" & res$value >= 4, ], digits = 3, row.names = FALSE)
cat("\nRossler x1-dot final-model term frequencies at 61 dB and Inf (frequency >= 0.05):\n")
print(fr[fr$sweep == "snr" & fr$value %in% c(61, Inf) & fr$frequency >= 0.05, ], row.names = FALSE)
cat("\nDadras x1-dot final-model term frequencies at n = 10^4 and 10^5 (frequency >= 0.05):\n")
print(fd[fd$sweep == "n" & fd$value %in% c(4, 5) & fd$frequency >= 0.05, ], row.names = FALSE)
cat("\nAizawa x3-dot final-model frequencies of x2^2x3 and x2^2x3^2:\n")
print(fa[fa$term %in% c("x2^2x3", "x2^2x3^2") & ((fa$sweep == "n" & fa$value %in% c(3.7, 4.4, 5)) | (fa$sweep == "snr" & fa$value %in% c(50, 55, 61, Inf))), ], row.names = FALSE)
