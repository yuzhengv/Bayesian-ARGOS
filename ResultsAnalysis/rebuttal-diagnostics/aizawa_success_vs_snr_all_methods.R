# Aizawa success rate vs SNR for Bayesian-ARGOS, ARGOS and SINDy from the
# stored 100-trial benchmark results (n = 5000), plus a per-term failure
# decomposition of the ARGOS third equation.  Reproduces the table behind the
# right-hand panel of the manuscript's Fig. 3a with the paper's own evaluation
# functions (MethodsEvaluation/, Pysindy/utils/), exactly as
# ResultsAnalysis/success-rate-plots/aizawa_results_compare.ipynb does.
#
# Outputs (this folder):
#   aizawa_success_vs_snr_all_methods.csv  -- success rate per SNR and method
#   aizawa_eq3_argos_vs_snr.csv            -- ARGOS x3-dot: exact / missing /
#                                             spurious rates, bootstrap-active terms
# Usage: Rscript aizawa_success_vs_snr_all_methods.R   (from any directory)

here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
root <- normalizePath(file.path(here, "..", ".."))
setwd(root)
suppressPackageStartupMessages({ library(tidyverse); library(stringr) })
source("./MethodsEvaluation/results_processing_for_analysis_ignore_errors.R")
source("./Pysindy/utils/pysindy_results_processing_for_analysis.R")

dynamical_system_name <- "aizawa"
true_terms_1 <- list("x2", "x1", "x1x3")
true_terms_2 <- list("x1", "x2", "x2x3")
true_terms_3 <- list("x3", "(Intercept)", "x1^3x3", "x3^3", "x1^2x3", "x1^2", "x2^2x3", "x2^2")
true_terms_1_py <- list("x1", "x0", "x0 x2")
true_terms_2_py <- list("x0", "x1", "x1 x2")
true_terms_3_py <- list("x2", "1", "x0^3 x2", "x2^3", "x0^2 x2", "x0^2", "x1^2 x2", "x1^2")

# --- success rate per SNR, all three methods (all three equations correct) ---
invisible(capture.output(tab <- generate_total_success_rate_table(
  root_path = root, experiment_name_list = list("argos-alasso", "bayesian-alasso-ro"),
  method_list = list("argos-alasso", "bayesian-alasso-ro"), function_number = 3, num_init = 100,
  dynamical_system_name = dynamical_system_name, exp_n = FALSE, typical_pattern = "n5000",
  true_terms_1 = true_terms_1, true_terms_2 = true_terms_2, true_terms_3 = true_terms_3,
  true_terms_4 = NULL, start = 1, number_step = 1)))
invisible(capture.output(res_py <- check_dynamical_system_py(
  root_path = file.path(root, "Pysindy"), dynamical_system_name = dynamical_system_name,
  exp_n = FALSE, typical_pattern = "n5000", true_terms_1 = true_terms_1_py,
  true_terms_2 = true_terms_2_py, true_terms_3 = true_terms_3_py, true_terms_4 = NULL,
  function_number = 3)))
invisible(capture.output(tab_py <- build_successful_rate_table_py(
  dynamics_identification_results = res_py, num_init = 100, exp_n = FALSE,
  start = 1, number_step = 1, method = "pysindy")))
all <- rbind(tab, tab_py)
all <- all[all$Condition == "Correct", c("snr", "Value", "Model")]
w <- reshape(all, idvar = "snr", timevar = "Model", direction = "wide")
w <- w[order(as.numeric(w$snr)), ]
names(w) <- sub("^Value\\.", "", names(w))
w$snr[w$snr == "62"] <- "Inf"   # the evaluation code labels the noiseless level 62
write.csv(w, file.path(here, "aizawa_success_vs_snr_all_methods.csv"), row.names = FALSE)

# --- ARGOS x3-dot: how it fails at each SNR -----------------------------------
load("Experiments/argos-alasso/aizawa/results/aizawa_3_snr1_snr61_bo2000_se100_alasso_ridge_n5000.RData")
truth <- c("(Intercept)", "x3", "x1^2", "x2^2", "x1^2x3", "x2^2x3", "x3^3", "x1^3x3")
rows <- list()
for (nm in names(snr_output)) {
  L <- snr_output[[nm]]
  d <- do.call(rbind, lapply(L, function(t) {
    o <- t$id_result$argos_output; if (is.null(o)) return(NULL)
    im <- o$identified_model[, 1]; ci <- o$ci
    active <- colnames(ci)[!(ci[1, ] == 0 & ci[2, ] == 0)]   # terms with a non-degenerate bootstrap interval
    sel <- names(im)[im != 0]
    data.frame(exact = setequal(sel, truth), missing = !all(truth %in% sel),
               miss_x1_2x3 = !("x1^2x3" %in% sel), miss_x2_2x3 = !("x2^2x3" %in% sel), miss_x1_3x3 = !("x1^3x3" %in% sel),
               spurious = length(setdiff(sel, truth)), n_active = length(active),
               n_active_spurious = length(setdiff(active, truth)))
  }))
  rows[[nm]] <- data.frame(snr = sub("snr=", "", nm), trials = nrow(d), success = mean(d$exact),
                           any_missing = mean(d$missing), miss_x1_2x3 = mean(d$miss_x1_2x3),
                           miss_x2_2x3 = mean(d$miss_x2_2x3), miss_x1_3x3 = mean(d$miss_x1_3x3),
                           any_spurious = mean(d$spurious > 0), mean_spurious = mean(d$spurious),
                           active_terms = mean(d$n_active), active_spurious = mean(d$n_active_spurious))
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(here, "aizawa_eq3_argos_vs_snr.csv"), row.names = FALSE)
keep <- out$snr %in% c("40", "45", "49", "50", "53", "55", "57", "59", "61", "Inf")
print(w[as.numeric(sub("Inf", "1000", w$snr)) >= 40, ], row.names = FALSE)
print(out[keep, ], digits = 2, row.names = FALSE)
