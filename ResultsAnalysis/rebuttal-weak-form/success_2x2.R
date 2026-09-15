# 2 x 2 comparison of response construction x selection stage:
#   Bayesian-ARGOS (SG)          Experiments/bayesian-alasso-ro/<system>/results        (stored benchmark)
#   Bayesian-ARGOS (Integration) Experiments/bayesian-alasso-ro-weak/<system>/results   (weak-form arm)
#   SINDy (SG)                   Pysindy/exp/<system>/results                            (stored benchmark)
#   SINDy (Integration)          Pysindy/weak/exp/<system>/results                       (weak-form arm)
# Success rate (all three equations exactly recovered) against n (49 dB) and
# against SNR (n = 5000), with the paper's own evaluation functions
# (MethodsEvaluation/, Pysindy/utils/), exactly as the success-rate notebooks.
# The evaluation code dispatches on a fixed list of method names, so the weak
# arms are loaded under the label of their SG counterpart and relabelled here.
#
# Usage: Rscript success_2x2.R [systems="aizawa,dadras,rossler"] [num_init=100]
# Output: results/success_2x2_n.csv, results/success_2x2_snr.csv (long format:
#         system, sweep, x, method, response, selection, success)

args <- commandArgs(trailingOnly = TRUE)
SYSTEMS_RUN <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else c("aizawa", "dadras", "rossler")
NUM_INIT <- if (length(args) >= 2) as.integer(args[2]) else 100L
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
root <- normalizePath(file.path(here, "..", ".."))
setwd(root)
suppressPackageStartupMessages({ library(tidyverse); library(stringr) })
source("./MethodsEvaluation/results_processing_for_analysis_ignore_errors.R")
source("./Pysindy/utils/pysindy_results_processing_for_analysis.R")
dir.create(file.path(here, "results"), showWarnings = FALSE)

TRUTH <- list(
  aizawa = list(r = list(list("x2", "x1", "x1x3"), list("x1", "x2", "x2x3"), list("x3", "(Intercept)", "x1^3x3", "x3^3", "x1^2x3", "x1^2", "x2^2x3", "x2^2")),
                py = list(list("x1", "x0", "x0 x2"), list("x0", "x1", "x1 x2"), list("x2", "1", "x0^3 x2", "x2^3", "x0^2 x2", "x0^2", "x1^2 x2", "x1^2"))),
  dadras = list(r = list(list("x2", "x1", "x2x3"), list("x2", "x1x3", "x3"), list("x1x2", "x3")),
                py = list(list("x1", "x0", "x1 x2"), list("x1", "x0 x2", "x2"), list("x0 x1", "x2"))),
  rossler = list(r = list(list("x2", "x3"), list("x1", "x2"), list("(Intercept)", "x1x3", "x3")),
                 py = list(list("x1", "x2"), list("x0", "x1"), list("1", "x0 x2", "x2")))
)
SWEEPS <- list(n = list(exp_n = TRUE, pattern = "snr49", start = 2, step = 0.1),
               snr = list(exp_n = FALSE, pattern = "n5000", start = 1, step = 1))

load_r <- function(system, experiment, sweep) {
  tt <- TRUTH[[system]]$r; s <- SWEEPS[[sweep]]
  out <- tryCatch({
    invisible(capture.output(tab <- generate_total_success_rate_table(
      root_path = root, experiment_name_list = list(experiment), method_list = list("bayesian-alasso-ro"),
      function_number = 3, num_init = NUM_INIT, dynamical_system_name = system, exp_n = s$exp_n, typical_pattern = s$pattern,
      true_terms_1 = tt[[1]], true_terms_2 = tt[[2]], true_terms_3 = tt[[3]], true_terms_4 = NULL, start = s$start, number_step = s$step)))
    tab
  }, error = function(e) { message(sprintf("  %s %s %s: %s", system, experiment, sweep, conditionMessage(e))); NULL })
  if (is.null(out)) return(NULL)
  out <- out[out$Condition == "Correct", ]
  data.frame(system = system, sweep = sweep, x = as.character(out[[1]]), success = as.numeric(out$Value),
             method = if (experiment == "bayesian-alasso-ro") "Bayesian-ARGOS (SG)" else "Bayesian-ARGOS (Integration)", stringsAsFactors = FALSE)
}
load_py <- function(system, rootpath, sweep) {
  tt <- TRUTH[[system]]$py; s <- SWEEPS[[sweep]]
  out <- tryCatch({
    invisible(capture.output(res <- check_dynamical_system_py(
      root_path = file.path(root, rootpath), dynamical_system_name = system, exp_n = s$exp_n, typical_pattern = s$pattern,
      true_terms_1 = tt[[1]], true_terms_2 = tt[[2]], true_terms_3 = tt[[3]], true_terms_4 = NULL, function_number = 3)))
    invisible(capture.output(tab <- build_successful_rate_table_py(dynamics_identification_results = res, num_init = NUM_INIT,
      exp_n = s$exp_n, start = s$start, number_step = s$step, method = "pysindy")))
    tab
  }, error = function(e) { message(sprintf("  %s %s %s: %s", system, rootpath, sweep, conditionMessage(e))); NULL })
  if (is.null(out)) return(NULL)
  out <- out[out$Condition == "Correct", ]
  data.frame(system = system, sweep = sweep, x = as.character(out[[1]]), success = as.numeric(out$Value),
             method = if (rootpath == "Pysindy") "SINDy (SG)" else "SINDy (Integration)", stringsAsFactors = FALSE)
}

rows <- list()
for (system in SYSTEMS_RUN) for (sweep in names(SWEEPS)) {
  rows[[length(rows) + 1]] <- load_r(system, "bayesian-alasso-ro", sweep)
  rows[[length(rows) + 1]] <- load_r(system, "bayesian-alasso-ro-weak", sweep)
  rows[[length(rows) + 1]] <- load_py(system, "Pysindy", sweep)
  rows[[length(rows) + 1]] <- load_py(system, "Pysindy/weak", sweep)
}
all <- do.call(rbind, rows)
all$response <- ifelse(grepl("Integration", all$method), "Integration", "SG")
all$selection <- ifelse(grepl("^Bayesian", all$method), "Bayesian-ARGOS", "SINDy")
# the evaluation code labels the noiseless level as snr_end + 1 (62)
all$x[all$sweep == "snr" & all$x == "62"] <- "Inf"
for (sw in names(SWEEPS)) {
  d <- all[all$sweep == sw, ]
  write.csv(d, file.path(here, "results", sprintf("success_2x2_%s.csv", sw)), row.names = FALSE)
  w <- reshape(d[, c("system", "x", "method", "success")], idvar = c("system", "x"), timevar = "method", direction = "wide")
  names(w) <- sub("^success\\.", "", names(w)); w <- w[order(w$system, suppressWarnings(as.numeric(sub("Inf", "1e9", w$x)))), ]
  cat("\n==", sw, "sweep\n"); print(w, row.names = FALSE, digits = 2)
}
