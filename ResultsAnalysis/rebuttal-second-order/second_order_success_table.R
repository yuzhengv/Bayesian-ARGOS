# Referee 2 Point 2 (b) and the editor's third requirement: numeric summary of
# the three second-order benchmark systems (Duffing, forced Duffing, van der Pol),
# identified in augmented first-order form by Bayesian-ARGOS, ARGOS and SINDy.
#
# Reads the stored 100-trial results with the evaluation functions the
# success-rate notebooks use (generate_total_success_rate_table and the PySINDy
# counterparts, unchanged) and writes
#   results/second_order_success_rates.csv  long table: system, sweep, method, x, success
#   tables/second_order_summary.md          80 % crossing, plateau, noiseless point,
#                                           success at reference points
# Usage: Rscript second_order_success_table.R   (from this folder or anywhere)

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
})

this_file <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
here <- dirname(this_file)
root_path <- normalizePath(file.path(here, "..", ".."))
setwd(root_path)
source("./MethodsEvaluation/results_processing_for_analysis.R")
source("./Pysindy/utils/pysindy_results_processing_for_analysis.R")
dir.create(file.path(here, "results"), showWarnings = FALSE)
dir.create(file.path(here, "tables"), showWarnings = FALSE)

method_label <- c(
    "bayesian-alasso-ro" = "Bayesian-ARGOS",
    "argos-alasso" = "ARGOS",
    "pysindy" = "SINDy"
)

# True terms as the notebooks name them (R pipeline and PySINDy naming)
systems <- list(
    duffing = list(
        r = list(t1 = list("x2"), t2 = list("x1", "x2", "x1^3")),
        py = list(t1 = list("x1"), t2 = list("x0", "x1", "x0^3"))
    ),
    `forced-duffing` = list(
        r = list(t1 = list("x2"), t2 = list("x1", "x2", "x1^3", "cos_2t")),
        py = list(t1 = list("x1"), t2 = list("x0", "x1", "x0^3", "u0"))
    ),
    vdp = list(
        r = list(t1 = list("x2"), t2 = list("x1", "x2", "x1^2x2")),
        py = list(t1 = list("x1"), t2 = list("x0", "x1", "x0^2 x1"))
    )
)

one_sweep <- function(sys, tt, exp_n) {
    tab_r <- generate_total_success_rate_table(
        root_path = root_path,
        experiment_name_list = list("argos-alasso", "bayesian-alasso-ro"),
        method_list = list("argos-alasso", "bayesian-alasso-ro"),
        function_number = 2,
        num_init = 100,
        dynamical_system_name = sys,
        exp_n = exp_n,
        typical_pattern = if (exp_n) "snr49" else "n5000",
        true_terms_1 = tt$r$t1, true_terms_2 = tt$r$t2, true_terms_3 = NULL, true_terms_4 = NULL,
        start = if (exp_n) 2 else 1,
        number_step = if (exp_n) 0.1 else 1
    )
    py_res <- check_dynamical_system_py(
        root_path = file.path(root_path, "Pysindy"),
        dynamical_system_name = sys,
        exp_n = exp_n,
        typical_pattern = if (exp_n) "snr49" else "n5000",
        true_terms_1 = tt$py$t1, true_terms_2 = tt$py$t2, true_terms_3 = NULL, true_terms_4 = NULL,
        function_number = 2
    )
    tab_py <- build_successful_rate_table_py(
        dynamics_identification_results = py_res,
        num_init = 100,
        exp_n = exp_n,
        start = if (exp_n) 2 else 1,
        number_step = if (exp_n) 0.1 else 1,
        method = "pysindy"
    )
    tab <- bind_rows(tab_r, tab_py)
    tab %>%
        transmute(
            system = sys,
            sweep = if (exp_n) "n" else "snr",
            method = unname(method_label[as.character(Model)]),
            x = if (exp_n) as.numeric(eta) else ifelse(as.numeric(snr) == 62, Inf, as.numeric(snr)), # log10(n); or SNR in dB (stored level 62 = noiseless)
            success = as.numeric(Value)
        )
}

long <- bind_rows(lapply(names(systems), function(sys) {
    bind_rows(one_sweep(sys, systems[[sys]], TRUE), one_sweep(sys, systems[[sys]], FALSE))
}))
long$method <- factor(long$method, levels = unname(method_label))
write.csv(long, file.path(here, "results", "second_order_success_rates.csv"), row.names = FALSE)

fmt_x <- function(v, sweep) {
    sweep <- rep_len(sweep, length(v))
    ifelse(is.na(v), "never", ifelse(sweep == "n", sprintf("10^%.1f", v), sprintf("%g dB", v)))
}

summ <- long %>%
    group_by(system, sweep, method) %>%
    arrange(x, .by_group = TRUE) %>%
    summarise(
        cross80 = { o <- order(x); i <- which(success[o] >= 0.8 & is.finite(x[o]))[1]; if (is.na(i)) NA_real_ else x[o][i] },
        plateau = if (first(sweep) == "n") mean(success[x >= 4]) else mean(success[is.finite(x) & x >= 49]),
        max_finite = max(success[is.finite(x)]),
        min_after_cross = { o <- order(x); i <- which(success[o] >= 0.8 & is.finite(x[o]))[1]; if (is.na(i)) NA_real_ else min(success[o][seq(i, sum(is.finite(x)))]) },
        inf = if (first(sweep) == "snr") success[is.infinite(x)][1] else NA_real_,
        .groups = "drop"
    )

ref_pts <- list(n = c(2.5, 3, 3.5, 4, 5), snr = c(13, 25, 37, 49, 61, Inf))
ref <- long %>%
    filter((sweep == "n" & abs(x - round(x * 2) / 2) < 1e-6 & x %in% ref_pts$n) |
        (sweep == "snr" & x %in% ref_pts$snr)) %>%
    mutate(xlab = fmt_x(x, sweep)) %>%
    select(system, sweep, method, xlab, success)
ref <- bind_rows(lapply(split(ref, ref$sweep), function(d) pivot_wider(d, names_from = xlab, values_from = success)))

md <- c(
    "# Second-order systems: numeric summary (internal reference)",
    "",
    sprintf("Generated %s by `second_order_success_table.R` from the stored results in `Experiments/{argos-alasso,bayesian-alasso-ro}/{duffing,forced-duffing,vdp}/results/` and `Pysindy/exp/{duffing,forced-duffing,vdp}/results/` (100 initial conditions per grid point; n sweep at 49 dB, SNR sweep at n = 5000; success = both equations exactly recovered, the notebooks' criterion).", format(Sys.Date())),
    "",
    "## 80 % crossing, plateau and noiseless point",
    "",
    "cross80 = first grid value at which success >= 0.8; plateau = mean success over n >= 1e4 (n sweep) or 49-61 dB (SNR sweep); min_after_cross = lowest success at or beyond the crossing (finite grid); inf = success at SNR = infinity.",
    ""
)
for (sw in c("n", "snr")) {
    for (sys in names(systems)) {
        s <- summ %>% filter(system == sys, sweep == sw) %>% arrange(method)
        md <- c(md, sprintf("### %s, %s sweep", sys, if (sw == "n") "n" else "SNR"), "",
            "| method | cross80 | plateau | max (finite) | min after crossing | inf |", "|---|---|---|---|---|---|",
            sprintf("| %s | %s | %.2f | %.2f | %s | %s |", s$method, fmt_x(s$cross80, sw), s$plateau, s$max_finite,
                ifelse(is.na(s$min_after_cross), "-", sprintf("%.2f", s$min_after_cross)),
                ifelse(is.na(s$inf), "-", sprintf("%.2f", s$inf))), "")
    }
}
md <- c(md, "## Success at reference points", "")
for (sw in c("n", "snr")) {
    for (sys in names(systems)) {
        r <- ref %>% filter(system == sys, sweep == sw) %>% arrange(method)
        cols <- setdiff(names(r), c("system", "sweep", "method"))
        md <- c(md, sprintf("### %s, %s sweep", sys, if (sw == "n") "n" else "SNR"), "",
            paste0("| method | ", paste(cols, collapse = " | "), " |"),
            paste0("|---|", paste(rep("---", length(cols)), collapse = "|"), "|"),
            apply(r, 1, function(row) paste0("| ", row[["method"]], " | ", paste(sprintf("%.2f", as.numeric(row[cols])), collapse = " | "), " |")), "")
    }
}
writeLines(md, file.path(here, "tables", "second_order_summary.md"))
cat(paste(md, collapse = "\n"), "\n")
