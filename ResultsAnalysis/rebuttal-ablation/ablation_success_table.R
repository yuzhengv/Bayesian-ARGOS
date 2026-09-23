# Referee 1 Point 1 (a, b): numeric summary of the screening-stage ablation.
#
# Reads the stored 100-trial results of the six screening variants
# (Experiments/bayesian-alasso-{ro,or,oo,rr,single-ols,single-ridge}) for the
# Lorenz, Rössler and Thomas systems with the evaluation function the
# success-rate notebooks use (generate_total_success_rate_table, unchanged),
# and writes
#   results/ablation_success_rates.csv   long table: system, sweep, arm, x, success
#   tables/ablation_summary.md           internal reference: 80 % crossing, reference
#                                        points, plateau, difference to the baseline
# Usage: Rscript ablation_success_table.R   (from this folder or anywhere)

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

experiment_name_list <- list(
    "bayesian-alasso-ro", "bayesian-alasso-or", "bayesian-alasso-oo",
    "bayesian-alasso-rr", "bayesian-alasso-single-ols", "bayesian-alasso-single-ridge"
)
arm_label <- c(
    "bayesian-alasso-ro" = "Ridge-OLS (baseline)",
    "bayesian-alasso-or" = "OLS-Ridge",
    "bayesian-alasso-oo" = "OLS-OLS",
    "bayesian-alasso-rr" = "Ridge-Ridge",
    "bayesian-alasso-single-ols" = "Single OLS",
    "bayesian-alasso-single-ridge" = "Single Ridge"
)

systems <- list(
    lorenz = list(t1 = list("x1", "x2"), t2 = list("x1", "x1x3", "x2"), t3 = list("x1x2", "x3")),
    rossler = list(t1 = list("x2", "x3"), t2 = list("x1", "x2"), t3 = list("(Intercept)", "x1x3", "x3")),
    thomas = list(t1 = list("sin_x2", "x1"), t2 = list("sin_x3", "x2"), t3 = list("sin_x1", "x3"))
)

one_sweep <- function(sys, tt, exp_n) {
    tab <- generate_total_success_rate_table(
        root_path = root_path,
        experiment_name_list = experiment_name_list,
        method_list = experiment_name_list,
        function_number = 3,
        num_init = 100,
        dynamical_system_name = sys,
        exp_n = exp_n,
        typical_pattern = if (exp_n) "snr49" else "n5000",
        true_terms_1 = tt$t1, true_terms_2 = tt$t2, true_terms_3 = tt$t3, true_terms_4 = NULL,
        start = if (exp_n) 2 else 1,
        number_step = if (exp_n) 0.1 else 1
    )
    tab %>%
        transmute(
            system = sys,
            sweep = if (exp_n) "n" else "snr",
            arm = unname(arm_label[as.character(Model)]),
            x = if (exp_n) as.numeric(eta) else ifelse(as.numeric(snr) == 62, Inf, as.numeric(snr)), # log10(n); or SNR in dB (the stored level 62 is the noiseless case)
            success = as.numeric(Value)
        )
}

long <- bind_rows(lapply(names(systems), function(sys) {
    bind_rows(one_sweep(sys, systems[[sys]], TRUE), one_sweep(sys, systems[[sys]], FALSE))
}))
long$arm <- factor(long$arm, levels = unname(arm_label))
write.csv(long, file.path(here, "results", "ablation_success_rates.csv"), row.names = FALSE)

# ---- summaries ---------------------------------------------------------------
first_cross <- function(x, y, level = 0.8) {
    ok <- which(y >= level & is.finite(x))
    if (length(ok) == 0) NA_real_ else x[min(ok)]
}
fmt_x <- function(v, sweep) {
    sweep <- rep_len(sweep, length(v)) # a scalar sweep would otherwise make ifelse() return one value
    ifelse(is.na(v), "never", ifelse(sweep == "n", sprintf("10^%.1f", v), sprintf("%g dB", v)))
}

summ <- long %>%
    group_by(system, sweep, arm) %>%
    arrange(x, .by_group = TRUE) %>%
    summarise(
        cross80 = { o <- order(x); i <- which(success[o] >= 0.8 & is.finite(x[o]))[1]; if (is.na(i)) NA_real_ else x[o][i] },
        plateau = if (first(sweep) == "n") mean(success[x >= 4]) else mean(success[is.finite(x) & x >= 49]),
        max_finite = max(success[is.finite(x)]),
        inf = if (first(sweep) == "snr") success[is.infinite(x)][1] else NA_real_,
        .groups = "drop"
    )

ref_pts <- list(n = c(2.5, 3, 3.5, 4, 5), snr = c(13, 25, 37, 49, 61, Inf))
ref <- long %>%
    filter((sweep == "n" & abs(x - round(x * 2) / 2) < 1e-6 & x %in% ref_pts$n) |
        (sweep == "snr" & x %in% ref_pts$snr)) %>%
    mutate(xlab = fmt_x(x, sweep)) %>%
    select(system, sweep, arm, xlab, success)
ref <- bind_rows(lapply(split(ref, ref$sweep), function(d) pivot_wider(d, names_from = xlab, values_from = success)))

# per grid point difference to the baseline, then the mean and the worst case
diff <- long %>%
    filter(is.finite(x)) %>%
    group_by(system, sweep, x) %>%
    mutate(d = success - success[arm == "Ridge-OLS (baseline)"]) %>%
    ungroup() %>%
    filter(arm != "Ridge-OLS (baseline)") %>%
    group_by(system, sweep, arm) %>%
    summarise(
        mean_diff = mean(d), min_diff = min(d), at_min = x[which.min(d)],
        n_points = n(), n_below = sum(d < -0.05), n_above = sum(d > 0.05),
        .groups = "drop"
    )

md <- c(
    "# Screening-stage ablation: numeric summary (internal reference)",
    "",
    sprintf("Generated %s by `ablation_success_table.R` from the stored results in `Experiments/bayesian-alasso-*/{lorenz,rossler,thomas}/results/` (100 initial conditions per grid point; n sweep at 49 dB, SNR sweep at n = 5000; success = all three equations exactly recovered, the notebooks' criterion).", format(Sys.Date())),
    "",
    "## 80 % crossing, plateau and noiseless point",
    "",
    "cross80 = first grid value at which success >= 0.8; plateau = mean success over n >= 1e4 (n sweep) or 49-61 dB (SNR sweep); inf = success at SNR = infinity.",
    ""
)
for (sw in c("n", "snr")) {
    for (sys in names(systems)) {
        s <- summ %>% filter(system == sys, sweep == sw) %>% arrange(arm)
        md <- c(md, sprintf("### %s, %s sweep", str_to_title(sys), if (sw == "n") "n" else "SNR"), "",
            "| arm | cross80 | plateau | max (finite) | inf |", "|---|---|---|---|---|",
            sprintf("| %s | %s | %.2f | %.2f | %s |", s$arm, fmt_x(s$cross80, sw), s$plateau, s$max_finite,
                ifelse(is.na(s$inf), "-", sprintf("%.2f", s$inf))), "")
    }
}
md <- c(md, "## Success at reference points", "")
for (sw in c("n", "snr")) {
    for (sys in names(systems)) {
        r <- ref %>% filter(system == sys, sweep == sw) %>% arrange(arm)
        cols <- setdiff(names(r), c("system", "sweep", "arm"))
        md <- c(md, sprintf("### %s, %s sweep", str_to_title(sys), if (sw == "n") "n" else "SNR"), "",
            paste0("| arm | ", paste(cols, collapse = " | "), " |"),
            paste0("|---|", paste(rep("---", length(cols)), collapse = "|"), "|"),
            apply(r, 1, function(row) paste0("| ", row[["arm"]], " | ", paste(sprintf("%.2f", as.numeric(row[cols])), collapse = " | "), " |")), "")
    }
}
md <- c(md, "## Difference to the Ridge-OLS baseline over the finite grid", "",
    "mean_diff / min_diff = mean and most negative (arm minus baseline) over all finite grid points; at_min = where the minimum occurs (log10 n or dB); n_below / n_above = grid points more than 0.05 below / above the baseline.", "",
    "| system | sweep | arm | mean_diff | min_diff | at_min | points | below | above |", "|---|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %+.3f | %+.2f | %s | %d | %d | %d |", diff$system, diff$sweep, diff$arm, diff$mean_diff, diff$min_diff,
        fmt_x(diff$at_min, diff$sweep), diff$n_points, diff$n_below, diff$n_above), "")
writeLines(md, file.path(here, "tables", "ablation_summary.md"))
cat(paste(md, collapse = "\n"), "\n")
