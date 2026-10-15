# Referee 1 Point 3 (R1.3): shared definitions for the spurious-terms evidence.
#
# Sourced by spurious_terms_table.R, make_fig6.R, make_si_term_selection.R and
# make_si_spurious_by_system.R. Holds the paths, the seven benchmark systems
# with their true terms in the two naming schemes (the R candidate library of
# ARGOS and Bayesian-ARGOS; the PySINDy library of SINDy), the three methods
# with their stored-experiment folders, the ring values of the circular plots,
# and a loader that returns the per-trial identified terms through the
# repository's own evaluation modules, unchanged.
#
# The four notebook modules (vary-n / vary-SNR for the ARGOS family and for
# PySINDy) define functions with the same names, so each is sourced into its
# own environment and called from there.

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
})

this_file <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
here <- dirname(this_file)
root_path <- normalizePath(file.path(here, "..", ".."))
manuscript_root <- normalizePath(file.path(root_path, "..", ".."), mustWork = FALSE)
dir.create(file.path(here, "results"), showWarnings = FALSE)
dir.create(file.path(here, "figures"), showWarnings = FALSE)

source_into_env <- function(f) {
    e <- new.env(parent = globalenv())
    suppressMessages(suppressWarnings(sys.source(f, envir = e)))
    e
}
module_env <- list()
get_module <- function(family, sweep) {
    key <- paste(family, sweep, sep = "_")
    if (is.null(module_env[[key]])) {
        f <- switch(key,
            r_n = "MethodsEvaluation/evaluation-functions-for-plots/functions_for_analysis_of_identified_terms_vary_n_enhanced.R",
            r_snr = "MethodsEvaluation/evaluation-functions-for-plots/functions_for_analysis_of_identified_terms_vary_snr_enhanced.R",
            py_n = "Pysindy/utils/utils-functions-for-plots/pysindy_functions_for_analysis_of_identified_terms_vary_n.R",
            py_snr = "Pysindy/utils/utils-functions-for-plots/pysindy_functions_for_analysis_of_identified_terms_vary_snr.R"
        )
        module_env[[key]] <<- source_into_env(file.path(root_path, f))
    }
    module_env[[key]]
}
quietly <- function(expr) { # the loaders print every file name and list length
    sink(tempfile()); on.exit(sink()); force(expr)
}

# ---------------------------------------------------------------- methods
METHODS <- c("Bayesian-ARGOS", "ARGOS", "SINDy")
METHOD_INFO <- list(
    "Bayesian-ARGOS" = list(family = "r", experiment = "bayesian-alasso-ro"),
    "ARGOS" = list(family = "r", experiment = "argos-alasso"),
    "SINDy" = list(family = "py", experiment = "pysindy")
)

# ---------------------------------------------------------------- systems
# True terms in the naming of the stored results: the R candidate library
# (ARGOS, Bayesian-ARGOS) and the PySINDy library (SINDy, zero-based state
# index), as in the success-rate notebooks ResultsAnalysis/success-rate-plots/.
SYSTEMS <- c("lorenz", "thomas", "rossler", "dadras", "aizawa", "sprott", "halvorsen")
SYSTEM_LABEL <- c(lorenz = "Lorenz", thomas = "Thomas", rossler = "Rössler", dadras = "Dadras",
                  aizawa = "Aizawa", sprott = "Sprott", halvorsen = "Halvorsen")
TRUE_TERMS <- list(
    lorenz = list(
        r = list(c("x1", "x2"), c("x1", "x1x3", "x2"), c("x1x2", "x3")),
        py = list(c("x0", "x1"), c("x0", "x0 x2", "x1"), c("x0 x1", "x2"))
    ),
    thomas = list(
        r = list(c("sin_x2", "x1"), c("sin_x3", "x2"), c("sin_x1", "x3")),
        py = list(c("sin(1 x1)", "x0"), c("sin(1 x2)", "x1"), c("sin(1 x0)", "x2"))
    ),
    rossler = list(
        r = list(c("x2", "x3"), c("x1", "x2"), c("(Intercept)", "x1x3", "x3")),
        py = list(c("x1", "x2"), c("x0", "x1"), c("1", "x0 x2", "x2"))
    ),
    dadras = list(
        r = list(c("x2", "x1", "x2x3"), c("x2", "x1x3", "x3"), c("x1x2", "x3")),
        py = list(c("x1", "x0", "x1 x2"), c("x1", "x0 x2", "x2"), c("x0 x1", "x2"))
    ),
    aizawa = list(
        r = list(c("x2", "x1", "x1x3"), c("x1", "x2", "x2x3"),
                 c("x3", "(Intercept)", "x1^3x3", "x3^3", "x1^2x3", "x1^2", "x2^2x3", "x2^2")),
        py = list(c("x1", "x0", "x0 x2"), c("x0", "x1", "x1 x2"),
                  c("x2", "1", "x0^3 x2", "x2^3", "x0^2 x2", "x0^2", "x1^2 x2", "x1^2"))
    ),
    sprott = list(
        r = list(c("x2", "x1x2", "x1x3"), c("(Intercept)", "x1^2", "x2x3"), c("x1", "x1^2", "x2^2")),
        py = list(c("x1", "x0 x1", "x0 x2"), c("1", "x0^2", "x1 x2"), c("x0", "x0^2", "x1^2"))
    ),
    halvorsen = list(
        r = list(c("x1", "x2", "x3", "x2^2"), c("x2", "x3", "x1", "x3^2"), c("x3", "x1", "x2", "x1^2")),
        py = list(c("x0", "x1", "x2", "x1^2"), c("x1", "x2", "x0", "x2^2"), c("x2", "x0", "x1", "x0^2"))
    )
)

# ---------------------------------------------------------------- rings
# Ring values of the circular plots. n sweep: one set for every system
# (log10 n). SNR sweep: the values the submitted SI figures used (the
# vary-snr notebooks' snr_list), Inf = the noiseless case.
RINGS_N <- c(2.8, 3, 3.5, 4, 5)
RINGS_SNR <- list(
    lorenz = c(20, 25, 30, 49, 61), thomas = c(17, 25, 30, 49, Inf),
    rossler = c(18, 20, 25, 37, 61, Inf), dadras = c(10, 13, 20, 25, 37, 61),
    aizawa = c(20, 25, 37, 50, 60, 61), sprott = c(24, 30, 37, 49, 61, Inf),
    halvorsen = c(12, 13, 20, 25, 37, 61)
)

# ---------------------------------------------------------------- loader
# Per-trial identified terms of the final models, through the notebooks'
# enhance_total_results_summary_df (R family) / _py (SINDy), unchanged.
# Returns: system, method, sweep, x (log10 n, or SNR in dB with Inf),
# equation (1..3), trial (1..100), terms (list of raw term names).
load_identified_terms <- function(system, method, sweep) {
    info <- METHOD_INFO[[method]]
    mod <- get_module(info$family, sweep)
    exp_n <- sweep == "n"
    pattern <- if (exp_n) "snr49" else "n5000"
    df <- if (info$family == "r") {
        quietly(mod$enhance_total_results_summary_df(
            root_path = root_path, experiment_name_list = list(info$experiment),
            method_list = list(info$experiment), function_number = 3, num_init = 100,
            dynamical_system_name = system, exp_n = exp_n, typical_pattern = pattern
        ))
    } else {
        quietly(mod$enhance_total_results_summary_df_py(
            root_path = file.path(root_path, "Pysindy"), experiment_name_list = list(info$experiment),
            method_list = list(info$experiment), num_init = 100,
            dynamical_system_name = system, exp_n = exp_n, typical_pattern = pattern
        ))
    }
    xcol <- if (exp_n) df$eta else df$snr
    tibble(
        system = system, method = method, sweep = sweep,
        x = as.numeric(xcol), equation = as.integer(df$function_id),
        trial = as.integer(sub(".*\\.", "", df$run_id)),
        terms = strsplit(df$identified_terms, ",\\s*")
    )
}

# Raw term name -> LaTeX label, through the module's own converter, so that
# the two naming schemes meet in one label ("$x_{1}x_{3}$"); the intercept is
# "$C$" in both.
term_label <- function(terms, family) {
    mod <- get_module(family, "n")
    out <- vapply(terms, function(t) mod$convert_term_to_latex(t), character(1))
    out[out == "(Intercept)"] <- "$C$"
    unname(out)
}

fmt_x <- function(x, sweep) {
    if (sweep == "n") sprintf("10^%g", x) else ifelse(is.infinite(x), "∞", sprintf("%g dB", x))
}
