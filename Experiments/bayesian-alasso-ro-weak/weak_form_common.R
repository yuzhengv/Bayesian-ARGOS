# Shared code for the weak-form (integral) benchmark arms:
#   Bayesian-ARGOS (Integration)  -- Experiments/bayesian-alasso-ro-weak/weak_{n,snr}.R
#   SINDy (Integration)           -- Pysindy/weak/weak_{n,snr}.R
#
# The Savitzky-Golay step of build_design_matrix is replaced by the weak-form
# library of PySINDy (WeakPDELibrary, Messenger & Bortz 2021):
#   row k:  b_k = -int phi_k' x dt      against   G_k = int phi_k Theta(x) dt
# on K subdomains of half-width H_XT (time units) with test functions
# (1 - s^2)^P.  Subdomain centres are drawn at random by PySINDy, so NumPy is
# seeded per task (the same seed in both arms gives both arms the same
# subdomains).  Everything downstream is unchanged: the two adaptive-lasso
# passes and stan_glm of bayesian_alasso_ro(), or STLSQ at the benchmark
# threshold for SINDy.
#
# Requirements: pysindy >= 1.7.3 (WeakPDELibrary with K / H_xt / p and
# convert_u_dot_integral); tested with pysindy 2.1.0.

suppressPackageStartupMessages({
  library(signal); library(magrittr); library(tidyverse); library(glmnet)
  library(Matrix); library(plyr); library(rstanarm); library(reticulate); library(parallel)
})

# Repository root: ARGOS_ROOT, else the Hamilton path of the other experiments,
# else three levels up from the system folder the job runs in.
find_argos_root <- function() {
  cand <- c(Sys.getenv("ARGOS_ROOT"), "/nobackup/qtzk83/Projects/Bayesian-ARGOS",
            normalizePath("../../..", mustWork = FALSE), normalizePath("../../../..", mustWork = FALSE))
  for (p in cand) if (nzchar(p) && file.exists(file.path(p, "R", "argos_files.R"))) return(normalizePath(p))
  stop("cannot locate the Bayesian-ARGOS repository root (set ARGOS_ROOT)")
}
ARGOS_ROOT <- find_argos_root()
setwd(ARGOS_ROOT)
if (nzchar(Sys.getenv("ARGOS_PYTHON"))) reticulate::use_python(Sys.getenv("ARGOS_PYTHON"), required = TRUE)
source("./R/argos_files.R")
source("./R/bayesian_alasso_ro.R")
source("./Pysindy/src/pysindy_exp_fun.R")      # process_coefficients (SINDy result format)
source_python("./DataGeneration/ode_auto.py")
ps <- import("pysindy")
np <- import("numpy")
PYSINDY_VERSION <- tryCatch(ps$`__version__`, error = function(e) "unknown")

# ------------------------------------------------------------ weak-form settings
H_XT <- as.numeric(Sys.getenv("WEAK_H", "0.25"))     # half-width of each subdomain, time units (support 2 H = 0.5 = 50 samples at dt 0.01)
P_TEST <- as.integer(Sys.getenv("WEAK_P", "8"))      # test-function degree (1 - s^2)^p
K_FRAC <- as.numeric(Sys.getenv("WEAK_KFRAC", "0.25"))  # number of subdomains K = ceil(K_FRAC * n): rows scale with n

weak_K <- function(n) as.integer(max(10L, ceiling(K_FRAC * n)))

# ------------------------------------------------------------ systems (verbatim from the benchmark drivers)
SYSTEMS <- list(
  aizawa = list(
    coeff = list(list(-3.5, -0.7, 1), list(3.5, -0.7, 1), list(0.95, 0.65, 0.1, -1 / 3, -0.25, -1, -0.25, -1)),
    names = list(list("x2", "x1", "x1x3"), list("x1", "x2", "x2x3"), list("x3", "", "x1^3x3", "x3^3", "x1^2x3", "x1^2", "x2^2x3", "x2^2")),
    init = function(num_init, seed) {           # Experiments/bayesian-alasso-ro/aizawa/system_{n,snr}.R
      set.seed(seed); out <- list()
      for (i in 1:num_init) { x <- runif(1, -2, 2); y <- runif(1, -2, 2); z <- runif(1, -1, 2); out[[i]] <- c(x, y, z) }
      out
    }),
  dadras = list(
    coeff = list(list(1, -3, 2.7), list(1.7, -1, 1), list(2, -9)),
    names = list(list("x2", "x1", "x2x3"), list("x2", "x1x3", "x3"), list("x1x2", "x3")),
    init = function(num_init, seed) {           # Experiments/bayesian-alasso-ro/dadras/dadras_{n,snr}.R
      set.seed(seed); out <- list()
      for (i in 1:num_init) out[[i]] <- runif(3, min = -4, max = 4)
      out
    }),
  rossler = list(
    coeff = list(list(-1, -1), list(1, 0.2), list(0.2, 1, -5.7)),
    names = list(list("x2", "x3"), list("x1", "x2"), list("", "x1x3", "x3")),
    init = function(num_init, seed) {           # Experiments/bayesian-alasso-ro/rossler/rossler_{n,snr}.R
      set.seed(seed); out <- list()
      for (i in 1:num_init) { x <- runif(1, -10, 10); y <- runif(1, -10, 10); z <- runif(1, 0, 20); out[[i]] <- c(x, y, z) }
      out
    })
)

simulate_system <- function(system, n, dt, init_conditions, snr) {
  s <- SYSTEMS[[system]]
  generate_noisy_dynamical_systems_pyversion(variable_coeff = s$coeff, variable_names = s$names,
                                             n = as.integer(n), dt = dt, init_conditions = as.list(init_conditions), snr = snr)
}

# ------------------------------------------------------------ weak library
# PySINDy weak library on the uniform time grid of the record.
weak_library <- function(n, dt, library_degree) {
  t_grid <- np$arange(0, n * dt - dt / 2, dt)
  poly <- ps$feature_library$PolynomialLibrary(degree = as.integer(library_degree), include_bias = TRUE)
  ps$feature_library$WeakPDELibrary(function_library = poly, spatiotemporal_grid = t_grid,
                                    K = weak_K(n), H_xt = H_XT, p = P_TEST)
}

# PySINDy feature name ("x0^2 x2") -> the R pipeline's name ("x1^2x3"):
# shift the variable index, sort the tokens alphabetically, concatenate.
py_to_r_name <- function(nm) {
  toks <- strsplit(nm, " ")[[1]]
  toks <- sapply(toks, function(tk) { v <- as.integer(sub("^x([0-9]+).*$", "\\1", tk)) + 1L
    pw <- if (grepl("\\^", tk)) sub("^.*\\^", "", tk) else NA
    if (is.na(pw)) paste0("x", v) else paste0("x", v, "^", pw) })
  paste(sort(unname(toks)), collapse = "")
}
monomial_degree_py <- function(nm) sum(sapply(strsplit(nm, " ")[[1]], function(tk) if (grepl("\\^", tk)) as.integer(sub("^.*\\^", "", tk)) else 1L))

# Weak-form design in the layout bayesian_alasso_ro() consumes (sorted_theta,
# monomial_orders, xdot_filtered): rows normalised by int phi_k dt so that the
# constant column is 1 (dropped; the intercept is the regression's own) and
# the coefficients keep their scale; columns renamed and ordered by degree.
build_weak_design <- function(x_t, dt, library_degree, np_seed) {
  n <- nrow(x_t)
  np$random$seed(as.integer(np_seed))
  lib <- weak_library(n, dt, library_degree)
  lib$fit(x_t)
  theta_w <- lib$transform(x_t)
  b_w <- lib$convert_u_dot_integral(x_t)
  names_py <- unlist(lib$get_feature_names())
  bias <- which(names_py == "1")
  stopifnot(length(bias) == 1)
  scale <- theta_w[, bias]
  theta_w <- theta_w[, -bias, drop = FALSE] / scale
  b_w <- b_w / scale
  names_py <- names_py[-bias]
  deg <- sapply(names_py, monomial_degree_py)
  ord <- order(deg, names_py)
  theta <- data.frame(theta_w[, ord, drop = FALSE])
  colnames(theta) <- sapply(names_py[ord], py_to_r_name)
  colnames(b_w) <- paste0("xdot", seq_len(ncol(b_w)))
  list(sorted_theta = theta, monomial_orders = unname(deg[ord]), xdot_filtered = b_w,
       K = nrow(theta_w), H_xt = H_XT, p = P_TEST)
}

# One Bayesian-ARGOS (Integration) fit, same return layout as the SG driver.
run_bayesian_argos_weak <- function(system, n, init_conditions, dt, snr, library_degree, library_type,
                                    state_var_deriv, ci_level, bayesian_ncpus, np_seed) {
  xn <- simulate_system(system, n, dt, init_conditions, snr)
  design_matrix <- build_weak_design(xn, dt, library_degree, np_seed)
  start_time <- Sys.time()
  out <- bayesian_alasso_ro(design_matrix = design_matrix, target_distribution = gaussian(), library_type = library_type,
                            state_var_deriv = state_var_deriv, ci_level = ci_level, bayesian_ncpus = bayesian_ncpus)
  list(argos_bi_output = out, run_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
       weak = list(K = design_matrix$K, H_xt = H_XT, p = P_TEST, np_seed = np_seed))
}

# One SINDy (Integration) fit: PySINDy with the weak library and the benchmark
# STLSQ threshold, on the raw samples; same return layout as the SG driver.
run_sindy_weak <- function(system, n, init_conditions, dt, snr, library_degree, threshold, np_seed) {
  xn <- simulate_system(system, n, dt, init_conditions, snr)
  np$random$seed(as.integer(np_seed))
  lib <- weak_library(n, dt, library_degree)
  start_time <- Sys.time()
  model <- ps$SINDy(optimizer = ps$STLSQ(threshold = threshold), feature_library = lib)
  model_fit <- model$fit(x = xn, t = dt)
  run_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  feature_names <- model_fit$get_feature_names()
  coefficients <- model_fit$coefficients()
  list(exp_reults = process_coefficients(coefficients, feature_names), run_time = run_time,
       weak = list(K = weak_K(n), H_xt = H_XT, p = P_TEST, np_seed = np_seed))
}

np_seed_for <- function(seed, cond_index, trial) as.integer(seed * 1000 + cond_index * 101 + trial)
