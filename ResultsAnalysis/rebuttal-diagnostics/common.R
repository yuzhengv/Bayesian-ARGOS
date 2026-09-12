# Shared helpers for the rebuttal diagnostics (Referee 1 Point 6, Referee 2
# minor point 2), R version.  Everything re-uses the code that produced the
# manuscript benchmarks: data generation through
# `DataGeneration/ode_auto.py` (scipy odeint, called via reticulate exactly as
# in Experiments/*/system_snr.R), `build_design_matrix` (Savitzky-Golay
# smoothing / differentiation + polynomial library), the two adaptive-lasso
# passes of `double_regression`, and the `stan_glm` + posterior-interval rule
# of `bayesian_alasso_ro` (`R/argos_files.R`, `R/bayesian_alasso_ro.R`).
#
# Source this file from the repository root or from this folder; it locates
# the root itself.

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(signal)
  library(Metrics)
  library(lmtest)
  library(parallel)
  library(reticulate)
})
if (requireNamespace("rstanarm", quietly = TRUE)) suppressPackageStartupMessages(library(rstanarm))  # needed by bayes_stage only
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
PYTHON_BIN <- Sys.getenv("ARGOS_PYTHON", "/home/yuzhengz/anaconda3/envs/env_dgx_spark/bin/python")

find_root <- function() {
  cand <- c(getwd(), file.path(getwd(), "..", ".."), dirname(dirname(getwd())))
  for (p in cand) if (file.exists(file.path(p, "R", "argos_files.R"))) return(normalizePath(p))
  stop("cannot locate the Bayesian-ARGOS repository root")
}
ROOT <- find_root()
source(file.path(ROOT, "R", "argos_files.R"))
source(file.path(ROOT, "R", "bayesian_alasso_ro.R"))

DT <- 0.01
LIBRARY_DEGREE <- 5
SG_POLY_ORDER <- 4
NCORES <- as.integer(Sys.getenv("ARGOS_NCORES", max(1L, min(16L, parallel::detectCores() - 4L))))

# ------------------------------------------------------------------ systems
# Coefficients / term names exactly as in Experiments/bayesian-alasso-ro/*.
SYSTEMS <- list(
  aizawa = list(
    # coefficient / term-name lists in the format of ode_auto.py, copied from
    # Experiments/bayesian-alasso-ro/aizawa/system_snr.R
    coeff = list(list(-3.5, -0.7, 1), list(3.5, -0.7, 1), list(0.95, 0.65, 0.1, -1 / 3, -0.25, -1, -0.25, -1)),
    names = list(list("x2", "x1", "x1x3"), list("x1", "x2", "x2x3"),
                 list("x3", "", "x1^3x3", "x3^3", "x1^2x3", "x1^2", "x2^2x3", "x2^2")),
    rhs = function(t, x, p) {
      x1 <- x[1]; x2 <- x[2]; x3 <- x[3]
      list(c(
        -3.5 * x2 - 0.7 * x1 + x1 * x3,
        3.5 * x1 - 0.7 * x2 + x2 * x3,
        0.95 * x3 + 0.65 + 0.1 * x1^3 * x3 - x3^3 / 3 - 0.25 * x1^2 * x3 - x1^2 - 0.25 * x2^2 * x3 - x2^2
      ))
    },
    ic_ranges = list(c(-2, 2), c(-2, 2), c(-1, 2)),
    # ground-truth support in the R naming of build_design_matrix ("x1^2x3")
    truth = list(
      `1` = list(intercept = FALSE, terms = c("x1", "x2", "x1x3")),
      `2` = list(intercept = FALSE, terms = c("x1", "x2", "x2x3")),
      `3` = list(intercept = TRUE,  terms = c("x3", "x1^2", "x2^2", "x1^2x3", "x2^2x3", "x1^3x3", "x3^3"),
                 beta = c(0.95, -1, -1, -0.25, -0.25, 0.1, -1 / 3), beta0 = 0.65)
    )
  ),
  rossler = list(
    # from Experiments/bayesian-alasso-ro/rossler/rossler_snr.R
    coeff = list(list(-1, -1), list(1, 0.2), list(0.2, 1, -5.7)),
    names = list(list("x2", "x3"), list("x1", "x2"), list("", "x1x3", "x3")),
    rhs = function(t, x, p) {
      x1 <- x[1]; x2 <- x[2]; x3 <- x[3]
      list(c(-x2 - x3, x1 + 0.2 * x2, 0.2 + x1 * x3 - 5.7 * x3))
    },
    ic_ranges = list(c(-10, 10), c(-10, 10), c(0, 20)),
    truth = list(
      `1` = list(intercept = FALSE, terms = c("x2", "x3"), beta = c(-1, -1), beta0 = 0),
      `2` = list(intercept = FALSE, terms = c("x1", "x2")),
      `3` = list(intercept = TRUE,  terms = c("x3", "x1x3"))
    )
  ),
  dadras = list(
    # from Experiments/bayesian-alasso-ro/dadras/dadras_n.R
    coeff = list(list(1, -3, 2.7), list(1.7, -1, 1), list(2, -9)),
    names = list(list("x2", "x1", "x2x3"), list("x2", "x1x3", "x3"), list("x1x2", "x3")),
    rhs = function(t, x, p) {
      x1 <- x[1]; x2 <- x[2]; x3 <- x[3]
      list(c(x2 - 3 * x1 + 2.7 * x2 * x3, 1.7 * x2 - x1 * x3 + x3, 2 * x1 * x2 - 9 * x3))
    },
    # dadras_n.R draws all three coordinates in one runif(3, -4, 4) call per
    # trial (see dadras_ic_pool in dadras_psis_loo.R), unlike the per-coordinate
    # loop of initial_conditions() below.
    ic_ranges = list(c(-4, 4), c(-4, 4), c(-4, 4)),
    truth = list(
      `1` = list(intercept = FALSE, terms = c("x1", "x2", "x2x3"), beta = c(-3, 1, 2.7), beta0 = 0),
      `2` = list(intercept = FALSE, terms = c("x2", "x3", "x1x3"), beta = c(1.7, 1, -1), beta0 = 0),
      `3` = list(intercept = FALSE, terms = c("x3", "x1x2"), beta = c(-9, 2), beta0 = 0)
    )
  )
)

rhs_matrix <- function(system, x) {
  f <- SYSTEMS[[system]]$rhs
  t(apply(x, 1, function(row) f(0, row, NULL)[[1]]))
}

# Same seeded pool of initial conditions as the benchmark scripts
# (system_n.R: set.seed(seed); for i in 1:100 { x <- runif(1, ...); y <- ...; z <- ... }).
initial_conditions <- function(system, num = 100, seed = 100) {
  r <- SYSTEMS[[system]]$ic_ranges
  set.seed(seed)
  out <- matrix(NA_real_, num, 3)
  for (i in seq_len(num)) {
    out[i, 1] <- runif(1, r[[1]][1], r[[1]][2])
    out[i, 2] <- runif(1, r[[2]][1], r[[2]][2])
    out[i, 3] <- runif(1, r[[3]][1], r[[3]][2])
  }
  out
}

# Data generation through the benchmark code path: ode_auto.py's
# generate_noisy_dynamical_systems_pyversion (scipy.integrate.odeint, additive
# Gaussian noise with sd = 10^(-SNR/20) * std of the clean coordinate).  Python
# is initialised lazily inside each process so that mclapply workers (forked)
# each own their interpreter.  The clean trajectory is the same call with
# snr = Inf (snr_volt = 0, no noise added).  NumPy is seeded before the noisy
# call so runs are reproducible; the generator itself is unchanged.
.py_ready <- FALSE
py_generator <- function() {
  if (!.py_ready) {
    reticulate::use_python(PYTHON_BIN, required = TRUE)
    reticulate::source_python(file.path(ROOT, "DataGeneration", "ode_auto.py"), envir = globalenv())
    assign(".np", reticulate::import("numpy"), envir = globalenv())
    assign(".py_ready", TRUE, envir = globalenv())
  }
  get("generate_noisy_dynamical_systems_pyversion", envir = globalenv())
}

simulate <- function(system, n, snr, ic, noise_seed) {
  gen <- py_generator()
  s <- SYSTEMS[[system]]
  ic <- as.list(as.numeric(ic))
  x_clean <- gen(variable_coeff = s$coeff, variable_names = s$names, n = as.integer(n), dt = DT,
                 init_conditions = ic, snr = Inf)
  x_clean <- as.matrix(x_clean)
  if (is.infinite(snr)) return(list(clean = x_clean, noisy = x_clean))
  .np$random$seed(as.integer(noise_seed))
  x_noisy <- as.matrix(gen(variable_coeff = s$coeff, variable_names = s$names, n = as.integer(n), dt = DT,
                           init_conditions = ic, snr = snr))
  list(clean = x_clean, noisy = x_noisy)
}

# The manuscript's library builder (degree 5, SG polynomial order 4).
library_dm <- function(x_noisy) {
  build_design_matrix(x_t = x_noisy, dt = DT, sg_poly_order = SG_POLY_ORDER,
                      library_degree = LIBRARY_DEGREE, library_type = "poly")
}

# Polynomial library evaluated on arbitrary state samples (no SG): all
# monomials x1^a x2^b x3^c with 1 <= a+b+c <= degree, named exactly as
# build_design_matrix names them ("x1^2x3": tokens sorted, power 1 implicit).
# Column order is by total degree; only the name set matters for VIF / kappa.
poly_features <- function(x, degree = LIBRARY_DEGREE) {
  x <- as.matrix(x)
  ex <- expand.grid(a = 0:degree, b = 0:degree, c = 0:degree)
  ex <- ex[rowSums(ex) >= 1 & rowSums(ex) <= degree, ]
  ex <- ex[order(rowSums(ex)), ]
  nm <- apply(ex, 1, function(p) {
    toks <- c(if (p[1] > 0) if (p[1] == 1) "x1" else paste0("x1^", p[1]),
              if (p[2] > 0) if (p[2] == 1) "x2" else paste0("x2^", p[2]),
              if (p[3] > 0) if (p[3] == 1) "x3" else paste0("x3^", p[3]))
    paste(sort(toks), collapse = "")
  })
  theta <- sapply(seq_len(nrow(ex)), function(i) x[, 1]^ex$a[i] * x[, 2]^ex$b[i] * x[, 3]^ex$c[i])
  colnames(theta) <- nm
  list(theta = theta, names = unname(nm), degree = rowSums(ex))
}

monomial_degree <- function(name) {
  toks <- regmatches(name, gregexpr("x[0-9](\\^[0-9])?", name))[[1]]
  sum(sapply(toks, function(t) if (grepl("\\^", t)) as.integer(sub(".*\\^", "", t)) else 1L))
}

# R naming "x1^2x3" -> the sklearn-style naming "x1^2 x3" used by the Python
# scripts and by make_figures.py, so both result sets share one schema.
to_py_name <- function(name) {
  sapply(name, function(nm) paste(regmatches(nm, gregexpr("x[0-9](\\^[0-9])?", nm))[[1]], collapse = " "),
         USE.NAMES = FALSE)
}

# ---------------------------------------------------- screening (no Stan)
# Reproduces double_regression(): pass 1 with ridge weights sets the degree
# cut, pass 2 with OLS weights selects the support.
screen <- function(dm, eq) {
  theta <- dm$sorted_theta
  target <- dm$xdot_filtered[, eq]
  data <- cbind.data.frame(target = target, theta)
  init <- alasso(data, weights_method = "ridge", ols_ps = TRUE)
  init[is.na(init)] <- 0
  nz_max <- max(which(init != 0))
  cut <- sum(dm$monomial_orders <= dm$monomial_orders[nz_max])
  if (is.na(cut) || cut == length(dm$monomial_orders)) {
    post <- data
  } else {
    post <- cbind.data.frame(target = target, data[-1][, 1:cut])
  }
  final <- alasso(post, weights_method = "ols", ols_ps = TRUE)
  final[is.na(final)] <- 0
  nz <- final != 0
  trimmed_names <- colnames(post)[-1]
  list(trimmed_names = trimmed_names,
       trimmed_theta = as.matrix(post[, -1, drop = FALSE]),
       selected_intercept = unname(nz[1]),
       selected_terms = trimmed_names[nz[-1]])
}

# ------------------------------------------------ collinearity diagnostics
# Centered VIF_j = 1 / (1 - R_j^2), R_j^2 from regressing standardised column j
# on all others (with intercept).  Computed by least squares per column so it
# is finite and >= 1 even when p ~ n.
vif <- function(X, names) {
  X <- as.matrix(X); sdv <- apply(X, 2, sd); keep <- sdv > 0
  out <- setNames(rep(NA_real_, length(names)), names)
  if (sum(keep) < 2) { out[keep] <- 1; return(out) }
  Z <- scale(X[, keep, drop = FALSE])
  n <- nrow(Z); v <- numeric(ncol(Z))
  for (j in seq_len(ncol(Z))) {
    others <- Z[, -j, drop = FALSE]
    fit <- lm.fit(others, Z[, j])
    r2 <- 1 - sum(fit$residuals^2) / (n - 1)
    v[j] <- 1 / max(1 - r2, 1e-16)
  }
  out[names[keep]] <- v
  out
}

# Belsley-Kuh-Welsch: unit-length scaling (incl. intercept), condition
# indices and variance-decomposition proportions (k x j).
belsley <- function(X, names, intercept = TRUE) {
  X <- as.matrix(X); cols <- names
  if (intercept) { X <- cbind(1, X); cols <- c("Intercept", cols) }
  Xs <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  sv <- svd(Xs)
  eta <- max(sv$d) / sv$d
  phi_raw <- (sv$v^2) / matrix(sv$d^2, nrow(sv$v), ncol(sv$v), byrow = TRUE)  # j x k
  phi <- t(phi_raw / rowSums(phi_raw))                                          # k x j
  colnames(phi) <- cols
  list(condition_number = max(eta), condition_indices = eta, proportions = phi)
}

near_dependencies <- function(bel, eta_min = 30, phi_min = 0.5) {
  groups <- list()
  for (k in seq_along(bel$condition_indices)) {
    if (bel$condition_indices[k] < eta_min) next
    members <- colnames(bel$proportions)[bel$proportions[k, ] > phi_min]
    if (length(members) >= 2) groups[[length(groups) + 1]] <- list(condition_index = bel$condition_indices[k], members = members)
  }
  groups
}

theta_columns <- function(theta, names, terms) as.matrix(theta)[, match(terms, names), drop = FALSE]

# ------------------------------------------ residual-structure diagnostics
residual_stats <- function(y, X, intercept, n_bins = 10) {
  X <- as.matrix(X)
  fit <- if (intercept) lm(y ~ X) else lm(y ~ 0 + X)
  e <- residuals(fit); f <- fitted(fit)
  bp <- lmtest::bptest(fit, varformula = ~ f, studentize = TRUE)           # aux. regressor = fitted value
  wh <- tryCatch(lmtest::bptest(fit, varformula = ~ f + I(f^2), studentize = TRUE),
                 error = function(e) list(statistic = NA, p.value = NA))
  ord <- order(f)
  bins <- split(e[ord], cut(seq_along(ord), n_bins, labels = FALSE))
  v <- sapply(bins, var)
  list(sigma = sd(e), bp_lm = unname(bp$statistic), bp_p = unname(bp$p.value), bp_r2 = unname(bp$statistic) / length(e),
       white_lm = unname(wh$statistic), white_p = unname(wh$p.value),
       var_ratio = max(v) / min(v), rho1 = cor(e[-length(e)], e[-1]),
       durbin_watson = sum(diff(e)^2) / sum(e^2),
       residuals = e, fitted = f, coef = coef(fit))
}

# ------------------------------------------------ Bayesian stage (package code path)
# Mirrors bayesian_alasso_ro(): rstanarm::stan_glm with its default priors on
# the given support, 90 % posterior interval, keep a term iff the interval
# excludes zero.  Returns the kept terms, the posterior-mean fit and residuals.
bayes_stage <- function(y, X, intercept, ci_level = 0.9, cores = 1) {
  X <- as.data.frame(X)
  dat <- cbind.data.frame(target = y, X)
  fml <- if (intercept) target ~ . else target ~ 0 + .
  glm_stan <- stan_glm(fml, family = gaussian(), data = dat, refresh = 0, cores = cores)
  CI <- posterior_interval(glm_stan, prob = ci_level)
  rownames(CI) <- gsub("`", "", rownames(CI))
  ci <- t(CI)
  cf <- coef(glm_stan); names(cf) <- gsub("`", "", names(cf))
  keep <- sapply(seq_along(cf), function(i) {
    nm <- names(cf)[i]
    ci[1, nm] <= cf[i] && ci[2, nm] >= cf[i] &&
      !((ci[1, nm] <= 0 && ci[2, nm] >= 0) || (ci[1, nm] >= 0 && ci[2, nm] <= 0))
  })
  kept <- names(cf)[keep]
  f <- as.numeric(fitted(glm_stan))            # posterior mean of the linear predictor
  list(kept_intercept = "(Intercept)" %in% kept, kept_terms = setdiff(kept, "(Intercept)"),
       fitted = f, residuals = y - f, coef = cf, ci = ci, sigma = sd(y - f), fit = glm_stan)
}

# Residual statistics from given residuals / fitted values (no refit)
residual_stats_from <- function(e, f, n_bins = 10) {
  aux <- lm(I(e^2) ~ f)
  bp_lm <- length(e) * summary(aux)$r.squared
  bp_p <- pchisq(bp_lm, df = 1, lower.tail = FALSE)
  aux2 <- lm(I(e^2) ~ f + I(f^2))
  white_lm <- length(e) * summary(aux2)$r.squared
  white_p <- pchisq(white_lm, df = 2, lower.tail = FALSE)
  ord <- order(f)
  bins <- split(e[ord], cut(seq_along(ord), n_bins, labels = FALSE))
  v <- sapply(bins, var)
  list(sigma = sd(e), bp_lm = bp_lm, bp_p = bp_p, bp_r2 = bp_lm / length(e),
       white_lm = white_lm, white_p = white_p,
       var_ratio = max(v) / min(v), rho1 = cor(e[-length(e)], e[-1]),
       durbin_watson = sum(diff(e)^2) / sum(e^2), residuals = e, fitted = f)
}

sg_params <- function(x_col) {
  comb <- sg_optimal_combination(x_col, DT, polyorder = SG_POLY_ORDER)[[2]]
  c(order = as.numeric(comb[[1]][1]), window = as.numeric(comb[[2]][1]))
}
sg_derivative <- function(x_col, order, window) sgolayfilt(x_col, p = order, n = window, m = 1, ts = DT)

log_grid <- function(start = 2, stop = 5, step = 0.2) {
  k <- round(seq(start, stop, by = step), 1)
  data.frame(log10n = k, n = round(10^k))
}
snr_label <- function(snr) if (is.infinite(snr)) "inf" else as.character(as.integer(snr))
