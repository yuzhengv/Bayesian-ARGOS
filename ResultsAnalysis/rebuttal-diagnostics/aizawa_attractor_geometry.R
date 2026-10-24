# Why the near-dependent chains of the Aizawa x3-dot candidate library differ by
# a factor of x3 (manuscript Section 2.3, Supplementary Note 4).
# For the 20 initial conditions of the benchmark pool, one clean trajectory each
# (n = 5000, dt = 0.01, through ode_auto.py as in exp1_collinearity.R):
#   * range, mean, fraction positive and mean frequency rms(xdot)/(2 pi sd(x))
#     of each state;
#   * |corr(theta, x_k * theta)| over all monomials of degree 1 to 4, for k = 1, 2, 3
#     (how much multiplying a column by a state changes its shape);
#   * |corr| of the exchanged pairs and of the chain of powers of x3, on the
#     trajectory, on uniform samples of the initial-condition box
#     [-2, 2] x [-2, 2] x [-1, 2], and on uniform samples with the box's x3 range
#     replaced by the attractor's (same sample size).
# Writes results/aizawa_attractor_geometry.csv (one row per trial) and prints
# medians over the trials. About 20 s.
suppressPackageStartupMessages(source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R")))
OUT <- file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "results")
SYSTEM <- "aizawa"; ICS <- initial_conditions(SYSTEM); N <- 5000
rhs <- SYSTEMS[[SYSTEM]]$rhs
PAIRS <- list(c("x1^2x3", "x1^2x3^2"), c("x2^2x3", "x2^2x3^2"), c("x3", "x3^2"), c("x3^2", "x3^3"), c("x3^3", "x3^4"))
pair_cor <- function(th) sapply(PAIRS, function(p) abs(cor(th[, p[1]], th[, p[2]])))
rows <- list()
for (trial in 0:19) {
  x <- simulate(SYSTEM, N, Inf, ICS[trial + 1, ], 0)$clean
  xd <- t(apply(x, 1, function(r) rhs(0, r, NULL)[[1]]))
  th <- poly_features(x, 4)$theta
  shape <- sapply(1:3, function(k) median(sapply(seq_len(ncol(th)), function(j) abs(cor(th[, j], th[, j] * x[, k]))), na.rm = TRUE))
  set.seed(100 + trial)
  xb <- cbind(runif(N, -2, 2), runif(N, -2, 2), runif(N, -1, 2))
  xb3 <- cbind(runif(N, -2, 2), runif(N, -2, 2), runif(N, min(x[, 3]), max(x[, 3])))
  r <- c(trial = trial,
         setNames(as.vector(sapply(1:3, function(j) c(min(x[, j]), max(x[, j]), mean(x[, j]), sd(x[, j]), mean(x[, j] > 0), sqrt(mean(xd[, j]^2)) / (2 * pi * sd(x[, j]))))),
                  as.vector(outer(c("min", "max", "mean", "sd", "frac_pos", "mean_freq"), 1:3, function(a, j) paste0("x", j, "_", a)))),
         setNames(shape, paste0("median_shape_cor_times_x", 1:3)),
         setNames(pair_cor(th), paste0("traj_", sapply(PAIRS, paste, collapse = "_vs_"))),
         setNames(pair_cor(poly_features(xb, 4)$theta), paste0("box_", sapply(PAIRS, paste, collapse = "_vs_"))),
         setNames(pair_cor(poly_features(xb3, 4)$theta), paste0("box_x3range_", sapply(PAIRS, paste, collapse = "_vs_"))))
  rows[[length(rows) + 1]] <- as.data.frame(as.list(r), check.names = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "aizawa_attractor_geometry.csv"), row.names = FALSE)
med <- apply(res[, -1], 2, median)
cat("medians over 20 trials (n = 5000, clean):\n"); print(round(med, 3))
