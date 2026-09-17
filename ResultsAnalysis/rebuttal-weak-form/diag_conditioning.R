# Diagnostic: conditioning of the candidate library under the weak transform.
# Aizawa, n = 5000, 49 dB, 5 trials: Belsley condition number (columns scaled
# to unit length, intercept included) and VIF of the seven true x3-dot terms,
# for the SG design (library on SG-smoothed states) and the weak design
# (PySINDy WeakPDELibrary, K_FRAC 0.25 and 0.04), on the full degree-5 library
# and on the degree-<=4 library of Fig. 3b.
# Usage: ARGOS_ROOT=<repo> ARGOS_PYTHON=<python> Rscript diag_conditioning.R
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "..", "..", "Experiments", "bayesian-alasso-ro-weak", "weak_form_common.R"))
TRUE3 <- c("x3","x1^2","x2^2","x1^2x3","x2^2x3","x1^3x3","x3^3")
kappa <- function(X) { X <- cbind(1, as.matrix(X)); Xs <- X / sqrt(colSums(X^2)); s <- svd(Xs, nu = 0, nv = 0)$d; max(s) / min(s) }
vif_true <- function(X, names) { X <- as.matrix(X); Z <- scale(X); C <- crossprod(Z) / nrow(Z); v <- diag(solve(C)); names(v) <- names; v[TRUE3] }
init <- SYSTEMS$aizawa$init(100, 100)
rows <- list()
for (j in 1:5) {
  xn <- simulate_system("aizawa", 5000, 0.01, init[[j]], 49)
  designs <- list(SG = build_design_matrix(x_t = xn, dt = 0.01, sg_poly_order = 4, library_degree = 5, library_type = "poly"))
  for (kf in c(0.25, 0.04)) { K_FRAC <<- kf; designs[[sprintf("weak K=%gn", kf)]] <- build_weak_design(xn, 0.01, 5, 1000 + j) }
  for (nm in names(designs)) { d <- designs[[nm]]; th <- as.data.frame(d$sorted_theta); deg <- d$monomial_orders; cn <- colnames(th)
    th4 <- th[, deg <= 4, drop = FALSE]
    v <- vif_true(th4, colnames(th4))
    rows[[length(rows) + 1]] <- data.frame(trial = j, design = nm, rows = nrow(th), kappa_full = kappa(th), kappa_deg4 = kappa(th4), kappa_true = kappa(th[, TRUE3]),
      vif_true_max = max(v), vif_true_med = median(v)) }
}
d <- do.call(rbind, rows)
s <- aggregate(cbind(rows, kappa_full, kappa_deg4, kappa_true, vif_true_max, vif_true_med) ~ design, d, median)
print(s, row.names = FALSE, digits = 3)
write.csv(d, file.path(here, "results", "diag_conditioning.csv"), row.names = FALSE)
