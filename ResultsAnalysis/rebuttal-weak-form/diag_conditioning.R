# Diagnostic: conditioning of the candidate library under the weak transform,
# at the paper's conventions (rewritten 2026-11-08; the 5-trial, degree-4
# version it replaces is kept in the git history).
# Aizawa, n = 5000, 49 dB, 20 trials: the same 20 trajectories and noise seeds
# as Supplementary Table S1 (table_s1_nested_vif_49db.R, exp1's ICS[1:20] and
# noise-seed scheme), the full degree-5 candidate library, and the Belsley
# condition number of ../rebuttal-diagnostics/common.R (columns scaled to unit
# length, intercept column included), for
#   SG    the method's design matrix (library on the SG-smoothed states; must
#         reproduce Table S1's "full" kappa trial by trial, checked below),
#   weak  the weak-form design of the 2 x 2 comparison on the same noisy
#         record (PySINDy WeakPDELibrary, K = n/4 subdomains, H = 0.25, p = 8,
#         Experiments/bayesian-alasso-ro-weak/weak_form_common.R), and
#   weak K=0.04n for reference (about 50 % overlap).
# Per design: kappa of the full library, of its degree <= 4 columns and of the
# seven true x3-dot terms alone, and the VIF of the true terms on the full
# library (max and median).  Output: results/diag_conditioning.csv (one row
# per trial x design) and medians with interquartile ranges on stdout.
# Usage: Rscript diag_conditioning.R   (ARGOS_ROOT and ARGOS_PYTHON default to the repository and the env_dgx_spark python; about a minute)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
Sys.setenv(ARGOS_PYTHON = Sys.getenv("ARGOS_PYTHON", "/home/yuzhengz/anaconda3/envs/env_dgx_spark/bin/python"))
if (!nzchar(Sys.getenv("ARGOS_ROOT"))) Sys.setenv(ARGOS_ROOT = normalizePath(file.path(here, "..", "..")))
source(file.path(here, "..", "..", "Experiments", "bayesian-alasso-ro-weak", "weak_form_common.R"))  # build_weak_design, K_FRAC, ps, np
source(file.path(here, "..", "rebuttal-diagnostics", "common.R"))   # sourced second: its SYSTEMS, simulate, library_dm, belsley take precedence
TRIALS <- 20L; N <- 5000L; SNR <- 49
TRUE3 <- c("x3", "x1^2", "x2^2", "x1^2x3", "x2^2x3", "x1^3x3", "x3^3")
ICS <- initial_conditions("aizawa")
noise_seed <- function(trial) 10000 * trial + round(log10(N) * 10) + 5000 + SNR   # exp1_collinearity.R, snr sweep
vif_of <- function(X) { Z <- scale(as.matrix(X)); diag(solve(crossprod(Z) / (nrow(Z) - 1))) }
rows <- list()
for (trial in seq_len(TRIALS) - 1) {
  s <- simulate("aizawa", N, SNR, ICS[trial + 1, ], noise_seed(trial))
  designs <- list(SG = library_dm(s$noisy))
  for (kf in c(0.25, 0.04)) { K_FRAC <<- kf; designs[[sprintf("weak K=%gn", kf)]] <- build_weak_design(s$noisy, DT, LIBRARY_DEGREE, 1000 + trial) }
  for (nm in names(designs)) {
    d <- designs[[nm]]; th <- as.data.frame(d$sorted_theta); deg <- d$monomial_orders; cn <- colnames(th)
    stopifnot(all(TRUE3 %in% cn), ncol(th) == 55)
    th4 <- th[, deg <= 4, drop = FALSE]; v <- vif_of(th); names(v) <- cn
    rows[[length(rows) + 1]] <- data.frame(trial = trial, design = nm, rows = nrow(th), p = ncol(th),
      kappa_full = belsley(th, cn)$condition_number, kappa_deg4 = belsley(th4, colnames(th4))$condition_number,
      kappa_true = belsley(th[, TRUE3], TRUE3)$condition_number, vif_true_max = max(v[TRUE3]), vif_true_med = median(v[TRUE3]))
  }
  cat(sprintf("trial %2d done\n", trial))
}
d <- do.call(rbind, rows)
write.csv(d, file.path(here, "results", "diag_conditioning.csv"), row.names = FALSE)
# check: the SG design must reproduce Table S1's full-library kappa trial by trial
ts1 <- read.csv(file.path(here, "..", "rebuttal-diagnostics", "results", "table_s1_kappa_49db.csv"))
ts1 <- ts1[ts1$design == "full", ]; sg <- d[d$design == "SG", ]
m <- merge(sg[, c("trial", "kappa_full")], ts1[, c("trial", "kappa")], by = "trial")
cat(sprintf("\nTable S1 check: %d trials matched, max relative difference in kappa_full = %.2e\n", nrow(m), max(abs(m$kappa_full / m$kappa - 1))))
q <- function(x) sprintf("%.3g [%.3g, %.3g]", median(x), quantile(x, .25), quantile(x, .75))
cat("\nmedian [q25, q75] over", TRIALS, "trials\n")
for (nm in unique(d$design)) { x <- d[d$design == nm, ]
  cat(sprintf("%-14s rows %5d  kappa_full %-22s kappa_deg4 %-22s kappa_true %-20s VIF_true max %-24s med %s\n", nm, median(x$rows), q(x$kappa_full), q(x$kappa_deg4), q(x$kappa_true), q(x$vif_true_max), q(x$vif_true_med))) }
