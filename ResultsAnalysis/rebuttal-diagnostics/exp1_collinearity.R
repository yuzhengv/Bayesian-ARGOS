# Experiment 1 (R) -- where does the Aizawa x3-dot collinearity come from?
#
# Same design as the Python exp1 in pyargos/results-diagnostics/rebuttal-
# diagnostics, but every step uses the manuscript's R pipeline
# (build_design_matrix, the two alasso passes).  Outputs share the Python CSV
# schema (term names converted to "x1^2 x3" style) so make_figures.py can be
# pointed at this folder.
#
# Part A: per-term VIF and Belsley condition number on four nested designs
#   full / trimmed (Fig. 3b's design) / selected (pass-2 support) / truth
#   over the n grid (49 dB) and the SNR grid (n = 5000).
# Part B: same monomials on the SG trajectory, the clean trajectory, uniform
#   box samples, 50 short transients (clean and noisy+SG); clean trajectory
#   over the n grid; Belsley near-dependency groups for one trial.
#
# Usage: Rscript exp1_collinearity.R [trials_snr=20] [trials_n=10]

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
args <- commandArgs(trailingOnly = TRUE)
TRIALS <- if (length(args) >= 1) as.integer(args[1]) else 20L
TRIALS_N <- if (length(args) >= 2) as.integer(args[2]) else 10L
OUT <- file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "results")
dir.create(OUT, showWarnings = FALSE)

EQ <- 3; SYSTEM <- "aizawa"
N_GRID <- log_grid(2, 5, 0.2)
SNR_GRID <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, Inf)
SNR_FIXED <- 49; N_FIXED <- 5000; K_SEGMENTS <- 50
TRUTH <- SYSTEMS[[SYSTEM]]$truth[[as.character(EQ)]]$terms
ICS <- initial_conditions(SYSTEM)

designs_for <- function(dm, scr) {
  names_all <- colnames(dm$sorted_theta)
  out <- list(full = list(dm$sorted_theta, names_all),
              trimmed = list(scr$trimmed_theta, scr$trimmed_names),
              truth = list(theta_columns(dm$sorted_theta, names_all, TRUTH), TRUTH))
  if (length(scr$selected_terms) >= 2)
    out$selected <- list(theta_columns(dm$sorted_theta, names_all, scr$selected_terms), scr$selected_terms)
  out
}

one_condition <- function(trial, k, n, snr, sweep) {
  t0 <- Sys.time()
  noise_seed <- 10000 * trial + round(k * 10) + if (sweep == "n") 0 else 5000 + if (is.infinite(snr)) 0 else snr
  s <- simulate(SYSTEM, n, snr, ICS[trial + 1, ], noise_seed)
  dm <- library_dm(s$noisy)
  scr <- screen(dm, EQ)
  vif_rows <- list(); kappa_rows <- list()
  for (design in names(designs_for(dm, scr))) {
    d <- designs_for(dm, scr)[[design]]; X <- d[[1]]; nm <- d[[2]]
    v <- vif(X, nm); bel <- belsley(X, nm)
    vif_rows[[design]] <- data.frame(trial = trial, sweep = sweep, log10n = k, n = n, snr = snr, design = design,
                                     term = to_py_name(nm), degree = sapply(nm, monomial_degree),
                                     is_truth = nm %in% TRUTH, vif = unname(v))
    kappa_rows[[design]] <- data.frame(trial = trial, sweep = sweep, log10n = k, n = n, snr = snr, design = design,
                                       p = length(nm), kappa = bel$condition_number,
                                       vif_max = max(v, na.rm = TRUE), vif_median = median(v, na.rm = TRUE),
                                       n_vif_gt10 = sum(v > 10, na.rm = TRUE))
  }
  meta <- data.frame(trial = trial, sweep = sweep, log10n = k, n = n, snr = snr,
                     trimmed_p = length(scr$trimmed_names),
                     trimmed_max_degree = max(sapply(scr$trimmed_names, monomial_degree)),
                     selected = paste(to_py_name(scr$selected_terms), collapse = ";"),
                     selected_intercept = scr$selected_intercept,
                     screen_exact_truth = setequal(scr$selected_terms, TRUTH) && scr$selected_intercept,
                     truth_in_selected = all(TRUTH %in% scr$selected_terms),
                     seconds = as.numeric(Sys.time() - t0, units = "secs"))
  list(vif = do.call(rbind, vif_rows), kappa = do.call(rbind, kappa_rows), meta = meta)
}

geometry <- function(trial) {
  set.seed(5000 + trial)
  rows <- list()
  record <- function(label, theta, nm) {
    for (design in c("full", "truth")) {
      X <- if (design == "full") theta else theta_columns(theta, nm, TRUTH)
      nmd <- if (design == "full") nm else TRUTH
      v <- vif(X, nmd); bel <- belsley(X, nmd)
      rows[[length(rows) + 1]] <<- data.frame(trial = trial, sampling = label, design = design, p = length(nmd),
                                              kappa = bel$condition_number, vif_max = max(v, na.rm = TRUE),
                                              vif_median = median(v, na.rm = TRUE))
    }
  }
  s <- simulate(SYSTEM, N_FIXED, SNR_FIXED, ICS[trial + 1, ], 777 + trial)
  dm <- library_dm(s$noisy)
  record("trajectory_sg", dm$sorted_theta, colnames(dm$sorted_theta))
  pf <- poly_features(s$clean); record("trajectory_clean", pf$theta, pf$names)
  lo <- apply(s$clean, 2, min); hi <- apply(s$clean, 2, max)
  box <- sapply(1:3, function(j) runif(N_FIXED, lo[j], hi[j]))
  pf <- poly_features(box); record("box_uniform", pf$theta, pf$names)
  r <- SYSTEMS[[SYSTEM]]$ic_ranges; seg_len <- N_FIXED %/% K_SEGMENTS
  segs_clean <- list(); segs_sg <- list()
  for (j in seq_len(K_SEGMENTS)) {
    ic_j <- c(runif(1, r[[1]][1], r[[1]][2]), runif(1, r[[2]][1], r[[2]][2]), runif(1, r[[3]][1], r[[3]][2]))
    sj <- simulate(SYSTEM, seg_len, SNR_FIXED, ic_j, 9000 + 100 * trial + j)
    segs_clean[[j]] <- sj$clean
    segs_sg[[j]] <- library_dm(sj$noisy)$sorted_theta
  }
  pf <- poly_features(do.call(rbind, segs_clean)); record("multi_transient", pf$theta, pf$names)
  record("multi_transient_sg", do.call(rbind, segs_sg), pf$names)
  do.call(rbind, rows)
}

clean_vs_n <- function(trial) {
  s <- simulate(SYSTEM, max(N_GRID$n), Inf, ICS[trial + 1, ], 0)
  rows <- list()
  for (i in seq_len(nrow(N_GRID))) {
    pf <- poly_features(s$clean[1:N_GRID$n[i], ])
    for (design in c("full", "truth")) {
      X <- if (design == "full") pf$theta else theta_columns(pf$theta, pf$names, TRUTH)
      nmd <- if (design == "full") pf$names else TRUTH
      v <- vif(X, nmd)
      rows[[length(rows) + 1]] <- data.frame(trial = trial, log10n = N_GRID$log10n[i], n = N_GRID$n[i], design = design,
                                             kappa = belsley(X, nmd)$condition_number,
                                             vif_max = max(v, na.rm = TRUE), vif_median = median(v, na.rm = TRUE))
    }
  }
  do.call(rbind, rows)
}

# ------------------------------------------------------------------ Part A
tasks <- list()
for (trial in seq_len(max(TRIALS, TRIALS_N)) - 1) {
  if (trial < TRIALS_N) for (i in seq_len(nrow(N_GRID))) tasks[[length(tasks) + 1]] <- list(trial, N_GRID$log10n[i], N_GRID$n[i], SNR_FIXED, "n")
  if (trial < TRIALS) for (snr in SNR_GRID) tasks[[length(tasks) + 1]] <- list(trial, log10(N_FIXED), N_FIXED, snr, "snr")
}
tasks <- tasks[order(-sapply(tasks, function(t) t[[3]]))]
cat("Part A:", length(tasks), "conditions on", NCORES, "cores\n")
t0 <- Sys.time()
res <- mclapply(tasks, function(t) one_condition(t[[1]], t[[2]], t[[3]], t[[4]], t[[5]]), mc.cores = NCORES)
ok <- !sapply(res, inherits, "try-error")
write.csv(do.call(rbind, lapply(res[ok], `[[`, "vif")), file.path(OUT, "exp1_vif_long.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(res[ok], `[[`, "kappa")), file.path(OUT, "exp1_kappa.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(res[ok], `[[`, "meta")), file.path(OUT, "exp1_screening_meta.csv"), row.names = FALSE)
cat("Part A done in", round(as.numeric(Sys.time() - t0, units = "secs")), "s;", sum(!ok), "failures\n")

# ------------------------------------------------------------------ Part B
t0 <- Sys.time()
geo <- mclapply(seq_len(TRIALS) - 1, geometry, mc.cores = NCORES)
write.csv(do.call(rbind, geo), file.path(OUT, "exp1_geometry.csv"), row.names = FALSE)
cvn <- mclapply(seq_len(TRIALS) - 1, clean_vs_n, mc.cores = NCORES)
write.csv(do.call(rbind, cvn), file.path(OUT, "exp1_clean_vs_n.csv"), row.names = FALSE)

s <- simulate(SYSTEM, N_FIXED, SNR_FIXED, ICS[1, ], 10000)
dm <- library_dm(s$noisy); scr <- screen(dm, EQ)
report <- list()
for (design in names(designs_for(dm, scr))) {
  d <- designs_for(dm, scr)[[design]]; X <- d[[1]]; nm <- d[[2]]
  bel <- belsley(X, nm)
  nd <- lapply(near_dependencies(bel), function(g) list(condition_index = g$condition_index, members = to_py_name(g$members)))
  report[[design]] <- list(p = length(nm), condition_number = bel$condition_number, near_dependencies = nd,
                           vif = as.list(setNames(vif(X, nm), to_py_name(nm))))
}
writeLines(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(OUT, "exp1_near_dependencies.json"))
cat("Part B done in", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
