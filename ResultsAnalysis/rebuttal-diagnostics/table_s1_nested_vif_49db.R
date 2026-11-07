# Supplementary Table S1 at the benchmark convention (n = 5000, SNR = 49 dB).
#
# exp1_collinearity.R computes the nested-design VIFs and condition numbers on
# its SNR grid (1, 5, 10, ..., 60 dB, Inf) at n = 5000, so Table S1 was first
# drawn from the 50 dB point.  This script runs the single condition
# n = 5000, SNR = 49 dB for the same 20 trajectories (initial conditions
# ICS[1:20, ], exp1's noise-seed scheme), reusing exp1's one_condition() and
# designs_for() unchanged: the function definitions of exp1_collinearity.R are
# evaluated without running its sweeps.
#
# Usage: Rscript table_s1_nested_vif_49db.R [trials=20]
# Outputs: results/table_s1_vif_49db.csv, results/table_s1_kappa_49db.csv,
#          results/table_s1_meta_49db.csv, tables/table_s1_49db.md (LaTeX rows).

HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(HERE, "common.R"))
args <- commandArgs(trailingOnly = TRUE)
TRIALS <- if (length(args) >= 1) as.integer(args[1]) else 20L
OUT <- file.path(HERE, "results"); TAB <- file.path(HERE, "tables")
dir.create(OUT, showWarnings = FALSE); dir.create(TAB, showWarnings = FALSE)

# exp1's constants (identical values) and its function definitions only
EQ <- 3; SYSTEM <- "aizawa"
TRUTH <- SYSTEMS[[SYSTEM]]$truth[[as.character(EQ)]]$terms
ICS <- initial_conditions(SYSTEM)
for (e in parse(file.path(HERE, "exp1_collinearity.R"))) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e)
}
stopifnot(exists("one_condition"), exists("designs_for"))

N_POINT <- 5000; SNR_POINT <- 49
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(TRIALS) - 1, function(trial) one_condition(trial, log10(N_POINT), N_POINT, SNR_POINT, "snr"), mc.cores = NCORES)
ok <- !sapply(res, inherits, "try-error")
vif <- do.call(rbind, lapply(res[ok], `[[`, "vif")); kap <- do.call(rbind, lapply(res[ok], `[[`, "kappa")); meta <- do.call(rbind, lapply(res[ok], `[[`, "meta"))
write.csv(vif, file.path(OUT, "table_s1_vif_49db.csv"), row.names = FALSE)
write.csv(kap, file.path(OUT, "table_s1_kappa_49db.csv"), row.names = FALSE)
write.csv(meta, file.path(OUT, "table_s1_meta_49db.csv"), row.names = FALSE)

fmt <- function(x) {
  if (is.na(x)) return("--")
  if (x >= 1e4) { e <- floor(log10(x)); return(sprintf("$%.2f\\times10^{%d}$", x / 10^e, e)) }
  sub("\\.$", "", formatC(signif(x, 3), digits = 3, format = "fg", flag = "#"))
}
cell <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return("--"); sprintf("%s [%s, %s]", fmt(median(v)), fmt(quantile(v, .25)), fmt(quantile(v, .75))) }
designs <- c(full = "Full candidate library", trimmed = "After the degree trim", selected = "Screened support", truth = "True terms only")
term_tex <- c("x3" = "$x_3$", "x1^2" = "$x_1^2$", "x2^2" = "$x_2^2$", "x1^2 x3" = "$x_1^2x_3$", "x2^2 x3" = "$x_2^2x_3$", "x1^3 x3" = "$x_1^3x_3$", "x3^3" = "$x_3^3$")
truth_py <- to_py_name(TRUTH)
lines <- c(sprintf("# Table S1 at n = %d, SNR = %d dB (%d trials, %s)", N_POINT, SNR_POINT, sum(ok), format(Sys.Date())),
           "median [q25, q75] over the trials in which the term is in the design; -- = never present", "")
for (tm in truth_py) {
  cells <- sapply(names(designs), function(d) cell(vif$vif[vif$design == d & vif$term == tm]))
  present <- sum(vif$design == "selected" & vif$term == tm)
  lines <- c(lines, sprintf("%-10s & %s \\\\   %% in screened support in %d of %d trials", term_tex[tm], paste(cells, collapse = " & "), present, sum(ok)))
}
lines <- c(lines, sprintf("$\\kappa$   & %s \\\\", paste(sapply(names(designs), function(d) cell(kap$kappa[kap$design == d])), collapse = " & ")))
sel_p <- kap$p[kap$design == "selected"]
lines <- c(lines, "", sprintf("screened support size: %d to %d terms, median %s (over %d trials with >= 2 selected terms); trimmed library: %s terms (median)",
                            min(sel_p), max(sel_p), format(median(sel_p)), length(sel_p), format(median(kap$p[kap$design == "trimmed"]))),
           sprintf("x3^3 on the screened support: %s", paste(fmt_vals <- sort(round(vif$vif[vif$design == "selected" & vif$term == "x3^3"], 1)), collapse = ", ")),
           sprintf("runtime: %.0f s on %d cores", as.numeric(Sys.time() - t0, units = "secs"), NCORES))
writeLines(lines, file.path(TAB, "table_s1_49db.md")); cat(lines, sep = "\n")
