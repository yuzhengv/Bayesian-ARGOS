# The three diagnostic insets of Fig. 1e (revised manuscript) from ONE
# successfully identified Aizawa trial, drawn in the house style
# (../plot-style/house_style.R; plotting standard block M1).
#
# Pipeline per try (exactly the calls of exp2_heteroscedasticity.R and
# dadras_psis_loo.R): benchmark data path (ode_auto.py through reticulate),
# build_design_matrix (Savitzky-Golay, degree-5 library), the two adaptive-
# lasso passes (screen), rstanarm::stan_glm on the selected support with the
# 90 % posterior-interval rule (bayes_stage), for all three equations.  Initial
# conditions are taken from the benchmark pool in order (rows 1, 2, ...), with
# the noise seed exp2 uses for that pool row (20000 + 100 * (row - 1) + SNR),
# until one trial recovers all three equations exactly (at most MAX_TRIES).
#
# For that trial, x3-dot: standardised residuals of the Stan fit against the
# fitted value (posterior mean); Pareto k-hat from the pointwise Gaussian
# log-likelihood of the posterior draws (loo::psis with relative_eff, chunked,
# as dadras_psis_loo.R); VIF of every candidate term of the design entering the
# second adaptive-lasso pass.
#
# Outputs: results/fig1_insets_trial.csv (the trial), fig1_insets_residuals.csv,
#   fig1_insets_khat.csv, fig1_insets_vif.csv; figures/fig1_inset_residuals.pdf,
#   fig1_inset_khat.pdf, fig1_inset_vif.pdf (final width 34 mm, aspect 0.95,
#   canvas 3.35 x 3.18 in at k = 2.5, cairo_pdf; no tag, no title: the drawing
#   carries them).  34 mm is the width the box of the drawing gives the inset
#   at the final figure width of 183 mm (compose_fig1_schematic.sh scales the
#   inset to 186 pt of the 1025 pt page), so that the inset text is 7 pt and
#   the tick labels 6.5 pt in the composed figure.
# Usage: Rscript fig1_insets.R [snr=57] [n=5000]     (one try is about 1 min)

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
suppressPackageStartupMessages(library(loo))
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "..", "plot-style", "house_style.R"))
source(file.path(here, "fig_helpers.R"))
args <- commandArgs(trailingOnly = TRUE)
SNR <- if (length(args) >= 1) as.numeric(args[1]) else 57
N <- if (length(args) >= 2) as.integer(args[2]) else 5000L
MAX_TRIES <- 20L; EQ <- 3L; SYSTEM <- "aizawa"
RES <- file.path(here, "results"); FIG <- file.path(here, "figures")
dir.create(RES, showWarnings = FALSE); dir.create(FIG, showWarnings = FALSE)
ICS <- initial_conditions(SYSTEM)

# Pointwise Pareto k-hat from the rstanarm Gaussian fit, chunked (dadras_psis_loo.R)
CHUNK <- 4000L
pareto_k_chunked <- function(fit, y, X, intercept) {
  draws <- as.matrix(fit); cn <- gsub("`", "", colnames(draws))
  sigma <- draws[, cn == "sigma"]
  beta_names <- c(if (intercept) "(Intercept)", colnames(X))
  beta <- draws[, match(beta_names, cn), drop = FALSE]
  Xd <- if (intercept) cbind(1, as.matrix(X)) else as.matrix(X)
  S <- nrow(draws); n_chain <- length(fit$stanfit@sim$samples)
  chain_id <- rep(seq_len(n_chain), each = S / n_chain)
  k <- numeric(length(y))
  for (start in seq(1, length(y), by = CHUNK)) {
    idx <- start:min(start + CHUNK - 1, length(y))
    mu <- beta %*% t(Xd[idx, , drop = FALSE])
    ll <- dnorm(matrix(y[idx], S, length(idx), byrow = TRUE), mu, sigma, log = TRUE)
    r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id, cores = 1)
    ps <- suppressWarnings(loo::psis(-ll, r_eff = r_eff, cores = 1))
    k[idx] <- loo::pareto_k_values(ps)
  }
  k
}

fit_equation <- function(dm, eq) {
  tr <- SYSTEMS[[SYSTEM]]$truth[[as.character(eq)]]
  nm <- colnames(dm$sorted_theta)
  scr <- screen(dm, eq); sel <- scr$selected_terms
  y <- dm$xdot_filtered[, eq]
  if (length(sel) == 0) return(list(eq = eq, ok = FALSE, scr = scr))
  Xs <- theta_columns(dm$sorted_theta, nm, sel); colnames(Xs) <- sel
  bs <- bayes_stage(y, Xs, scr$selected_intercept)
  ok <- setequal(bs$kept_terms, tr$terms) && bs$kept_intercept == tr$intercept
  list(eq = eq, ok = ok, scr = scr, bs = bs, y = y, Xs = Xs, kept = c(if (bs$kept_intercept) "(Intercept)", bs$kept_terms))
}

RNGkind("L'Ecuyer-CMRG"); set.seed(1)        # reproducible Stan seeds in the forked workers
t0 <- Sys.time(); chosen <- NULL; log_rows <- list()
for (row in seq_len(MAX_TRIES)) {
  seed <- 20000 + 100 * (row - 1) + SNR
  s <- simulate(SYSTEM, N, SNR, ICS[row, ], seed)
  dm <- library_dm(s$noisy)
  fits <- mclapply(1:3, function(eq) fit_equation(dm, eq), mc.cores = 3)     # one Stan fit per forked worker
  ok <- sapply(fits, `[[`, "ok")
  log_rows[[row]] <- data.frame(pool_row = row, noise_seed = seed, n = N, snr = SNR,
                                eq1_ok = ok[1], eq2_ok = ok[2], eq3_ok = ok[3],
                                eq3_final = paste(fits[[3]]$kept, collapse = "|"))
  cat(sprintf("pool row %d: equations recovered %s  (%.0f s)\n", row, paste(ok, collapse = " "), as.numeric(Sys.time() - t0, units = "secs")))
  if (all(ok)) { chosen <- list(row = row, seed = seed, s = s, dm = dm, fits = fits); break }
}
write.csv(do.call(rbind, log_rows), file.path(RES, "fig1_insets_trial.csv"), row.names = FALSE)
if (is.null(chosen)) stop("no trial recovered all three equations within ", MAX_TRIES, " tries")

f3 <- chosen$fits[[EQ]]; bs <- f3$bs; scr <- f3$scr
e <- bs$residuals; f <- bs$fitted; z <- e / sd(e)
tt <- (seq_along(f3$y) - 1) * DT
k <- pareto_k_chunked(bs$fit, f3$y, f3$Xs, scr$selected_intercept)
v <- vif(scr$trimmed_theta, scr$trimmed_names)
truth3 <- SYSTEMS[[SYSTEM]]$truth[[as.character(EQ)]]$terms
vif_df <- data.frame(index = seq_along(v), term = names(v), vif = unname(v), is_truth = names(v) %in% truth3,
                     in_pass2 = names(v) %in% scr$selected_terms)
write.csv(data.frame(index = seq_along(e) - 1, fitted = f, residual = e, standardised = z), file.path(RES, "fig1_insets_residuals.csv"), row.names = FALSE)
write.csv(data.frame(index = seq_along(k) - 1, t = tt, k = k), file.path(RES, "fig1_insets_khat.csv"), row.names = FALSE)
write.csv(vif_df, file.path(RES, "fig1_insets_vif.csv"), row.names = FALSE)
cat(sprintf("chosen pool row %d (seed %d); x3-dot: pass-2 support %d terms, trimmed design %d terms, max k-hat %.3f, %d points with k-hat > 0.7, largest VIF of a true term %.1f (%s), largest VIF overall %.3g\n",
            chosen$row, chosen$seed, length(scr$selected_terms), length(v), max(k), sum(k > 0.7),
            max(v[truth3]), names(which.max(v[truth3])), max(v)))

# ------------------------------------------------------------------ insets
inset_theme <- house_theme() + theme(plot.margin = margin(t = 0.08, r = 0.12, b = 0.05, l = 0.05, unit = "in"))
INSET_W <- 34; INSET_ASPECT <- 0.95

# (i) standardised residuals against the fitted value, binned mean +- 2 sd in 12 equal-count bins (make_figures.R, resid_panel)
ord <- order(f); bins <- split(seq_along(f)[ord], cut(seq_along(f), 12, labels = FALSE))
b <- data.frame(f = sapply(bins, function(i) mean(f[i])), m = sapply(bins, function(i) mean(z[i])), s = sapply(bins, function(i) sd(z[i])))
p_res <- ggplot() +
  geom_point(data = data.frame(f, z), aes(f, z), colour = HOUSE_PRIMARY, alpha = 0.25, size = 0.9) +
  geom_ribbon(data = b, aes(f, ymin = m - 2 * s, ymax = m + 2 * s), fill = HOUSE_HIGHLIGHT, alpha = 0.2) +
  geom_line(data = b, aes(f, m), colour = HOUSE_HIGHLIGHT, linewidth = HOUSE_LINE) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.44) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 3)) +
  labs(x = "fitted value (posterior mean)", y = "standardised residual") + inset_theme +
  theme(axis.title.x = element_text(size = HOUSE_PT$tick))
save_figure(p_res, file.path(FIG, "fig1_inset_residuals.pdf"), width_mm = INSET_W, aspect = INSET_ASPECT, device = cairo_pdf)

# (ii) Pareto k-hat against time, dashed 0.7 line
kd <- data.frame(t = tt, k = k, flag = k > 0.7)
p_k <- ggplot(kd, aes(t, k)) +
  geom_point(data = kd[!kd$flag, ], colour = HOUSE_PRIMARY, size = 0.6, alpha = 0.5) +
  geom_point(data = kd[kd$flag, ], colour = HOUSE_HIGHLIGHT, shape = 17, size = 2.5) +
  geom_hline(yintercept = 0.7, linetype = 2, colour = HOUSE_HIGHLIGHT, linewidth = 0.44) +
  scale_y_continuous(limits = c(min(k) - 0.02, max(0.8, max(k) * 1.05)), breaks = seq(-1, 3, 0.25)) +
  scale_x_continuous(breaks = pretty_breaks(4)) +
  labs(x = "time", y = expression("Pareto " * hat(italic(k))[italic(i)])) + inset_theme
save_figure(p_k, file.path(FIG, "fig1_inset_khat.pdf"), width_mm = INSET_W, aspect = INSET_ASPECT, device = cairo_pdf)

# (iii) VIF of every candidate term of the design entering pass 2; true terms coloured and
# labelled (ggrepel leaders), others grey
suppressPackageStartupMessages(library(ggrepel))
vif_df$colour <- ifelse(vif_df$is_truth, AIZAWA_TRUTH_COL[vif_df$term], HOUSE_GREY)
vt <- vif_df[vif_df$is_truth, ]; vt$label <- term_label(vt$term)
p_v <- ggplot(vif_df, aes(index, vif)) +
  geom_col(data = vif_df[!vif_df$is_truth, ], fill = HOUSE_GREY, width = 0.8) +
  geom_col(data = vt, aes(fill = term), width = 0.8) +
  geom_hline(yintercept = 10, linetype = 2, linewidth = 0.44) +
  geom_text_repel(data = vt, aes(index, vif, label = label), colour = vt$colour, parse = TRUE, size = HOUSE_PT$tick / .pt,
                  family = HOUSE_FONT, nudge_y = 0.8, direction = "both", min.segment.length = 0, segment.size = 0.3,
                  segment.colour = "grey50", box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf, seed = 3, force = 4, force_pull = 0.3) +
  scale_fill_manual(values = AIZAWA_TRUTH_COL, guide = "none") +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.8)) +
  scale_y_log10_powers("VIF", limits = c(1, 10^(ceiling(log10(max(v))) + 1)), step = 1) +
  labs(x = "candidate term") + inset_theme
save_figure(p_v, file.path(FIG, "fig1_inset_vif.pdf"), width_mm = INSET_W, aspect = INSET_ASPECT, device = cairo_pdf)
cat("insets ->", FIG, "; total", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
