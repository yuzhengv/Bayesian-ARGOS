# Figures for the reply to Referee 1 Point 6, drawn in the manuscript's house
# style (the ggplot theme, palette, fonts and broken-SNR-axis construction of
# ResultsAnalysis/success-rate-plots/*_results_compare.ipynb) and following
# Referee 1's figure requirements:
#   * Point 6  -- highlight only the curves that matter (true terms in colour,
#                 every other library term in grey, no markers);
#   * minor 4  -- one notation for the noiseless case (the infinity symbol on a
#                 broken axis) and SNR ticks in steps of 10 dB;
#   * minor 5  -- one visual style: sans-serif, no grid, two axis borders only,
#                 same tick placement, same colour scheme as Fig. 3a;
#   * minor 6  -- subplots that are compared share the same vertical scale and
#                 are labelled with system and condition.
#
# Reads results/*.csv produced by exp1_collinearity.R and
# exp2_heteroscedasticity.R and writes figures/*.pdf.
# Usage: Rscript make_figures.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(latex2exp); library(ggh4x); library(patchwork)
})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)

# ---------------------------------------------------------------- house style
ggplot_theme0 <- theme(
  axis.line = element_line(colour = "black"),
  axis.ticks.length = unit(.25, "cm"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_blank(), panel.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 24), legend.title = element_text(face = "bold", size = 26),
  axis.text = element_text(size = 28), axis.title.x = element_text(size = 28),
  axis.title.y = element_text(size = 28, angle = 90, vjust = 0.5),
  plot.title = element_text(size = 26, hjust = 0.5), plot.tag = element_text(size = 28, face = "bold"),
  legend.position = "none", legend.text.align = 0,
  plot.margin = margin(t = 0.75, r = 0.3, b = 0.25, l = 0.1, unit = "in")
)
ggplot_theme1 <- ggplot_theme0 + theme(plot.background = element_rect(fill = 0, colour = 0))
HOUSE <- c("#d74d3d", "#505ed5", "#8eb63d", "#FA5524", "#EF60AD")   # ARGOS red, Bayesian-ARGOS blue, SINDy green, ...
GREY <- "grey75"

create_separators <- function(x, extra_x, y, extra_y, angle = 45, scale = 1, length = .1) {
  add_y <- length * sin(angle * pi / 180) / 2; add_x <- length * cos(angle * pi / 180)
  list(x = x - add_x * scale, xend = x + add_x * scale + extra_x,
       y = rep(y - add_y * scale - extra_y, length(x)), yend = rep(y + add_y * scale - extra_y / 2, length(x)))
}

# SNR axis as in Fig. 3a: finite values in dB, infinity at x = 73 behind an axis break
INF_X <- 73
snr_x <- function(s) ifelse(is.infinite(s), INF_X, s)
snr_axis <- function(p, y_break, ticks = seq(10, 60, 10), extra_y = 0.1) {
  xstart <- 65.5; xend <- 69.5; extra_x <- 1
  seg <- create_separators(c(xstart, xend), extra_x = extra_x, y = y_break, extra_y = extra_y, angle = 75)
  p + scale_x_continuous(limits = c(min(ticks) - 4, INF_X + 1), breaks = c(ticks, INF_X),
                         labels = c(ticks, TeX("$\\infty$"))) +
    guides(x = guide_axis_truncated(trunc_lower = c(-Inf, xend + extra_x / 2), trunc_upper = c(xstart + extra_x / 2, Inf))) +
    annotate("segment", x = seg$x, xend = seg$xend, y = seg$y, yend = seg$yend) +
    coord_cartesian(clip = "off")
}
log_n_scale <- scale_x_continuous(breaks = seq(2, 5, 0.5),
  labels = c(expression(10^2), expression(10^2.5), expression(10^3), expression(10^3.5), expression(10^4), expression(10^4.5), expression(10^5)))
log10_y <- function(name, limits = NULL, step = 1) scale_y_continuous(name = name, trans = "log10", breaks = 10^seq(-6, 16, step),
  labels = trans_format("log10", math_format(10^.x)), limits = limits)

# term label "x1^2 x3" -> plotmath x[1]^2 * x[3]
term_pm <- function(t) {
  toks <- strsplit(t, " ")[[1]]
  paste(sapply(toks, function(k) { v <- sub("\\^.*", "", k); p <- if (grepl("\\^", k)) sub(".*\\^", "", k) else NULL
    paste0("x[", substring(v, 2), "]", if (!is.null(p)) paste0("^", p) else "") }), collapse = "*")
}
TRUTH <- c("x3", "x1^2", "x2^2", "x1^2 x3", "x2^2 x3", "x1^3 x3", "x3^3")
TRUTH_COL <- setNames(c(HOUSE, "#009E9E", "#8B5A2B"), TRUTH)

# ======================================================= Fig. 3b replacement
vif <- read.csv(file.path(RES, "exp1_vif_long.csv"))
med <- vif %>% filter(design == "trimmed") %>% group_by(sweep, term, log10n, snr) %>%
  summarise(vif = median(vif, na.rm = TRUE), .groups = "drop") %>% mutate(is_truth = term %in% TRUTH)
vif_panel <- function(sw, xvar, xlab) {
  d <- med %>% filter(sweep == sw) %>% mutate(x = if (sw == "snr") snr_x(snr) else log10n)
  others <- d %>% filter(!is_truth); truth <- d %>% filter(is_truth) %>% mutate(term = factor(term, levels = TRUTH))
  p <- ggplot() +
    geom_line(data = others, aes(x, vif, group = term), colour = GREY, linewidth = 0.7) +
    geom_line(data = truth, aes(x, vif, colour = term), linewidth = 1.6) +
    geom_hline(yintercept = 10, linetype = 2) +
    scale_colour_manual(values = TRUTH_COL, labels = sapply(TRUTH, function(t) parse(text = term_pm(t))[[1]]), name = "true terms") +
    log10_y("VIF", step = if (sw == "n") 4 else 1) + labs(x = xlab) + ggplot_theme1 +
    theme(legend.position = if (sw == "n") c(0.8, 0.72) else "none", legend.text = element_text(size = 22), legend.title = element_text(size = 22)) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(linewidth = 3)))
  if (sw == "snr") snr_axis(p, y_break = 1, extra_y = 0) else p + log_n_scale
}
pa <- vif_panel("n", "log10n", expression(italic("n"))) + labs(tag = "a", title = "SNR = 49 dB")
pb <- vif_panel("snr", "snr", TeX("SNR (dB)")) + labs(tag = "b", title = expression(italic("n") == 5000))
ggsave(file.path(FIG, "fig_vif_redesign.pdf"), pa + pb, width = 18, height = 6.5, dpi = 300, units = "in")

# ================================================= collinearity: source and n
kappa <- read.csv(file.path(RES, "exp1_kappa.csv")); geo <- read.csv(file.path(RES, "exp1_geometry.csv"))
DES <- c(full = "full library", trimmed = "degree-trimmed library (Fig. 3b)", selected = "screened support", truth = "true terms only")
DES_COL <- c(full = HOUSE[1], trimmed = HOUSE[4], selected = HOUSE[3], truth = HOUSE[2]); DES_SHAPE <- c(full = 15, trimmed = 17, selected = 16, truth = 18)
ks <- kappa %>% group_by(sweep, design, log10n, snr) %>%
  summarise(m = median(kappa), lo = quantile(kappa, .25), hi = quantile(kappa, .75), .groups = "drop") %>%
  mutate(design = factor(design, levels = names(DES)))
kappa_panel <- function(sw, xlab) {
  d <- ks %>% filter(sweep == sw) %>% mutate(x = if (sw == "snr") snr_x(snr) else log10n)
  p <- ggplot(d, aes(x, m, colour = design, fill = design, shape = design)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 1.4) + geom_point(size = 4) +
    scale_colour_manual(values = DES_COL, labels = DES, name = NULL) + scale_fill_manual(values = DES_COL, labels = DES, name = NULL) +
    scale_shape_manual(values = DES_SHAPE, labels = DES, name = NULL) +
    log10_y(expression("condition number " * kappa), step = if (sw == "n") 2 else 1) + labs(x = xlab) + ggplot_theme1
  if (sw == "snr") snr_axis(p, y_break = 1, extra_y = 0) else p + log_n_scale
}
pk1 <- kappa_panel("n", expression(italic("n"))) + labs(tag = "a", title = "SNR = 49 dB") +
  theme(legend.position = c(0.62, 0.85), legend.text = element_text(size = 20))
pk2 <- kappa_panel("snr", TeX("SNR (dB)")) + labs(tag = "b", title = expression(italic("n") == 5000))
SAMP <- c(trajectory_sg = "one\ntrajectory\n(noisy, SG)", trajectory_clean = "one\ntrajectory\n(clean)",
          multi_transient_sg = "50\ntransients\n(noisy, SG)", multi_transient = "50\ntransients\n(clean)", box_uniform = "uniform\nbox\nsamples")
gs <- geo %>% filter(design %in% c("full", "truth")) %>% group_by(sampling, design) %>%
  summarise(m = median(kappa), lo = quantile(kappa, .25), hi = quantile(kappa, .75), .groups = "drop") %>%
  mutate(sampling = factor(sampling, levels = names(SAMP)), design = factor(design, levels = c("full", "truth")))
pk3 <- ggplot(gs, aes(sampling, m, fill = design)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8), width = 0.25, linewidth = 0.8) +
  scale_fill_manual(values = DES_COL[c("full", "truth")], labels = DES[c("full", "truth")], name = NULL) +
  scale_x_discrete(labels = SAMP) + log10_y(expression("condition number " * kappa)) +
  labs(x = NULL, tag = "c", title = expression(italic("n") == 5000)) + ggplot_theme1 +
  theme(axis.text.x = element_text(size = 19, lineheight = 0.9), legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 22))
# Reply / SI version: panels a and b only (library designs on the single trajectory).
# Panel c (one trajectory vs 50 transients vs uniform box samples) is WITHHELD from the
# rebuttal by decision of 2026-09-11 (see notes/rebuttal-memos/reviewer1/
# transient-experiment-withheld.md in the manuscript repo); it is kept for internal use
# in the *_internal.pdf file only.
ggsave(file.path(FIG, "fig_collinearity_source.pdf"), pk1 + pk2,
       width = 19, height = 7.5, dpi = 300, units = "in")
ggsave(file.path(FIG, "fig_collinearity_source_internal.pdf"), pk1 + pk2 + pk3 + plot_layout(widths = c(1, 1, 1.35)),
       width = 28, height = 7.5, dpi = 300, units = "in")

# ====================================================== residual structure
stats <- read.csv(file.path(RES, "exp2_stats.csv"))
ex <- read.csv(file.path(RES, "exp2_examples.csv"), stringsAsFactors = FALSE)
SER <- c(aizawa_true = "Aizawa, true support", aizawa_sel = "Aizawa, screened support", rossler_true = "Rössler, true support")
SER_COL <- c(aizawa_true = HOUSE[2], aizawa_sel = HOUSE[4], rossler_true = HOUSE[1]); SER_SHAPE <- c(aizawa_true = 16, aizawa_sel = 17, rossler_true = 15)
long <- bind_rows(
  stats %>% transmute(series = paste0(system, "_true"), snr, rho1, var_ratio),
  stats %>% filter(system == "aizawa") %>% transmute(series = "aizawa_sel", snr, rho1 = sel_rho1, var_ratio = sel_var_ratio)) %>%
  group_by(series, snr) %>% summarise(across(c(rho1, var_ratio), list(m = median, lo = ~quantile(.x, .25), hi = ~quantile(.x, .75))), .groups = "drop") %>%
  mutate(x = snr_x(snr), series = factor(series, levels = names(SER)))
stat_panel <- function(ycol, ylab, ylog = FALSE, y_break = 0) {
  p <- ggplot(long, aes(x, .data[[paste0(ycol, "_m")]], colour = series, fill = series, shape = series)) +
    geom_ribbon(aes(ymin = .data[[paste0(ycol, "_lo")]], ymax = .data[[paste0(ycol, "_hi")]]), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.4) + geom_point(size = 4) +
    scale_colour_manual(values = SER_COL, labels = SER, name = NULL) + scale_fill_manual(values = SER_COL, labels = SER, name = NULL) +
    scale_shape_manual(values = SER_SHAPE, labels = SER, name = NULL) + labs(x = TeX("SNR (dB)")) + ggplot_theme1
  p <- if (ylog) p + log10_y(ylab) else p + scale_y_continuous(name = ylab, limits = c(0, 1), breaks = seq(0, 1, 0.25))
  snr_axis(p, y_break = y_break, extra_y = if (ylog) 0 else 0.05)
}
ph1 <- stat_panel("rho1", expression("lag-1 residual autocorrelation " * rho[1]), FALSE, y_break = 0) + labs(tag = "a") +
  theme(legend.position = c(0.42, 0.25), legend.text = element_text(size = 21))
ph2 <- stat_panel("var_ratio", "residual variance ratio (deciles)", TRUE, y_break = 1) + labs(tag = "b")
comp <- stats %>% filter(system == "aizawa") %>% transmute(snr, d = var_d / var_e, g = var_g / var_e, h = var_h / var_e) %>%
  pivot_longer(-snr, names_to = "part") %>% group_by(snr, part) %>% summarise(m = median(value), .groups = "drop") %>%
  mutate(x = snr_x(snr), part = factor(part, levels = c("d", "g", "h")))
PART <- c(d = "derivative estimation", g = "state smoothing", h = "coefficient error"); PART_COL <- c(d = HOUSE[3], g = HOUSE[5], h = "#8B5A2B")
ph3 <- snr_axis(ggplot(comp, aes(x, m, colour = part, shape = part)) + geom_line(linewidth = 1.4) + geom_point(size = 4) +
  scale_colour_manual(values = PART_COL, labels = PART, name = NULL) + scale_shape_manual(values = c(16, 17, 15), labels = PART, name = NULL) +
  log10_y("share of residual variance", limits = c(1e-5, 3)) + labs(x = TeX("SNR (dB)"), tag = "c", title = "Aizawa, true support") + ggplot_theme1 +
  theme(legend.position = c(0.45, 0.6), legend.direction = "horizontal", legend.text = element_text(size = 19)), y_break = 1e-5, extra_y = 0)

# residual-vs-fitted panels, standardised residuals so that the three panels share one scale
arr <- function(system, snr, key) {
  s_num <- if (identical(snr, "inf")) Inf else as.numeric(snr)
  v <- ex %>% filter(.data$system == !!system, .data$snr == s_num, .data$key == !!key) %>% arrange(index)
  as.numeric(v$value)
}
resid_panel <- function(system, snr, title, tag) {
  e <- arr(system, snr, "residuals"); f <- arr(system, snr, "fitted"); z <- e / sd(e)
  ord <- order(f); bins <- split(seq_along(f)[ord], cut(seq_along(f), 12, labels = FALSE))
  b <- data.frame(f = sapply(bins, function(i) mean(f[i])), m = sapply(bins, function(i) mean(z[i])), s = sapply(bins, function(i) sd(z[i])))
  ggplot() + geom_point(data = data.frame(f, z), aes(f, z), colour = HOUSE[2], alpha = 0.25, size = 0.9) +
    geom_ribbon(data = b, aes(f, ymin = m - 2 * s, ymax = m + 2 * s), fill = HOUSE[1], alpha = 0.2) +
    geom_line(data = b, aes(f, m), colour = HOUSE[1], linewidth = 1.6) + geom_hline(yintercept = 0, linetype = 3) +
    scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 3)) +
    labs(x = "fitted value (posterior mean)", y = "standardised residual", title = title, tag = tag) + ggplot_theme1
}
ph4 <- resid_panel("aizawa", "60", "Aizawa, SNR = 60 dB", "d")
ph5 <- resid_panel("aizawa", "inf", expression("Aizawa, SNR = " * infinity), "e")
ph6 <- resid_panel("rossler", "inf", expression("Rössler, SNR = " * infinity), "f")
ggsave(file.path(FIG, "fig_heteroscedasticity.pdf"), (ph1 + ph2 + ph3) / (ph4 + ph5 + ph6), width = 27, height = 14, dpi = 300, units = "in")

# ================================= residual vs derivative-estimation error (time window)
win <- 1001:1600
overlay <- function(system, snr, title, tag) {
  e <- arr(system, snr, "residuals"); d <- arr(system, snr, if (snr == "inf") "d_trunc" else "d"); s <- sd(e)
  dd <- data.frame(t = rep((win - 1) * 0.01, 2), v = c(e[win], d[win]) / s,
                   what = rep(c("residual", if (snr == "inf") "SG truncation error" else "derivative-estimation error"), each = length(win)))
  ggplot(dd, aes(t, v, colour = what, linetype = what)) + geom_line(linewidth = 1.2) +
    scale_colour_manual(values = c(HOUSE[2], HOUSE[1]), name = NULL) + scale_linetype_manual(values = c(1, 2), name = NULL) +
    scale_y_continuous(limits = c(-4, 4)) + labs(x = "time", y = "value / sd(residual)", title = title, tag = tag) + ggplot_theme1 +
    theme(legend.position = c(0.75, 0.9), legend.text = element_text(size = 20))
}
po1 <- overlay("aizawa", "50", "Aizawa, SNR = 50 dB", "a"); po2 <- overlay("aizawa", "inf", expression("Aizawa, SNR = " * infinity), "b")
ggsave(file.path(FIG, "fig_residual_vs_derivative_error.pdf"), po1 + po2, width = 18, height = 6.5, dpi = 300, units = "in")
cat("figures ->", FIG, "\n")
