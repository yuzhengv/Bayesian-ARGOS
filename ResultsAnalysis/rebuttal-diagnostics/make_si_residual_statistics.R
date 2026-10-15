# Supplementary figure S-D (plotting standard): residual structure against SNR
# and the Aizawa residual panels, the panels a to e of make_figures.R's
# fig_heteroscedasticity restyled with ../plot-style/house_style.R (the old
# panel f, Rossler at infinity, is in Fig. 5 of the manuscript).
#   a  lag-1 residual autocorrelation against SNR (medians, interquartile
#      ribbons; Aizawa x3-dot true support, Aizawa x3-dot screened support,
#      Rossler x1-dot true support; results/exp2_stats.csv)
#   b  max / min residual variance over deciles of the fitted value against SNR
#   c  share of the residual variance explained by derivative estimation, state
#      smoothing and coefficient error (Aizawa x3-dot, true support)
#   d  Aizawa x3-dot standardised residuals against the fitted value at 60 dB
#   e  the same at infinity (results/exp2_examples.csv, example trial)
# In-panel legends in the same corner (top left) of a, b and c.  Canvas
# 18 x 12.6 in at k = 2.5 (final width 183 mm, aspect 0.7), cairo_pdf.  Writes
# figures/residual_statistics.pdf and copies it to <manuscript>/si-plots and
# <manuscript>/reply-to-reviewers/Figures.
# Usage: Rscript make_si_residual_statistics.R      (seconds)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(patchwork) })
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "..", "plot-style", "house_style.R"))
source(file.path(here, "fig_helpers.R"))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)
OUT <- file.path(FIG, "residual_statistics.pdf")

stats <- read.csv(file.path(RES, "exp2_stats.csv"), stringsAsFactors = FALSE)
ex <- read.csv(file.path(RES, "exp2_examples.csv"), stringsAsFactors = FALSE)
SER <- c(aizawa_true = 'paste("Aizawa ", dot(italic(x))[3], ", true support")',
         aizawa_sel = 'paste("Aizawa ", dot(italic(x))[3], ", screened support")',
         rossler_true = 'paste("Rössler ", dot(italic(x))[1], ", true support")')
SER_COL <- c(aizawa_true = HOUSE_PRIMARY, aizawa_sel = HOUSE_EXTRA[1], rossler_true = HOUSE_HIGHLIGHT)
SER_SHAPE <- c(aizawa_true = 16, aizawa_sel = 17, rossler_true = 15)
long <- bind_rows(
  stats %>% transmute(series = paste0(system, "_true"), snr, rho1, var_ratio),
  stats %>% filter(system == "aizawa") %>% transmute(series = "aizawa_sel", snr, rho1 = sel_rho1, var_ratio = sel_var_ratio)) %>%
  group_by(series, snr) %>%
  summarise(across(c(rho1, var_ratio), list(m = median, lo = ~quantile(.x, .25), hi = ~quantile(.x, .75))), .groups = "drop") %>%
  filter(snr >= 6) %>% mutate(x = snr_x(snr), series = factor(series, levels = names(SER)))
inside_legend <- theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
                       legend.key.width = unit(0.5, "in"), legend.key.height = unit(0.28, "in"), legend.text = element_text(size = HOUSE_PT$text))
stat_panel <- function(ycol, ylab, ylog = FALSE, y_break = 0, ylim = NULL) {
  fin <- long %>% filter(is.finite(snr))
  p <- ggplot(long, aes(x, .data[[paste0(ycol, "_m")]], colour = series, fill = series, shape = series)) +
    geom_ribbon(data = fin, aes(ymin = .data[[paste0(ycol, "_lo")]], ymax = .data[[paste0(ycol, "_hi")]]), alpha = 0.15, colour = NA) +
    geom_line(data = fin, linewidth = HOUSE_LINE) + geom_point(size = HOUSE_MARKER) +
    geom_linerange(data = long %>% filter(is.infinite(snr)), aes(ymin = .data[[paste0(ycol, "_lo")]], ymax = .data[[paste0(ycol, "_hi")]]), linewidth = 0.6, alpha = 0.6) +
    scale_colour_manual(values = SER_COL, labels = lapply(SER, function(s) parse(text = s)[[1]]), name = NULL) +
    scale_fill_manual(values = SER_COL, labels = lapply(SER, function(s) parse(text = s)[[1]]), name = NULL) +
    scale_shape_manual(values = SER_SHAPE, labels = lapply(SER, function(s) parse(text = s)[[1]]), name = NULL) +
    labs(title = expression(italic(n) == 5000)) + house_theme("inside") + inside_legend
  p <- if (ylog) p + scale_y_log10_powers(ylab, limits = ylim) else
    p + scale_y_continuous(name = ylab, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02)), oob = oob_keep)
  if (ylog) snr_axis_log(p, ylim) else snr_axis(p, y_break = y_break, extra_y = 0.05)
}
pa <- stat_panel("rho1", expression("lag-1 residual autocorrelation " * italic(rho)[1]), FALSE, y_break = 0) + labs(tag = "a")
pb <- stat_panel("var_ratio", "residual variance ratio (deciles)", TRUE, y_break = 1, ylim = c(1, 1e6)) + labs(tag = "b")

comp <- stats %>% filter(system == "aizawa", snr >= 6) %>% transmute(snr, d = var_d / var_e, g = var_g / var_e, h = var_h / var_e) %>%
  pivot_longer(-snr, names_to = "part") %>% group_by(snr, part) %>% summarise(m = median(value), .groups = "drop") %>%
  mutate(x = snr_x(snr), part = factor(part, levels = c("d", "g", "h")))
PART <- c(d = "derivative estimation", g = "state smoothing", h = "coefficient error")
PART_COL <- c(d = HOUSE_COLOURS[["SINDy"]], g = HOUSE_EXTRA[2], h = HOUSE_EXTRA[4]); PART_SHAPE <- c(d = 16, g = 17, h = 15)
C_LIM <- c(1e-5, 1e2)
pc <- snr_axis_log(ggplot(comp, aes(x, m, colour = part, shape = part)) +
  geom_line(data = comp %>% filter(is.finite(snr)), linewidth = HOUSE_LINE) + geom_point(size = HOUSE_MARKER) +
  scale_colour_manual(values = PART_COL, labels = PART, name = NULL) + scale_shape_manual(values = PART_SHAPE, labels = PART, name = NULL) +
  scale_y_log10_powers("share of residual variance", limits = C_LIM) +
  labs(tag = "c", title = eq_title("Aizawa", 3, quote("true support, " * italic(n) == 5000))) + house_theme("inside") + inside_legend,
  C_LIM)

arr <- function(system, snr, key) {                       # snr: a number or "inf" (read.csv parses "inf" as Inf)
  s_num <- if (identical(snr, "inf")) Inf else as.numeric(snr)
  v <- ex[ex$system == system & ex$snr == s_num & ex$key == key, ]
  as.numeric(v$value[order(v$index)])
}
pd <- resid_panel(arr("aizawa", "60", "residuals"), arr("aizawa", "60", "fitted"), eq_title("Aizawa", 3, snr_cond(60)), "d")
pe <- resid_panel(arr("aizawa", "inf", "residuals"), arr("aizawa", "inf", "fitted"), eq_title("Aizawa", 3, snr_cond(Inf)), "e")

p <- pa + pb + pc + pd + pe + plot_layout(design = "AABBCC\nDDDEEE", heights = c(1, 1))
save_figure(p, OUT, width_mm = 183, aspect = 0.7, device = cairo_pdf)
copy_to_manuscript(OUT, "si-plots")
cat("wrote", OUT, "\n")
