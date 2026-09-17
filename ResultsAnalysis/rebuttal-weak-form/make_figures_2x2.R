# Figures for the derivative-vs-integral comparison: success rate against n
# (49 dB) and against SNR (n = 5000) for Aizawa, Dadras and Rossler.
#   fig_sg_vs_weak.pdf                 Bayesian-ARGOS with the Savitzky-Golay step
#                                      vs with the weak-form design (the SI figure)
#   fig_success_2x2_internal.pdf       the same plus the two SINDy arms, internal
#                                      reference only (decision 2026-09-16)  House style of ../rebuttal-diagnostics/make_figures.R (Fig. 2-3 of
# the manuscript: sans-serif, two axis borders, no grid, ARGOS palette, the
# noiseless level on a broken SNR axis).  Reads results/success_2x2_{n,snr}.csv
# written by success_2x2.R.  Usage: Rscript make_figures_2x2.R
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork); library(latex2exp); library(ggh4x) })
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)

th <- theme(axis.line = element_line(colour = "black"), axis.ticks.length = unit(.2, "cm"),
            panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
            legend.key = element_blank(), legend.title = element_blank(), legend.position = "top",
            text = element_text(size = 15), plot.tag = element_text(face = "bold", size = 18),
            plot.title = element_text(size = 15, hjust = 0.5))
METHODS <- c("Bayesian-ARGOS (SG)", "Bayesian-ARGOS (Integration)", "SINDy (SG)", "SINDy (Integration)")
LABELS <- c("Bayesian-ARGOS (SG)", "Bayesian-ARGOS (Integration)", "SINDy (SG)", "SINDy (Integration; PySINDy weak library, baseline STLSQ settings)")
COL <- c("Bayesian-ARGOS (SG)" = "#505ed5", "Bayesian-ARGOS (Integration)" = "#505ed5", "SINDy (SG)" = "#8eb63d", "SINDy (Integration)" = "#8eb63d")
LTY <- c("Bayesian-ARGOS (SG)" = "solid", "Bayesian-ARGOS (Integration)" = "dashed", "SINDy (SG)" = "solid", "SINDy (Integration)" = "dashed")
SHP <- c("Bayesian-ARGOS (SG)" = 16, "Bayesian-ARGOS (Integration)" = 1, "SINDy (SG)" = 17, "SINDy (Integration)" = 2)
CASE <- c(aizawa = "Aizawa", dadras = "Dadras", rossler = "Rössler")

dn_all <- read.csv(file.path(RES, "success_2x2_n.csv")) %>% mutate(x = as.numeric(x), method = factor(method, METHODS))
ds_all <- read.csv(file.path(RES, "success_2x2_snr.csv")) %>% mutate(snr = suppressWarnings(as.numeric(x)), method = factor(method, METHODS))
INF_X <- 73
ds_all$xpos <- ifelse(is.finite(ds_all$snr), ds_all$snr, INF_X)

# SNR axis as in Fig. 3a: finite values in dB, infinity behind an axis break
create_separators <- function(x, extra_x, y, extra_y, angle = 45, scale = 1, length = .1) {
  add_y <- length * sin(angle * pi / 180) / 2; add_x <- length * cos(angle * pi / 180)
  list(x = x - add_x * scale, xend = x + add_x * scale + extra_x,
       y = rep(y - add_y * scale - extra_y, length(x)), yend = rep(y + add_y * scale - extra_y / 2, length(x)))
}
snr_axis <- function(p, ticks = seq(10, 60, 10)) {
  seg <- create_separators(c(65.5, 69.5), extra_x = 1, y = 0, extra_y = 0.05, angle = 75)
  p + scale_x_continuous(limits = c(0, INF_X + 1), breaks = c(ticks, INF_X), labels = c(ticks, TeX("$\\infty$"))) +
    guides(x = guide_axis_truncated(trunc_lower = c(-Inf, 70), trunc_upper = c(66, Inf))) +
    annotate("segment", x = seg$x, xend = seg$xend, y = seg$y, yend = seg$yend) + coord_cartesian(clip = "off")
}

panel_n <- function(sys) {
  ggplot(dn[dn$system == sys, ], aes(x, success, colour = method, linetype = method, shape = method)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
    scale_colour_manual(values = COL, labels = LABELS) + scale_linetype_manual(values = LTY, labels = LABELS) + scale_shape_manual(values = SHP, labels = LABELS) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
    labs(x = TeX("$\\log_{10} n$"), y = "success rate", title = paste0(CASE[sys], ", 49 dB")) + th
}
panel_snr <- function(sys) {
  p <- ggplot(ds[ds$system == sys, ], aes(xpos, success, colour = method, linetype = method, shape = method)) +
    geom_line(data = ds[ds$system == sys & is.finite(ds$snr), ], linewidth = 0.9) + geom_point(size = 2.2) +
    scale_colour_manual(values = COL, labels = LABELS) + scale_linetype_manual(values = LTY, labels = LABELS) + scale_shape_manual(values = SHP, labels = LABELS) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
    labs(x = "SNR (dB)", y = "success rate", title = paste0(CASE[sys], ", n = 5000")) + th
  snr_axis(p)
}
draw <- function(keep, file) {
  dn <<- dn_all[dn_all$method %in% keep, ]; ds <<- ds_all[ds_all$method %in% keep, ]
  p <- (panel_n("aizawa") + panel_n("dadras") + panel_n("rossler")) /
       (panel_snr("aizawa") + panel_snr("dadras") + panel_snr("rossler")) +
       plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
  ggsave(file.path(FIG, file), p, width = 16, height = 9.5, device = cairo_pdf)
  cat("wrote", file.path(FIG, file), "\n")
}
draw(METHODS[1:2], "fig_sg_vs_weak.pdf")
draw(METHODS, "fig_success_2x2_internal.pdf")
