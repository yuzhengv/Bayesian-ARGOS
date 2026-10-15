# Supplementary figure S-E (plotting standard): Bayesian-ARGOS with the
# Savitzky-Golay derivative against Bayesian-ARGOS with the weak-form design,
# success rate against n (49 dB, top row) and against SNR (n = 5000, bottom
# row) for Aizawa, Dadras and Rossler, from results/success_2x2_{n,snr}.csv
# (success_2x2.R).  The house-style counterpart of make_figures_2x2.R's
# fig_sg_vs_weak.pdf: theme, palette, n axis with 10^k labels, SNR axis with
# the broken infinity position, the 0.8 line and band (standard 4.4, 4.6):
# Savitzky-Golay arm solid line and filled circle, weak-form arm dashed line
# and open circle, both in the Bayesian-ARGOS blue; one legend row below.  The
# SINDy arms stay internal (make_figures_2x2.R).  Canvas 18 x 11.7 in at
# k = 2.5 (final width 183 mm, aspect 0.65), cairo_pdf.  Writes
# figures/sg_vs_integration.pdf and copies it to <manuscript>/si-plots and
# <manuscript>/reply-to-reviewers/Figures.
# Usage: Rscript make_figures_sg_vs_weak_house.R      (seconds)

suppressPackageStartupMessages({ library(dplyr); library(patchwork) })
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "..", "plot-style", "house_style.R"))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)
OUT <- file.path(FIG, "sg_vs_integration.pdf")

ARMS <- c("Bayesian-ARGOS (SG)", "Bayesian-ARGOS (Integration)")
LABELS <- c("Bayesian-ARGOS, Savitzky-Golay derivative", "Bayesian-ARGOS, weak-form design")
COL <- setNames(rep(HOUSE_COLOURS[["Bayesian-ARGOS"]], 2), ARMS)
LTY <- setNames(c("solid", "dashed"), ARMS); SHP <- setNames(c(16, 1), ARMS)
CASE <- c(aizawa = "Aizawa", dadras = "Dadras", rossler = "Rössler")

dn <- read.csv(file.path(RES, "success_2x2_n.csv"), stringsAsFactors = FALSE) %>% filter(method %in% ARMS) %>%
  mutate(x = as.numeric(x), method = factor(method, ARMS))
ds <- read.csv(file.path(RES, "success_2x2_snr.csv"), stringsAsFactors = FALSE) %>% filter(method %in% ARMS) %>%
  mutate(snr = suppressWarnings(as.numeric(x)), method = factor(method, ARMS)) %>% filter(snr >= 6) %>% mutate(xpos = snr_x(snr))

arm_scales <- list(scale_colour_manual(values = COL, labels = LABELS, name = NULL),
                   scale_linetype_manual(values = LTY, labels = LABELS, name = NULL),
                   scale_shape_manual(values = SHP, labels = LABELS, name = NULL))
panel_n <- function(sys, tag) {
  ggplot(dn[dn$system == sys, ], aes(x, success, colour = method, linetype = method, shape = method)) +
    success_rate_layers() + geom_line(linewidth = HOUSE_CONTEXT_LINE) + geom_point(size = HOUSE_POINT) + arm_scales +
    scale_x_log_n() + scale_y_success() + labs(title = paste0(CASE[sys], ", SNR = 49 dB"), tag = tag) + house_theme()
}
panel_snr <- function(sys, tag) {
  d <- ds[ds$system == sys, ]
  p <- ggplot(d, aes(xpos, success, colour = method, linetype = method, shape = method)) +
    success_rate_layers() + geom_line(data = d[is.finite(d$snr), ], linewidth = HOUSE_CONTEXT_LINE) + geom_point(size = HOUSE_POINT) + arm_scales +
    scale_y_success() + labs(title = as.expression(substitute(S * ", " * italic(n) == 5000, list(S = CASE[[sys]]))), tag = tag) + house_theme()
  snr_axis(p, y_break = 0, extra_y = 0.05)
}
p <- (panel_n("aizawa", "a") + panel_n("dadras", "b") + panel_n("rossler", "c")) /
     (panel_snr("aizawa", "d") + panel_snr("dadras", "e") + panel_snr("rossler", "f")) +
     plot_layout(guides = "collect") &
     theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.width = unit(0.6, "in"),
           legend.text = element_text(size = HOUSE_PT$text), legend.spacing.x = unit(0.3, "in"))
save_figure(p, OUT, width_mm = 183, aspect = 0.65, device = cairo_pdf)
mr <- Sys.getenv("ARGOS_MANUSCRIPT_DIR", normalizePath(file.path(here, "..", "..", "..", "..")))
for (d in c(file.path(mr, "si-plots"), file.path(mr, "reply-to-reviewers", "Figures")))
  if (dir.exists(d)) { file.copy(OUT, file.path(d, basename(OUT)), overwrite = TRUE); cat("copied to", file.path(d, basename(OUT)), "\n") }
cat("wrote", OUT, "\n")
