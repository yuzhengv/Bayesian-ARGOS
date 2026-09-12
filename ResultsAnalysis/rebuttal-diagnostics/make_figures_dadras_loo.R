# Figure for the reply to Referee 1 Point 7 (Dadras PSIS-LOO), house style as
# in make_figures.R.  Reads results/dadras_loo_*.csv (dadras_psis_loo.R) and
# results/dadras_stored_intercept_vs_n.csv (stored 100-trial benchmark results)
# and writes figures/fig_dadras_psis_loo.pdf:
#   a  Pareto k-hat against time for one trial at n = 1e4 and n = 1e5 (same
#      vertical scale, time on the x axis; flagged points k > 0.7 in red)
#   b  the flagged observations of all trials at n = 1e5 on the x1-x3
#      projection of the attractor (clean trajectory of trial 1 in grey)
#   c  the spurious intercept in the stored benchmark: posterior mean of the
#      intercept in the trials that admitted it (median, min-max) against the
#      half-width of the 90 % posterior interval 1.645 sigma / sqrt(n); the
#      number of trials (of 100) keeping the intercept is printed above.
# Usage: Rscript make_figures_dadras_loo.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales); library(latex2exp); library(patchwork)
})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "common.R"))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)

ggplot_theme0 <- theme(
  axis.line = element_line(colour = "black"), axis.ticks.length = unit(.25, "cm"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_blank(), panel.background = element_blank(), legend.key = element_blank(),
  legend.text = element_text(size = 24), legend.title = element_text(face = "bold", size = 26),
  axis.text = element_text(size = 28), axis.title.x = element_text(size = 28),
  axis.title.y = element_text(size = 28, angle = 90, vjust = 0.5),
  plot.title = element_text(size = 26, hjust = 0.5), plot.tag = element_text(size = 28, face = "bold"),
  legend.position = "none", legend.text.align = 0,
  plot.margin = margin(t = 0.75, r = 0.3, b = 0.25, l = 0.1, unit = "in"),
  plot.background = element_rect(fill = 0, colour = 0)
)
RED <- "#d74d3d"; BLUE <- "#505ed5"; GREY <- "grey75"

# ---------------------------------------------------------------- a: k-hat vs time
ser <- read.csv(file.path(RES, "dadras_loo_khat_series.csv")) %>% dplyr::filter(kind == "full", trial == 1, log10n %in% c(4, 5))
ymax <- max(ser$k) * 1.05
khat_panel <- function(l10, tag) {
  s <- ser %>% dplyr::filter(log10n == l10) %>% mutate(flag = k > 0.7)
  ggplot(s, aes(t, k)) +
    geom_point(data = dplyr::filter(s, !flag), colour = BLUE, size = 0.6, alpha = 0.5) +
    geom_point(data = dplyr::filter(s, flag), colour = RED, shape = 17, size = 3.2) +
    geom_hline(yintercept = 0.7, linetype = "dashed", colour = RED, linewidth = 0.8) +
    scale_y_continuous(limits = c(min(ser$k), ymax), breaks = seq(-0.5, 2.5, 0.5)) +
    scale_x_continuous(breaks = pretty_breaks(5)) +
    labs(x = "time", y = TeX("Pareto $\\hat{k}_i$"), tag = tag,
         title = if (l10 == 4) TeX("$n = 10^4$ ($t \\leq 100$)") else TeX("$n = 10^5$ ($t \\leq 1000$)")) + ggplot_theme0
}
pa1 <- khat_panel(4, "a"); pa2 <- khat_panel(5, "")

# ---------------------------------------------------------------- b: flagged points on the attractor
fl <- read.csv(file.path(RES, "dadras_loo_flagged.csv")) %>% dplyr::filter(k > 0.7, log10n == 5)
set.seed(100); ic <- t(sapply(1:100, function(i) runif(3, -4, 4)))[1, ]
traj <- simulate("dadras", 1e5, Inf, ic, 1)$clean
traj <- as.data.frame(traj[seq(1, nrow(traj), by = 2), ]); names(traj) <- c("x1", "x2", "x3")
pb <- ggplot() +
  geom_path(data = traj, aes(x1, x3), colour = GREY, linewidth = 0.25) +
  geom_point(data = fl, aes(x1, x3), colour = RED, shape = 17, size = 3.5) +
  labs(x = TeX("$x_1$"), y = TeX("$x_3$"), tag = "b",
       title = TeX(sprintf("$\\hat{k}_i > 0.7$, %d trials at $n = 10^5$ (%d points)", length(unique(fl$trial)), nrow(fl)))) + ggplot_theme0

# ---------------------------------------------------------------- c: intercept bias vs interval half-width (stored results)
st <- read.csv(file.path(RES, "dadras_stored_intercept_vs_n.csv")) %>% dplyr::filter(log10n >= 4.6)
pc <- ggplot(st, aes(log10n)) +
  geom_ribbon(aes(ymin = int_min, ymax = int_max), fill = alpha(RED, 0.2)) +
  geom_line(aes(y = int_median, colour = "intercept"), linewidth = 1.2) + geom_point(aes(y = int_median, colour = "intercept"), size = 3.5) +
  geom_line(aes(y = halfwidth_90, colour = "halfwidth"), linewidth = 1.2, linetype = "longdash") +
  geom_text(aes(y = int_max * 1.25, label = n_intercept), size = 8, colour = RED) +
  scale_colour_manual(values = c(intercept = RED, halfwidth = BLUE),
                      labels = c(intercept = "posterior mean of the admitted intercept", halfwidth = TeX("$1.645\\,\\sigma / \\sqrt{n}$")), name = NULL) +
  scale_x_continuous(breaks = seq(4.6, 5, 0.1), labels = c(expression(10^4.6), expression(10^4.7), expression(10^4.8), expression(10^4.9), expression(10^5))) +
  scale_y_continuous(limits = c(0, 0.017), breaks = seq(0, 0.016, 0.004)) +
  labs(x = TeX("$n$"), y = "intercept", tag = "c", title = "stored benchmark: admitted intercept (count above)") + ggplot_theme0 +
  theme(legend.position = c(0.55, 0.9), legend.text = element_text(size = 22))

ggsave(file.path(FIG, "fig_dadras_psis_loo.pdf"), (pa1 / pa2) | pb | pc, width = 30, height = 11, dpi = 300, units = "in")
cat("wrote", file.path(FIG, "fig_dadras_psis_loo.pdf"), "\n")
