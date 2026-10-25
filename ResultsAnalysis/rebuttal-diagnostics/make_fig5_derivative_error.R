# Revised Fig. 5 of the manuscript (plotting standard block M5): mechanism II,
# the deterministic error of the Savitzky-Golay derivative made significant by a
# narrowing credible interval (Rossler x1-dot as the noise vanishes, Dadras
# x1-dot as n grows). The panels follow the order of the main text, the
# Rossler (noise) route first and the Dadras (n) route second.
#   row 1 (symptom)    a  term-selection frequency of the final models for Rossler
#                         x1-dot at 61 dB and infinity (n = 5000; results/rossler_
#                         eq1_term_frequencies.csv), b  the same for Dadras x1-dot
#                         at n = 10^4 and 10^5 (49 dB; results/dadras_eq1_term_
#                         frequencies.csv); terms with a frequency of at least
#                         0.05 at any shown condition; true terms blue, the
#                         spurious intercept red (named in the key and on the axis)
#   row 2 (diagnostic) c, d  Rossler x1-dot standardised residuals against the
#                         fitted value at 61 dB and at infinity (results/exp2_
#                         examples.csv, example trial), one scale, binned mean
#                         +- 2 sd in 12 equal-count bins, dotted zero line
#   row 3 (diagnostic) e, f  Dadras x1-dot Pareto k-hat against time as the block
#                         maximum over 100 observations (trial 6, n = 10^4 and
#                         10^5, one vertical scale, results/dadras_loo_khat_
#                         series.csv), blocks with a flagged point k-hat > 0.7 as
#                         red triangles, dashed 0.7 line
#   row 4              g  the flagged observations of all 20 trials at n = 10^5
#                         on the x1-x3 projection of the attractor (clean
#                         trajectory of trial 1 in grey, simulated through
#                         common.R as in make_figures_dadras_loo.R; results/
#                         dadras_loo_flagged.csv)
#       (mechanism)    h  the admitted intercept in the stored benchmark against n
#                         (posterior mean in the trials that admit it: median and
#                         range, the count of such trials printed above) against
#                         the 90 % half-width 1.645 sigma / sqrt(n) (results/
#                         dadras_stored_intercept_vs_n.csv)
# One vertical key per row at the left of the two panels (standard 4.7.2).
# Canvas 18 x 20.7 in at k = 2.5 (final width 183 mm, aspect 1.15), cairo_pdf.
# Writes figures/diagnostics_derivative_error.pdf and copies it to
# <manuscript>/imgs and <manuscript>/reply-to-reviewers/Figures.
# Usage: Rscript make_fig5_derivative_error.R      (about 20 s; simulates one
#   clean Dadras trajectory through ode_auto.py)

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
suppressPackageStartupMessages({ library(patchwork) })
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "..", "plot-style", "house_style.R"))
source(file.path(here, "fig_helpers.R"))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)
OUT <- file.path(FIG, "diagnostics_derivative_error.pdf")
K_FLAG <- 0.7

# ------------------------------------------------------------- row 1: symptom
fd <- read.csv(file.path(RES, "dadras_eq1_term_frequencies.csv"), stringsAsFactors = FALSE) %>%
  dplyr::filter(sweep == "n", value %in% c(4, 5)) %>% mutate(cond = factor(value, levels = c(4, 5)))
fr <- read.csv(file.path(RES, "rossler_eq1_term_frequencies.csv"), stringsAsFactors = FALSE) %>%
  dplyr::filter(sweep == "snr", value %in% c(61, Inf)) %>% mutate(cond = factor(value, levels = c(61, Inf)))
DAD_TRUTH <- c("x1", "x2", "x2x3"); ROS_TRUTH <- c("x2", "x3")
blue <- function(terms) setNames(rep(HOUSE_PRIMARY, length(terms)), terms)
# the spurious intercept is identified by its colour (key at the left) and its tick label
pa1 <- skyline_panel(fd, DAD_TRUTH, blue(DAD_TRUTH), label_terms = unique(fd$term), min_freq = 0.05,
                     cond_labels = c("italic(n) == 10^4", "italic(n) == 10^5"), order_by = 5,
                     intercept_col = HOUSE_HIGHLIGHT, intercept_first = FALSE, bar_width = 0.6, label_angle = 0) +
  labs(title = eq_title("Dadras", 1, quote("SNR = 49 dB")), tag = "b")
pa2 <- skyline_panel(fr, ROS_TRUTH, blue(ROS_TRUTH), label_terms = unique(fr$term), min_freq = 0.05,
                     cond_labels = c('paste("SNR = 61 dB")', 'paste("SNR = ", infinity)'), order_by = Inf,
                     intercept_col = HOUSE_HIGHLIGHT, intercept_first = FALSE, bar_width = 0.6, label_angle = 0) +
  labs(title = eq_title("Rössler", 1, quote(italic(n) == 5000)), tag = "a")
key_1 <- key_panel(data.frame(label = c("true term", "spurious intercept"), colour = c(HOUSE_PRIMARY, HOUSE_HIGHLIGHT), type = "fill"))

# ------------------------------------------- row 2: Rossler residuals (exp2)
ex <- read.csv(file.path(RES, "exp2_examples.csv"), stringsAsFactors = FALSE)
arr <- function(system, snr, key) {                       # snr: a number or "inf" (read.csv parses "inf" as Inf)
  s_num <- if (identical(snr, "inf")) Inf else as.numeric(snr)
  v <- ex[ex$system == system & ex$snr == s_num & ex$key == key, ]
  as.numeric(v$value[order(v$index)])
}
pb1 <- resid_panel(arr("rossler", "61", "residuals"), arr("rossler", "61", "fitted"), eq_title("Rössler", 1, snr_cond(61)), "c")
pb2 <- resid_panel(arr("rossler", "inf", "residuals"), arr("rossler", "inf", "fitted"), eq_title("Rössler", 1, snr_cond(Inf)), "d")
key_2 <- key_panel(data.frame(label = c("residual", "binned mean \u00b1 2 sd\n(12 bins)", "zero line"),
                              colour = c(HOUSE_PRIMARY, HOUSE_HIGHLIGHT, "black"), type = c("point", "line", "dotted"),
                              shape = c(16, NA, NA), band = c(FALSE, TRUE, FALSE)), point_size = 2, slots = 8)

# ---------------------------------------------- row 3: Dadras k-hat vs time
# Block maximum of k-hat over 100 consecutive observations, stored for every trial
# of the PSIS-LOO run, for a trial with flagged observations at both n. The full
# per-observation series is stored for trial 1 only, which has none at n = 10^4.
SERIES_TRIAL <- 6L
ser <- read.csv(file.path(RES, "dadras_loo_khat_series.csv")) %>% dplyr::filter(kind == "block_max", trial == SERIES_TRIAL, log10n %in% c(4, 5))
k_lim <- c(min(ser$k) - 0.05, max(ser$k) * 1.05)
khat_panel <- function(l10, tag) {
  s <- ser %>% dplyr::filter(log10n == l10) %>% mutate(flag = k > K_FLAG)
  tmax <- if (l10 == 4) 100 else 1000
  cat(sprintf("k-hat panel: trial %d, n = 1e%g, %d blocks, %d with a flagged observation\n", SERIES_TRIAL, l10, nrow(s), sum(s$flag)))
  ggplot(s, aes(t, k)) +
    geom_point(data = dplyr::filter(s, !flag), colour = HOUSE_PRIMARY, size = 1, alpha = 0.6) +
    geom_point(data = dplyr::filter(s, flag), colour = HOUSE_HIGHLIGHT, shape = 17, size = 3) +
    geom_hline(yintercept = K_FLAG, linetype = 2, colour = HOUSE_HIGHLIGHT, linewidth = 0.6) +
    scale_y_continuous(limits = k_lim, breaks = seq(-1, 8, 1), expand = expansion(mult = 0.02)) +
    scale_x_continuous(breaks = pretty_breaks(5), limits = c(0, tmax), expand = expansion(mult = 0.02)) +
    labs(x = "time", y = expression("Pareto " * hat(italic(k))[italic(i)] * " (block max.)"), tag = tag,
         title = as.expression(substitute("Dadras " * dot(italic(x))[1] * ", " * paste(italic(n) == 10^L, " (", italic(t) <= T, ")"), list(L = l10, T = tmax)))) +
    house_theme()
}
pc1 <- khat_panel(4, "e"); pc2 <- khat_panel(5, "f")
key_3 <- key_panel(data.frame(label = c('paste(hat(italic(k))[italic(i)] <= 0.7)', 'paste(hat(italic(k))[italic(i)] > 0.7, " (flagged)")', 'paste("threshold ", hat(italic(k))[italic(i)] == 0.7)'),
                              colour = c(HOUSE_PRIMARY, HOUSE_HIGHLIGHT, HOUSE_HIGHLIGHT), type = c("point", "point", "dashed"),
                              shape = c(16, 17, NA), parse = TRUE), point_size = 3)

# ------------------------------- row 4: flagged points on the attractor, intercept vs n
fl <- read.csv(file.path(RES, "dadras_loo_flagged.csv")) %>% dplyr::filter(k > K_FLAG, log10n == 5)
set.seed(100); ic <- t(sapply(1:100, function(i) runif(3, -4, 4)))[1, ]       # the benchmark pool, row 1 (dadras_psis_loo.R)
traj <- simulate("dadras", 1e5, Inf, ic, 1)$clean
traj <- as.data.frame(traj[seq(1, nrow(traj), by = 2), ]); names(traj) <- c("x1", "x2", "x3")
n_tr <- length(unique(fl$trial)); n_pts <- nrow(fl)
cat(sprintf("flagged observations panel: %d trials, %d points with k-hat > %.1f at n = 1e5\n", n_tr, n_pts, K_FLAG))
pd1 <- ggplot() +
  geom_path(data = traj, aes(x1, x3), colour = HOUSE_GREY, linewidth = 0.25) +
  geom_point(data = fl, aes(x1, x3), colour = HOUSE_HIGHLIGHT, shape = 17, size = 3) +
  labs(x = expression(italic(x)[1]), y = expression(italic(x)[3]), tag = "g",
       title = as.expression(substitute("Dadras, " * italic(n) == 10^5 * ", " * K * " trials", list(K = n_tr)))) + house_theme()

st <- read.csv(file.path(RES, "dadras_stored_intercept_vs_n.csv")) %>% dplyr::filter(log10n >= 4.6)
pd2 <- ggplot(st, aes(log10n)) +
  geom_ribbon(aes(ymin = int_min, ymax = int_max), fill = alpha(HOUSE_HIGHLIGHT, 0.2)) +
  geom_line(aes(y = int_median), colour = HOUSE_HIGHLIGHT, linewidth = HOUSE_LINE) +
  geom_point(aes(y = int_median), colour = HOUSE_HIGHLIGHT, size = HOUSE_MARKER) +
  geom_line(aes(y = halfwidth_90), colour = HOUSE_PRIMARY, linewidth = HOUSE_LINE, linetype = "longdash") +
  geom_text(aes(y = int_max + 0.0012, label = n_intercept), size = HOUSE_PT$text / .pt, colour = HOUSE_HIGHLIGHT, family = HOUSE_FONT) +
  scale_x_continuous(breaks = seq(4.6, 5, 0.1), labels = c(expression(10^4.6), expression(10^4.7), expression(10^4.8), expression(10^4.9), expression(10^5))) +
  scale_y_continuous(limits = c(0, 0.017), breaks = seq(0, 0.016, 0.004), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(italic(n)), y = "intercept", tag = "h", title = eq_title("Dadras", 1, quote("SNR = 49 dB"))) + house_theme()
key_4 <- key_panel(data.frame(label = c("clean trajectory\n(trial 1)", 'paste("flagged, ", hat(italic(k))[italic(i)] > 0.7)',
                                        "admitted intercept\n(median and range)", 'paste(1.645, " ", italic(sigma) / sqrt(italic(n)))'),
                              colour = c(HOUSE_GREY, HOUSE_HIGHLIGHT, HOUSE_HIGHLIGHT, HOUSE_PRIMARY),
                              type = c("line", "point", "line+point", "dashed"), shape = c(NA, 17, 16, NA),
                              band = c(FALSE, FALSE, TRUE, FALSE), parse = c(FALSE, TRUE, FALSE, TRUE)), point_size = 3, slots = 7)

# ------------------------------------------------------------------ assemble
p <- key_1 + pa2 + pa1 + key_2 + pb1 + pb2 + key_3 + pc1 + pc2 + key_4 + pd1 + pd2 +
  plot_layout(design = "ABC\nDEF\nGHI\nJKL", widths = c(0.48, 1, 1), heights = c(1.1, 1, 1, 1))
save_figure(p, OUT, width_mm = 183, aspect = 1.15, device = cairo_pdf)
copy_to_manuscript(OUT, "imgs")
cat("wrote", OUT, "\n")
