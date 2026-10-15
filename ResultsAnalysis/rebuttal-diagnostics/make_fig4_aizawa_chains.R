# Revised Fig. 4 of the manuscript (plotting standard block M4): mechanism I,
# near-collinear chains of terms in the candidate library, Aizawa x3-dot.
#   row a (symptom)    term-selection frequency of the final models (fraction of
#                      the 100 stored benchmark trials, results/aizawa_eq3_term_
#                      frequencies.csv) as one bar per candidate term and one
#                      facet per condition: n = 10^3.7, 10^4.4, 10^5 at 49 dB
#                      (left) and SNR = 50, 55, 61 dB and infinity at n = 5000
#                      (right); terms with a frequency of at least 0.1 at any
#                      shown condition; the eight true terms in their colours
#                      (intercept black), every other term grey; the chain
#                      neighbour x2^2 x3^2 annotated; the key at the left names
#                      the true terms in the order of the bars (no tick labels)
#   row b (diagnostic) median VIF (20 trials) of every term of the design entering
#                      the second adaptive-lasso pass against n and against SNR
#                      (results/exp1_vif_long.csv), true terms coloured, others
#                      grey, dashed VIF = 10 line, both panels on one log scale
#   row c (mechanism)  condition number of the four nested designs, median and
#                      interquartile range (results/exp1_kappa.csv), one log scale
# One vertical key per row at the left of the two compared panels (standard
# 4.7.2).  Canvas 18 x 19.8 in at k = 2.5 (final width 183 mm, aspect 1.1),
# cairo_pdf.  Writes figures/diagnostics_aizawa_chains.pdf and copies it to
# <manuscript>/imgs and <manuscript>/reply-to-reviewers/Figures.
# Usage: Rscript make_fig4_aizawa_chains.R        (seconds)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(patchwork) })
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(here, "..", "plot-style", "house_style.R"))
source(file.path(here, "fig_helpers.R"))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)
OUT <- file.path(FIG, "diagnostics_aizawa_chains.pdf")
TITLE_N <- "SNR = 49 dB"; TITLE_SNR <- expression(italic(n) == 5000)
CHAIN <- "x2^2x3^2"

# ------------------------------------------------------------- row a: symptom
fa <- read.csv(file.path(RES, "aizawa_eq3_term_frequencies.csv"), stringsAsFactors = FALSE)
N_SHOW <- c(3.7, 4.4, 5); SNR_SHOW <- c(50, 55, 61, Inf)
fn <- fa %>% filter(sweep == "n", value %in% N_SHOW) %>% mutate(cond = factor(value, levels = N_SHOW))
fs <- fa %>% filter(sweep == "snr", value %in% SNR_SHOW) %>% mutate(cond = factor(value, levels = SNR_SHOW))
label_terms <- character(0)      # the key names the eight true terms (same order as the bars); no room for 45-54 tick labels
chain_label <- 'paste(italic(x)[2]^2*italic(x)[3]^2, " replaces ", italic(x)[2]^2*italic(x)[3])'
ann_n <- data.frame(term = CHAIN, cond = factor(5, levels = N_SHOW), label = chain_label, x_off = 4, hjust = 0,
                    y_text = 0.72, y_arrow = fn$frequency[fn$term == CHAIN & fn$value == 5] + 0.03)
ann_s <- data.frame(term = CHAIN, cond = factor(Inf, levels = SNR_SHOW), label = chain_label, x_off = 4, hjust = 0,
                    y_text = 0.74, y_arrow = fs$frequency[fs$term == CHAIN & is.infinite(fs$value)] + 0.03)
pa1 <- skyline_panel(fn, AIZAWA_TRUTH, AIZAWA_TRUTH_COL, label_terms,
                     cond_labels = c("italic(n) == 10^3.7", "italic(n) == 10^4.4", "italic(n) == 10^5"), order_by = 5,
                     annotate_df = ann_n) + labs(title = TITLE_N, tag = "a")
pa2 <- skyline_panel(fs, AIZAWA_TRUTH, AIZAWA_TRUTH_COL, label_terms,
                     cond_labels = c('paste("SNR = 50 dB")', 'paste("SNR = 55 dB")', 'paste("SNR = 61 dB")', 'paste("SNR = ", infinity)'),
                     order_by = Inf, annotate_df = ann_s) + labs(title = TITLE_SNR, tag = "b")
key_a <- key_panel(data.frame(label = c(term_label(AIZAWA_TRUTH), "intercept", "other candidate term"),
                              colour = c(unname(AIZAWA_TRUTH_COL), "black", HOUSE_GREY), type = "fill",
                              parse = c(rep(TRUE, 7), FALSE, FALSE)),
                   title = 'paste("Aizawa ", dot(italic(x))[3])', title_parse = TRUE)

# ---------------------------------------------------------- row b: diagnostic
vif <- read.csv(file.path(RES, "exp1_vif_long.csv"), stringsAsFactors = FALSE)
med <- vif %>% filter(design == "trimmed") %>% group_by(sweep, term, log10n, snr) %>%
  summarise(vif = median(vif, na.rm = TRUE), .groups = "drop") %>%
  mutate(term = to_r_name(term), is_truth = term %in% AIZAWA_TRUTH)
VIF_LIM <- c(1, 1e16)
vif_panel <- function(sw) {
  d <- med %>% filter(sweep == sw) %>% mutate(x = if (sw == "snr") snr_x(snr) else log10n)
  if (sw == "snr") d <- d %>% filter(snr >= 6)                     # the axis starts at 6 dB
  fin <- if (sw == "snr") d %>% filter(is.finite(snr)) else d
  others <- fin %>% filter(!is_truth); truth <- d %>% filter(is_truth) %>% mutate(term = factor(term, levels = AIZAWA_TRUTH))
  p <- ggplot() +
    geom_line(data = others, aes(x, vif, group = term), colour = HOUSE_GREY, linewidth = HOUSE_CONTEXT_LINE) +
    geom_line(data = truth %>% filter(x < INF_X), aes(x, vif, colour = term), linewidth = HOUSE_LINE) +
    geom_point(data = truth, aes(x, vif, colour = term, shape = term), size = HOUSE_MARKER) +
    geom_hline(yintercept = 10, linetype = 2, linewidth = 0.44) +
    scale_colour_manual(values = AIZAWA_TRUTH_COL, guide = "none") +
    scale_shape_manual(values = AIZAWA_TRUTH_SHAPE, guide = "none") +
    scale_y_log10_powers("VIF", limits = VIF_LIM, step = 4) + house_theme()
  if (sw == "snr") {
    p <- p + geom_point(data = d %>% filter(!is_truth, is.infinite(snr)), aes(x, vif), colour = HOUSE_GREY, size = 1.2)
    snr_axis_log(p, VIF_LIM)
  } else p + scale_x_log_n()
}
pb1 <- vif_panel("n") + labs(title = TITLE_N, tag = "c")
pb2 <- vif_panel("snr") + labs(title = TITLE_SNR, tag = "d")
key_b <- key_panel(data.frame(label = c(term_label(AIZAWA_TRUTH), "other candidate term"),
                              colour = c(unname(AIZAWA_TRUTH_COL), HOUSE_GREY), type = c(rep("line+point", 7), "line"),
                              shape = c(unname(AIZAWA_TRUTH_SHAPE), 16), parse = c(rep(TRUE, 7), FALSE)))

# ----------------------------------------------------------- row c: mechanism
kappa <- read.csv(file.path(RES, "exp1_kappa.csv"), stringsAsFactors = FALSE)
DES <- c(full = "full candidate\nlibrary", trimmed = "candidate library\nafter the degree trim", selected = "screened support", truth = "true terms only")
DES_COL <- c(full = HOUSE_HIGHLIGHT, trimmed = HOUSE_EXTRA[1], selected = HOUSE_COLOURS[["SINDy"]], truth = HOUSE_PRIMARY)
DES_SHAPE <- c(full = 15, trimmed = 17, selected = 16, truth = 18)
ks <- kappa %>% group_by(sweep, design, log10n, snr) %>%
  summarise(m = median(kappa), lo = quantile(kappa, .25), hi = quantile(kappa, .75), .groups = "drop") %>%
  mutate(design = factor(design, levels = names(DES)))
KAPPA_LIM <- c(1, 1e13)
kappa_panel <- function(sw) {
  d <- ks %>% filter(sweep == sw) %>% mutate(x = if (sw == "snr") snr_x(snr) else log10n)
  if (sw == "snr") d <- d %>% filter(snr >= 6)
  fin <- d %>% filter(x < INF_X)
  p <- ggplot(d, aes(x, m, colour = design, fill = design, shape = design)) +
    geom_ribbon(data = fin, aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(data = fin, linewidth = HOUSE_LINE) + geom_point(size = HOUSE_MARKER) +
    scale_colour_manual(values = DES_COL, guide = "none") + scale_fill_manual(values = DES_COL, guide = "none") +
    scale_shape_manual(values = DES_SHAPE, guide = "none") +
    scale_y_log10_powers(expression("condition number " * italic(kappa)), limits = KAPPA_LIM, step = 2) + house_theme()
  if (sw == "snr") {
    p <- p + geom_linerange(data = d %>% filter(x == INF_X), aes(ymin = lo, ymax = hi), linewidth = 0.6, alpha = 0.6)
    snr_axis_log(p, KAPPA_LIM)
  } else p + scale_x_log_n()
}
pc1 <- kappa_panel("n") + labs(title = TITLE_N, tag = "e")
pc2 <- kappa_panel("snr") + labs(title = TITLE_SNR, tag = "f")
key_c <- key_panel(data.frame(label = unname(DES), colour = unname(DES_COL), type = "line+point", shape = unname(DES_SHAPE), band = TRUE), slots = 9)

# ------------------------------------------------------------------ assemble
p <- key_a + pa1 + pa2 + key_b + pb1 + pb2 + key_c + pc1 + pc2 +
  plot_layout(design = "ABC\nDEF\nGHI", widths = c(0.48, 1, 1), heights = c(1.35, 1, 1))
save_figure(p, OUT, width_mm = 183, aspect = 1.1, device = cairo_pdf)
copy_to_manuscript(OUT, "imgs")
cat("wrote", OUT, "\n")
