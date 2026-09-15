# Figure for the weak-form check (Referee 1 Points 6-7), house style of
# ../rebuttal-diagnostics/make_figures.R.  Reads results/weak_form_trials.csv.
#   a  Dadras x1-dot: OLS intercept on the true support (median, IQR over
#      trials) against n, SG derivative vs weak form -- the bias of Point 7;
#   b  Dadras x1-dot: fraction of trials in which the spurious intercept
#      survives the 90 % interval rule against n;
#   c  Aizawa x3-dot / Rossler x3-dot: number of spurious terms kept after the
#      Bayesian stage at finite SNR and in the noiseless limit;
#   d  the same cases: residual sd of the true model (median, range) at finite
#      SNR and in the noiseless limit -- the collapse by 4-5 orders of magnitude
#      that makes any deterministic error "significant" happens for both
#      estimators.
# Usage: Rscript make_figures_weak_form.R [primary="weak/raw/m=1w/p=8"]
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(latex2exp) })
args <- commandArgs(trailingOnly = TRUE)
PRIMARY <- if (length(args) >= 1) args[1] else "weak/raw/m=2w/p=8"
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
RES <- file.path(here, "results"); FIG <- file.path(here, "figures"); dir.create(FIG, showWarnings = FALSE)
d <- read.csv(file.path(RES, "weak_form_trials.csv")) %>% filter(arm %in% c("sg", PRIMARY)) %>%
  mutate(method = ifelse(arm == "sg", "Savitzky-Golay derivative", "weak form"),
         snr_lab = ifelse(is.infinite(snr), "SNR = ∞", paste0("SNR = ", snr, " dB")),
         case = paste0(c(dadras = "Dadras", aizawa = "Aizawa", rossler = "Rössler")[system], " ", c(dadras = "x₁", aizawa = "x₃", rossler = "x₃")[system]))
th <- theme(axis.line = element_line(colour = "black"), axis.ticks.length = unit(.2, "cm"),
            panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
            legend.key = element_blank(), legend.position = "top", legend.title = element_blank(),
            text = element_text(size = 15), plot.tag = element_text(face = "bold", size = 18),
            strip.background = element_blank(), strip.text = element_text(size = 14))
COL <- c("Savitzky-Golay derivative" = "#d74d3d", "weak form" = "#505ed5")

dd <- d %>% filter(system == "dadras")
a <- dd %>% group_by(method, log10n) %>% summarise(med = median(ols_int_true_support), lo = quantile(ols_int_true_support, .25), hi = quantile(ols_int_true_support, .75), .groups = "drop") %>%
  ggplot(aes(log10n, med, colour = method, fill = method)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, colour = NA) + geom_line(linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  labs(x = TeX("$\\log_{10} n$"), y = "OLS intercept, true support", tag = "a", title = "Dadras x₁, 49 dB") + th
b <- dd %>% group_by(method, log10n) %>% summarise(frac = mean(kept_intercept), .groups = "drop") %>%
  ggplot(aes(log10n, frac, colour = method)) + geom_line(linewidth = 1) + geom_point(size = 2.5) + scale_colour_manual(values = COL) +
  scale_y_continuous(limits = c(0, 1)) + labs(x = TeX("$\\log_{10} n$"), y = "spurious intercept kept (fraction)", tag = "b", title = "Dadras x₁, 49 dB") + th
da <- d %>% filter(system != "dadras") %>% mutate(snr_lab = factor(snr_lab, levels = unique(snr_lab[order(is.infinite(snr), snr)])))
c <- da %>% group_by(case, snr_lab, method) %>% summarise(med = median(n_spurious_kept), lo = quantile(n_spurious_kept, .25), hi = quantile(n_spurious_kept, .75), .groups = "drop") %>%
  ggplot(aes(snr_lab, med, colour = method)) + geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = .5), size = .7) +
  facet_wrap(~case, scales = "free_x") + scale_colour_manual(values = COL) + labs(x = NULL, y = "spurious terms kept", tag = "c") + th
e <- da %>% group_by(case, snr_lab, method) %>% summarise(med = median(sigma_true_support), lo = min(sigma_true_support), hi = max(sigma_true_support), .groups = "drop") %>%
  ggplot(aes(snr_lab, med, colour = method)) + geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = .5), size = .7) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + facet_wrap(~case, scales = "free_x") + scale_colour_manual(values = COL) +
  labs(x = NULL, y = "residual sd, true support", tag = "d") + th
noleg <- theme(legend.position = "none")
p <- (a + (b + noleg)) / ((c + noleg) + (e + noleg))
ggsave(file.path(FIG, "fig_weak_form_check.pdf"), p, width = 14, height = 11, device = cairo_pdf)
cat("wrote", file.path(FIG, "fig_weak_form_check.pdf"), "primary arm:", PRIMARY, "\n")
