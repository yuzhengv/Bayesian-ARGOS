# Referee 1 Point 3 (R1.3): revised Fig. 6, shrinkage and spurious terms.
#   figures/fig6_block_a.pdf        block a: 2 x 3 aggregate panels (ggplot2,
#                                   canvas 18 x 10.8 in = k 2.5 of 183 x 110 mm)
#   figures/circ_lorenz_n_<m>.pdf   block b: Lorenz, n sweep, one circular plot
#                                   per method (12 x 12 in canvas, shown at 59 mm)
#   figures/shrinkage_spurious_terms.tex / .pdf   the composite, 183 mm wide,
#                                   copied to imgs/ and reply-to-reviewers/Figures/
# Block a: rows = sweep (n at SNR = 49 dB; SNR at n = 5000), columns = spurious
# terms per equation, true terms missed per equation, mean selection frequency
# of the spurious terms; each curve is the mean over the seven systems, the
# ribbon the minimum to maximum over the systems; tags a to f, one legend row.
# Block b: rings 10^5 (inner) to 10^2.8 (outer); tags g to i.
# Usage: Rscript make_fig6.R

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
source(file.path(here, "panel_helpers.R"))
source(file.path(here, "circular_plot_house.R"))
source(file.path(here, "tex_compose.R"))
suppressPackageStartupMessages(library(patchwork))
t_start <- Sys.time()
figdir <- file.path(here, "figures")

# ---- block a ------------------------------------------------------------------
s <- read.csv(file.path(here, "results", "spurious_summary.csv"))
agg <- s %>%
    pivot_longer(c(spurious_per_eq, missed_per_eq, mean_freq_spurious), names_to = "quantity", values_to = "v") %>%
    group_by(sweep, quantity, method, x) %>%
    summarise(y = mean(v), ymin = min(v), ymax = max(v), n_systems = n(), .groups = "drop")
stopifnot(all(agg$n_systems == 7))
write.csv(agg %>% select(-n_systems), file.path(here, "results", "spurious_summary_aggregate.csv"), row.names = FALSE)
shown <- agg %>% filter(sweep == "n" | x >= SNR_XMIN)
ymax <- c(spurious_per_eq = max(shown$ymax[shown$quantity == "spurious_per_eq"]) * 1.05,
          missed_per_eq = max(shown$ymax[shown$quantity == "missed_per_eq"]) * 1.05,
          mean_freq_spurious = 1)
quantities <- c("spurious_per_eq", "missed_per_eq", "mean_freq_spurious")
panels <- list(); k <- 0
for (sw in c("n", "snr")) {
    for (q in quantities) {
        k <- k + 1
        d <- agg %>% filter(sweep == sw, quantity == q) %>% select(x, method, y, ymin, ymax)
        panels[[k]] <- sweep_panel(d, sw, q, ymax[[q]], title = if (sw == "n") QUANTITY_TITLE[[q]] else NULL,
                                   ylab = QUANTITY_YLAB[[q]], xlab = TRUE, tag = letters[k], marker = HOUSE_POINT)
    }
}
row1 <- row_label("SNR = 49 dB") + panels[[1]] + panels[[2]] + panels[[3]] + plot_layout(widths = c(0.05, 1, 1, 1))
row2 <- row_label("n = 5000") + panels[[4]] + panels[[5]] + panels[[6]] + plot_layout(widths = c(0.05, 1, 1, 1))
fig_a <- (row1 / row2) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(),
          legend.text = element_text(size = HOUSE_PT$text), legend.key.width = unit(0.6, "in"))
block_a <- file.path(figdir, "fig6_block_a.pdf")
save_figure(fig_a, block_a, width_mm = 183, aspect = 0.60)

# ---- block b ------------------------------------------------------------------
f <- read.csv(file.path(here, "results", "term_frequencies.csv"))
rings <- rev(RINGS_N)   # innermost first
circ_width <- 59
circ_files <- c(); label_stats <- list()
for (m in METHODS) {
    key <- c("Bayesian-ARGOS" = "ba", "ARGOS" = "argos", "SINDy" = "sindy")[[m]]
    file <- file.path(figdir, sprintf("circ_lorenz_n_%s.pdf", key))
    d <- f %>% filter(system == "lorenz", sweep == "n", method == m)
    res <- draw_circular_pdf(file, term_freq = d, rings = rings, ring_labels = ring_label_n(rings),
                             ring_heading = expression(italic(n)), hue = HOUSE_COLOURS[[m]],
                             min_freq_for_label = 0.1, final_width_mm = circ_width)
    circ_files[m] <- basename(file)
    label_stats[[m]] <- res$labels %>% mutate(method = m)
}
label_stats <- bind_rows(label_stats)
write.csv(label_stats, file.path(here, "results", "fig6_circular_label_stats.csv"), row.names = FALSE)

# ---- composite ----------------------------------------------------------------
sz <- pdf_size_mm(block_a); h_a <- 183 * sz[2] / sz[1]
gap <- (183 - 3 * circ_width) / 2
xs <- c(0, circ_width + gap, 2 * (circ_width + gap))
y_head <- h_a + 1.5; y_title <- y_head + 4.5; y_img <- y_title + 3.6; y_leg <- y_img + circ_width + 2
body <- c(
    tex_image(0, 0, 183, basename(block_a)),
    tex_text(91.5, y_head, "Lorenz, SNR = 49 dB, 100 trials per ring", anchor = "north"),
    unlist(lapply(seq_along(METHODS), function(i) c(
        tex_text(xs[i], y_title - 0.4, letters[6 + i], anchor = "north west", style = "Tag"),
        tex_text(xs[i] + circ_width / 2, y_title, METHODS[i], anchor = "north"),
        tex_image(xs[i], y_img, circ_width, circ_files[[METHODS[i]]])
    ))),
    tex_fill_legend(31, y_leg)
)
out <- compose_tex(file.path(figdir, "shrinkage_spurious_terms.tex"), body)
for (dest in c(file.path(manuscript_root, "imgs"), file.path(manuscript_root, "reply-to-reviewers", "Figures"))) {
    if (dir.exists(dest)) file.copy(out, file.path(dest, basename(out)), overwrite = TRUE)
}
cat("written:", out, " size (mm):", paste(round(pdf_size_mm(out), 1), collapse = " x "),
    " runtime:", format(Sys.time() - t_start), "\n")
