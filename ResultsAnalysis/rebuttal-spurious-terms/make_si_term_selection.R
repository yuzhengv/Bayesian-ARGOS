# Referee 1 Point 3 (R1.3): Supplementary term-selection figures, one per
# system, every method under both sweeps (the submitted SI showed two pairwise
# comparisons per system).
#   figures/circ_<system>_<sweep>_<method>.pdf   42 circular plots (12 x 12 in
#                                                canvas, text sized for PLOT_MM)
#   figures/term_selection_<system>.tex / .pdf   composites, copied to si-plots/
#                                                and reply-to-reviewers/Figures/
# Layout: 3 rows (Bayesian-ARGOS, ARGOS, SINDy) x 2 columns (n sweep at
# SNR = 49 dB, rings 10^2.8 ... 10^5; SNR sweep at n = 5000, rings as the
# submitted SI figures used them, the noiseless ring labelled with the infinity
# symbol); rings run from the hardest condition outside (smallest n, lowest
# SNR) to the easiest inside; tags a to f row-major, method as row title,
# sweep as column title, one legend line. Three stacked square plots give an
# aspect near 1.5, so with PLOT_MM = 76 the figure is 161 x 243 mm, inside the
# standard's 247 mm page depth; 88 mm plots would make it 280 mm tall.
# Usage: Rscript make_si_term_selection.R

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
source(file.path(root_path, "ResultsAnalysis", "plot-style", "house_style.R"))   # HOUSE_COLOURS
source(file.path(here, "circular_plot_house.R"))
source(file.path(here, "tex_compose.R"))
t_start <- Sys.time()
figdir <- file.path(here, "figures")
PLOT_MM <- 76; ROW_TITLE_MM <- 5; GUTTER_MM <- 4; ROW_GAP_MM <- 3
METHOD_KEY <- c("Bayesian-ARGOS" = "ba", "ARGOS" = "argos", "SINDy" = "sindy")

f <- read.csv(file.path(here, "results", "term_frequencies.csv"))
stats <- list()
for (sys in SYSTEMS) {
    body <- c(); tag <- 0
    x_img <- c(ROW_TITLE_MM, ROW_TITLE_MM + PLOT_MM + GUTTER_MM)
    y0 <- 4.5
    for (sw in c("n", "snr")) {
        body <- c(body, tex_text(x_img[match(sw, c("n", "snr"))] + PLOT_MM / 2, 0,
                                 if (sw == "n") "Against \\textit{n} at SNR = 49 dB" else "Against SNR at \\textit{n} = 5000", anchor = "north"))   # \textit, not $n$: TeX math would load a serif font
    }
    for (i in seq_along(METHODS)) {
        m <- METHODS[i]
        y_img <- y0 + (i - 1) * (PLOT_MM + ROW_GAP_MM)
        body <- c(body, sprintf("\\node[anchor=center, inner sep=0pt, rotate=90] at (%.3f,%.3f) {\\Lab{%s}};",
                                ROW_TITLE_MM / 2 - 0.5, y_img + PLOT_MM / 2, m))
        for (sw in c("n", "snr")) {
            tag <- tag + 1
            rings <- if (sw == "n") rev(RINGS_N) else rev(RINGS_SNR[[sys]])   # innermost first
            labels <- if (sw == "n") ring_label_n(rings) else ring_label_snr(rings)
            heading <- if (sw == "n") expression(italic(n)) else "SNR (dB)"
            file <- file.path(figdir, sprintf("circ_%s_%s_%s.pdf", sys, sw, METHOD_KEY[[m]]))
            d <- f %>% filter(system == sys, sweep == sw, method == m)
            res <- draw_circular_pdf(file, term_freq = d, rings = rings, ring_labels = labels, ring_heading = heading,
                                     hue = HOUSE_COLOURS[[m]], min_freq_for_label = 0.1, final_width_mm = PLOT_MM)
            stats[[length(stats) + 1]] <- res$labels %>% mutate(system = sys, method = m, sweep = sw)
            j <- match(sw, c("n", "snr"))
            body <- c(body,
                      tex_text(x_img[j] - 1.5, y_img - 0.5, letters[tag], anchor = "north west", style = "Tag"),
                      tex_image(x_img[j], y_img, PLOT_MM, basename(file)))
        }
    }
    y_leg <- y0 + 3 * PLOT_MM + 2 * ROW_GAP_MM + 1.5
    width <- ROW_TITLE_MM + 2 * PLOT_MM + GUTTER_MM
    body <- c(body, tex_fill_legend(ROW_TITLE_MM + (2 * PLOT_MM + GUTTER_MM - 121) / 2, y_leg))
    out <- compose_tex(file.path(figdir, sprintf("term_selection_%s.tex", sys)), body)
    for (dest in c(file.path(manuscript_root, "si-plots"), file.path(manuscript_root, "reply-to-reviewers", "Figures"))) {
        if (dir.exists(dest)) file.copy(out, file.path(dest, basename(out)), overwrite = TRUE)
    }
    cat(sprintf("%-10s %s  %s mm\n", sys, basename(out), paste(round(pdf_size_mm(out), 1), collapse = " x ")))
}
stats <- bind_rows(stats)
write.csv(stats, file.path(here, "results", "si_circular_label_stats.csv"), row.names = FALSE)
cat("runtime:", format(Sys.time() - t_start), "\n")
