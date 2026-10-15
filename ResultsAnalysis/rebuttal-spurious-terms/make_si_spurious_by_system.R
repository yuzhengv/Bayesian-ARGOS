# Referee 1 Point 3 (R1.3): Supplementary figure, spurious terms per equation
# for every system and method against n (SNR = 49 dB) and against SNR
# (n = 5000); the per-system detail behind block a of the revised Fig. 6.
#   figures/spurious_terms_by_system.pdf  (copied to si-plots/ and reply-to-reviewers/Figures/)
# Canvas 18 x 23.4 in (k = 2.5 route), final width 183 mm, aspect 1.3; 7 rows
# (systems) x 2 columns (sweeps), tags a to n, system name per row, column
# titles once, one legend row below. All 14 panels share the y axis.
# Usage: Rscript make_si_spurious_by_system.R

source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "common.R"))
source(file.path(here, "panel_helpers.R"))
suppressPackageStartupMessages(library(patchwork))
t_start <- Sys.time()

s <- read.csv(file.path(here, "results", "spurious_summary.csv"))
ymax <- max(s$spurious_per_eq[s$sweep == "n" | s$x >= SNR_XMIN]) * 1.05
panels <- list(); tags <- letters
k <- 0
for (i in seq_along(SYSTEMS)) {
    sys <- SYSTEMS[i]
    for (sw in c("n", "snr")) {
        k <- k + 1
        d <- s %>% filter(system == sys, sweep == sw) %>% transmute(x, method, y = spurious_per_eq)
        panels[[k]] <- sweep_panel(d, sw, "spurious_per_eq", ymax,
            title = if (i == 1) c(n = "SNR = 49 dB", snr = "n = 5000")[[sw]] else NULL,
            ylab = if (sw == "n") "Spurious terms\nper equation" else NULL,
            xlab = i == length(SYSTEMS), tag = tags[k], marker = HOUSE_POINT)
    }
}
rows <- lapply(seq_along(SYSTEMS), function(i) {
    row_label(SYSTEM_LABEL[[SYSTEMS[i]]]) + panels[[2 * i - 1]] + panels[[2 * i]] + plot_layout(widths = c(0.05, 1, 1))
})
fig <- Reduce(`/`, rows) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(),
          legend.text = element_text(size = HOUSE_PT$text), legend.key.width = unit(0.6, "in"))
out <- file.path(here, "figures", "spurious_terms_by_system.pdf")
save_figure(fig, out, width_mm = 183, aspect = 1.3)
for (dest in c(file.path(manuscript_root, "si-plots"), file.path(manuscript_root, "reply-to-reviewers", "Figures"))) {
    if (dir.exists(dest)) file.copy(out, file.path(dest, basename(out)), overwrite = TRUE)
}
cat("written:", out, " runtime:", format(Sys.time() - t_start), "\n")
