# Referee 1 Point 3 (R1.3): the circular term-selection plot in the house style.
#
# A copy of plot_freq_of_identified_terms_circular_equations_as_sectors_refined
# (MethodsEvaluation/evaluation-functions-for-plots/
# functions_for_analysis_of_identified_terms_vary_n_enhanced.R, from line 2322),
# modified here only; the package copy is untouched. The layout is the same:
# one sector per equation, one ring per sweep value, one bar per selected term
# with the bar height equal to its selection frequency over the 100 trials,
# dark fill = true term, light tint of the same hue = spurious term. Changes:
#   - the data come from results/term_frequencies.csv (one frequency per
#     equation, ring and term) instead of being recounted from the summary
#     data frame, so that the plot and the aggregate panel share one source
#   - sans-serif text (the pdf device's Helvetica) at sizes computed from the
#     final printed width: term labels 5 pt, ring labels 6.5 pt, headings 7 pt
#   - rings labelled with the sweep values in a vertical column at the top
#     (10^k, or SNR in dB with the noiseless ring labelled by the infinity
#     symbol); the gap before the first sector is widened to hold the column
#   - sector labels are the derivatives (x-dot with subscript) on the grey
#     outer arc; no radial axis and no grid beyond the ring boundaries
#   - rings are drawn with a frequency of 1 reaching bar_max of the ring
#     thickness so that the labels have room; a label is drawn for every term
#     whose frequency reaches min_freq_for_label, and when more labels are
#     requested in a cell than fit at the minimum size, the labelled bars are
#     spread over the cell and the most frequent terms keep their labels (the
#     function reports what was dropped)
#   - one hue per method, passed in; no bold, no title, no legend (the
#     composite carries them)

suppressPackageStartupMessages(library(circlize))

# "$x_{1}^{2}x_{3}$" -> plotmath expression x[1]^2 * x[3]; "$C$" -> C;
# "$sin(x_{2})$" -> sin(x[2]). The substitution chain of the original function.
latex_term_to_plotmath <- function(label) {
    m <- gsub("^\\$|\\$$", "", label)
    if (m == "(Intercept)") m <- "C"
    m <- gsub("_\\{([^}]+)\\}", "[\\1]", m)
    m <- gsub("_([0-9a-zA-Z])", "[\\1]", m)
    m <- gsub("\\^\\{([^}]+)\\}", "^{\\1}", m)
    m <- gsub("\\\\([a-zA-Z]+)\\(([^)]+)\\)", "\\1(\\2)", m)
    m <- gsub("([a-zA-Z][\\[0-9\\]]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9\\]]+)", "\\1*\\2", m, perl = TRUE)
    m <- gsub("(\\][0-9]*(?:\\^\\{?[0-9a-zA-Z]+\\}?)?)([a-zA-Z][\\[0-9])", "\\1*\\2", m)
    m <- gsub("(\\^\\{?[0-9a-zA-Z]+\\}?)([a-zA-Z][\\[0-9])", "\\1*\\2", m)
    m <- gsub("\\*\\*+", "*", m)
    m <- gsub("\\*$", "", m)
    tryCatch(parse(text = m), error = function(e) label)   # an expression, as circos.text indexes labels[i]
}

# term_freq: data frame with equation (1..3), x (ring value), term (LaTeX
#   label), frequency (0..1), is_truth (logical); rows with frequency 0 are
#   ignored.
# rings: ring values from the innermost ring outward; ring_labels: a list of
#   expressions or strings in the same order; ring_heading: expression drawn
#   above the column of ring labels.
# Sizes are final sizes in points at final_width_mm; the canvas is canvas_in
#   square with the given pointsize, so cex = pt * (canvas_in * 25.4 /
#   final_width_mm) / pointsize.
circular_term_plot_house <- function(
    term_freq, rings, ring_labels, ring_heading = NULL, hue = "#505ed5",
    min_freq_for_label = 0.1, final_width_mm = 60, canvas_in = 12, pointsize = 12,
    label_pt = 5, text_pt = 7, tick_pt = 6.5, bar_max = 0.35, inner_hole = 0.14,
    eq_track = 0.07, gap_deg = c(4, 4, 30), tint_alpha = 0.4, ring_border = "grey80",
    canvas_lim = 1.0, sector_labels = list(expression(dot(x)[1]), expression(dot(x)[2]), expression(dot(x)[3]))
) {
    scale <- canvas_in * 25.4 / final_width_mm
    cex_label <- label_pt * scale / pointsize
    cex_text <- text_pt * scale / pointsize
    cex_tick <- tick_pt * scale / pointsize
    n_rings <- length(rings)
    n_eq <- length(sector_labels)
    sectors <- paste0("Eq", seq_len(n_eq))
    track_margin <- 0.003
    ring_h <- (1 - eq_track - inner_hole - 2 * track_margin * (n_rings + 1)) / n_rings
    y_max <- 1 / bar_max
    R_pt <- canvas_in * 72 / 2 / canvas_lim       # device points per unit radius
    sector_rad <- (360 - sum(rep_len(gap_deg, n_eq))) / n_eq * pi / 180
    label_h_pt <- label_pt * scale * 1.15           # tangential footprint of a radial label (sub- and superscripts)
    label_log <- list()
    yplot <- vector("list", n_rings)

    term_freq <- term_freq[term_freq$frequency > 0, ]
    circos.clear()
    circos.par(start.degree = 90 - gap_deg[length(gap_deg)] / 2, gap.degree = gap_deg,
               track.margin = c(track_margin, track_margin), cell.padding = c(0, 0, 0, 0),
               canvas.xlim = c(-canvas_lim, canvas_lim), canvas.ylim = c(-canvas_lim, canvas_lim),
               points.overflow.warning = FALSE)
    circos.initialize(factors = sectors, xlim = c(0, 1))

    # outer arc with the equation labels
    # light band with black labels: the dot accent of x-dot rises above the
    # band and must stay visible on the white page
    circos.track(factors = sectors, ylim = c(0, 1), track.height = eq_track, bg.col = "grey88",
                 bg.border = "white", panel.fun = function(x, y) {
        k <- CELL_META$sector.numeric.index
        circos.text(0.5, 0.45, labels = sector_labels[[k]], cex = cex_text, col = "black",
                    facing = "inside", niceFacing = TRUE)
    })

    # rings, created from the outermost inward (circlize adds tracks inward)
    for (j in rev(seq_len(n_rings))) {
        ring_val <- rings[j]
        circos.track(factors = sectors, ylim = c(0, y_max), track.height = ring_h, bg.col = NA,
                     bg.border = ring_border, bg.lwd = 0.6, panel.fun = function(x, y) {
            eq <- CELL_META$sector.numeric.index
            if (eq == 1) yplot[[j]] <<- CELL_META$yplot
            d <- term_freq[term_freq$equation == eq & term_freq$x == ring_val, ]
            if (nrow(d) == 0) return(invisible())
            d <- d[order(!d$is_truth, -d$frequency, d$term), ]
            n <- nrow(d)
            w <- min(0.84 / (n * 1.2), 0.1)
            pitch <- w * 1.2
            # tangential room for a label at the label radius, as a fraction of the sector
            r_lab <- CELL_META$yplot[1] + bar_max * (CELL_META$yplot[2] - CELL_META$yplot[1]) * 0.5
            pitch_needed <- label_h_pt / (sector_rad * r_lab * R_pt)
            # few bars: spread them so that every label has room
            if (n * pitch_needed <= 0.84 && pitch < pitch_needed) pitch <- pitch_needed
            slots <- 0.08 + w / 2 + (seq_len(n) - 1) * pitch
            pitch_pt <- pitch * sector_rad * r_lab * R_pt
            m <- max(1L, ceiling(label_h_pt / pitch_pt))
            cand <- which(d$frequency >= min_freq_for_label)
            k_fit <- floor((n - 1) / m) + 1
            k <- min(length(cand), k_fit)
            lab_rows <- head(cand, k)                      # the most frequent (true first) keep their labels
            if (k < n) {
                lab_slots <- 1 + (seq_len(k) - 1) * m
                other_slots <- setdiff(seq_len(n), lab_slots)
                order_rows <- integer(n)
                order_rows[lab_slots] <- lab_rows
                order_rows[other_slots] <- setdiff(seq_len(n), lab_rows)
                d <- d[order_rows, ]
                lab_rows <- lab_slots
            }
            label_log[[paste(eq, ring_val)]] <<- data.frame(equation = eq, ring = ring_val, n_terms = n,
                n_candidates = length(cand), n_labelled = k, min_freq_labelled = if (k > 0) min(d$frequency[lab_rows]) else NA)
            for (i in seq_len(n)) {
                col <- if (d$is_truth[i]) hue else adjustcolor(hue, alpha.f = tint_alpha)
                circos.rect(slots[i] - w / 2, 0, slots[i] + w / 2, d$frequency[i], col = col, border = col, lwd = 0.5)
                if (i %in% lab_rows) {
                    circos.text(slots[i], d$frequency[i] + 0.04 * y_max, labels = latex_term_to_plotmath(d$term[i]),
                                cex = cex_label, col = "black", facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                }
            }
        })
    }
    # ring labels in the gap at the top, and the heading
    for (j in seq_len(n_rings)) {
        r_mid <- mean(yplot[[j]])
        text(0, r_mid, labels = ring_labels[[j]], cex = cex_tick, col = "black", adj = c(0.5, 0.5))
    }
    if (!is.null(ring_heading)) {
        text(0, 1 - eq_track / 2 - track_margin, labels = ring_heading, cex = cex_text, col = "black", adj = c(0.5, 0.5))
    }
    circos.clear()
    invisible(list(labels = do.call(rbind, label_log), cex = c(label = cex_label, text = cex_text, tick = cex_tick),
                   label_pt_final = label_pt, ring_h = ring_h))
}

# Draws one circular plot to a vector PDF (canvas_in square, Helvetica).
draw_circular_pdf <- function(file, ..., canvas_in = 12, pointsize = 12) {
    pdf(file, width = canvas_in, height = canvas_in, family = "Helvetica", pointsize = pointsize, useDingbats = FALSE)
    on.exit(dev.off(), add = TRUE)
    par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA, family = "sans")
    circular_term_plot_house(..., canvas_in = canvas_in, pointsize = pointsize)
}

ring_label_n <- function(v) lapply(v, function(k) parse(text = paste0("10^", k)))
ring_label_snr <- function(v) lapply(v, function(s) if (is.infinite(s)) expression(infinity) else as.character(s))
