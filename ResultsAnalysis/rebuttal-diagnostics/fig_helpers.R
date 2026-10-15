# Small helpers shared by the figure scripts of this folder (fig1_insets.R,
# make_fig4_aizawa_chains.R, make_fig5_derivative_error.R,
# make_si_residual_statistics.R): term labels in plotmath, the fixed colours of
# the seven true terms of the Aizawa x3-dot equation (standard 4.5.4: the same
# colour in every panel of a figure), vertical keys drawn as panels (standard
# 4.7.2) and the term-selection "skyline" panels (decision D6).  The theme,
# palette, axes and save routine come from ../plot-style/house_style.R, which
# every figure script sources first.

# "(Intercept)" / "x1^2x3" (the naming of build_design_matrix and of the stored
# models) or "x1^2 x3" (the sklearn-style naming of the results CSVs) -> plotmath
term_label <- function(term) {
  sapply(term, function(t) {
    if (t == "(Intercept)") return("intercept")
    toks <- regmatches(t, gregexpr("x[0-9](\\^[0-9])?", t))[[1]]
    paste(sapply(toks, function(k) {
      v <- sub("\\^.*", "", k); p <- if (grepl("\\^", k)) sub(".*\\^", "", k) else NULL
      paste0("italic(x)[", substring(v, 2), "]", if (!is.null(p)) paste0("^", p) else "")
    }), collapse = "*")
  }, USE.NAMES = FALSE)
}
term_expr <- function(term) lapply(term_label(term), function(s) parse(text = s)[[1]])
# "x1^2 x3" -> "x1^2x3"
to_r_name <- function(term) gsub(" ", "", term)

# the seven non-constant true terms of Aizawa x3-dot, in the order of the keys,
# with their fixed colours (HOUSE_SERIES: blue, red, green, orange, pink, teal,
# brown) and a marker shape each (standard 4.5.8: colour doubled by a shape)
AIZAWA_TRUTH <- c("x3", "x1^2", "x2^2", "x1^2x3", "x2^2x3", "x1^3x3", "x3^3")
AIZAWA_TRUTH_COL <- setNames(HOUSE_SERIES[seq_along(AIZAWA_TRUTH)], AIZAWA_TRUTH)
AIZAWA_TRUTH_SHAPE <- setNames(c(16, 15, 17, 18, 1, 0, 2), AIZAWA_TRUTH)

# panel titles: system, equation and condition (standard 4.8.2)
eq_title <- function(system, eq, condition = NULL) {
  base <- substitute(S * " " * dot(italic(x))[E], list(S = system, E = eq))
  if (is.null(condition)) return(as.expression(base))
  as.expression(substitute(B * ", " * C, list(B = base, C = condition)))
}
snr_cond <- function(snr) if (is.infinite(snr)) quote("SNR = " * infinity) else substitute("SNR = " * S * " dB", list(S = snr))

# snr_axis() for a log y axis with limits lim: the separators keep the size
# they have on a 0-1 axis (house_style.R, sep_y_scale).
snr_axis_log <- function(p, lim, ...) {
  R <- log10(lim[2] / lim[1])
  snr_axis(p, y_break = lim[1], extra_y = 0.05 * R, ylog = TRUE, sep_y_scale = R, ...)
}

# A vertical key drawn as its own panel (standard 4.7.2: once per row, at the
# left of the row, no background, no frame, text 7 pt final).  One entry per
# row of `entries`: label (a plotmath string when parse is TRUE), colour, type
# ("line", "point", "line+point", "fill", "dashed", "dotted"), optional shape,
# linetype and band (TRUE draws the interquartile ribbon behind the line).  The
# entries are centred vertically on a grid of `slots` rows so that the line
# spacing is the same whatever the number of entries.
key_panel <- function(entries, title = NULL, title_parse = FALSE, slots = 13, x_sym = c(0.02, 0.20), x_text = 0.27,
                      text_pt = HOUSE_PT$text, point_size = HOUSE_MARKER, lineheight = 0.9) {
  stopifnot(is.data.frame(entries))
  n <- nrow(entries)
  if (is.null(entries$parse)) entries$parse <- FALSE
  if (is.null(entries$shape)) entries$shape <- 16
  if (is.null(entries$band)) entries$band <- FALSE
  if (is.null(entries$alpha)) entries$alpha <- 0.2
  n_title <- if (is.null(title)) 0 else 1
  top <- slots / 2 + (n + n_title) / 2
  entries$y <- top - n_title - seq_len(n) + 0.5
  xm <- mean(x_sym)
  p <- ggplot() + theme_void() + theme(text = element_text(family = HOUSE_FONT), plot.margin = margin(0, 0, 0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, slots), expand = FALSE, clip = "off")
  rb <- entries[entries$band, ]
  if (nrow(rb)) p <- p + geom_rect(data = rb, aes(xmin = x_sym[1], xmax = x_sym[2], ymin = y - 0.3, ymax = y + 0.3), fill = alpha(rb$colour, rb$alpha))
  fl <- entries[entries$type == "fill", ]
  if (nrow(fl)) p <- p + geom_rect(data = fl, aes(xmin = x_sym[1] + 0.03, xmax = x_sym[2] - 0.03, ymin = y - 0.32, ymax = y + 0.32), fill = fl$colour)
  ln <- entries[entries$type %in% c("line", "line+point", "dashed", "dotted"), ]
  if (nrow(ln)) p <- p + geom_segment(data = ln, aes(x = x_sym[1], xend = x_sym[2], y = y, yend = y), colour = ln$colour,
                                       linetype = ifelse(ln$type == "dashed", "dashed", ifelse(ln$type == "dotted", "dotted", "solid")),
                                       linewidth = ifelse(ln$type %in% c("dashed", "dotted"), 0.6, HOUSE_LINE))
  pt <- entries[entries$type %in% c("point", "line+point"), ]
  if (nrow(pt)) p <- p + geom_point(data = pt, aes(x = xm, y = y), colour = pt$colour, shape = pt$shape, size = point_size)
  tx <- entries[!entries$parse, ]
  if (nrow(tx)) p <- p + geom_text(data = tx, aes(x = x_text, y = y, label = label), hjust = 0, vjust = 0.5, size = text_pt / .pt, family = HOUSE_FONT, lineheight = lineheight)
  tp <- entries[entries$parse, ]
  if (nrow(tp)) p <- p + geom_text(data = tp, aes(x = x_text, y = y, label = label), hjust = 0, vjust = 0.5, size = text_pt / .pt, family = HOUSE_FONT, parse = TRUE)
  if (!is.null(title)) p <- p + annotate("text", x = x_sym[1], y = top - 0.5, label = title, hjust = 0, vjust = 0.5, size = text_pt / .pt,
                                         family = HOUSE_FONT, fontface = "bold", parse = title_parse)
  p
}

# The term-selection "skyline": one bar per candidate term, one facet per
# condition stacked vertically (decision D6), true terms in colour, the other
# terms grey, tick labels only for the named terms.
# freq: data frame with term, frequency, cond (factor, in display order);
# truth_col: named colours of the true terms; the intercept takes intercept_col
# and goes first (intercept_first) or last; label_terms: the terms that receive
# a tick label; order_by: the condition (a level of cond) whose frequencies
# order the grey terms (decreasing); annotate_df: optional rows (term, cond,
# label, x_off, y_text, y_arrow, hjust) drawn as an arrow and a parsed label.
skyline_panel <- function(freq, truth_terms, truth_col, label_terms, cond_labels, order_by, min_freq = 0.1,
                          intercept_col = "black", intercept_first = TRUE, ylab = "selection frequency",
                          annotate_df = NULL, bar_width = 0.8, label_angle = 90) {
  shown <- unique(freq$term[freq$frequency >= min_freq])
  freq <- freq[freq$term %in% shown, ]
  last <- freq[as.character(freq$cond) == as.character(order_by), ]
  others <- setdiff(shown, c("(Intercept)", truth_terms))
  if (length(others)) others <- others[order(-vapply(others, function(t) { v <- last$frequency[last$term == t]; if (length(v)) v else 0 }, numeric(1)), others)]
  int <- intersect("(Intercept)", shown)
  order_terms <- if (intercept_first) c(int, intersect(truth_terms, shown), others) else c(intersect(truth_terms, shown), int, others)
  freq$index <- match(freq$term, order_terms)
  freq$colour <- ifelse(freq$term == "(Intercept)", intercept_col, ifelse(freq$term %in% truth_terms, truth_col[freq$term], HOUSE_GREY))
  lab_idx <- match(intersect(order_terms, label_terms), order_terms)
  lab <- order_terms[lab_idx]
  labels_df <- data.frame(cond = factor(levels(freq$cond), levels = levels(freq$cond)), label = cond_labels)
  p <- ggplot(freq, aes(index, frequency)) +
    geom_col(fill = freq$colour, width = bar_width) +
    geom_text(data = labels_df, aes(x = Inf, y = 1.02, label = label), parse = TRUE, hjust = 1, vjust = 1,
              size = HOUSE_PT$text / .pt, family = HOUSE_FONT) +
    facet_grid(cond ~ ., axes = "all", axis.labels = "margins") +
    scale_x_continuous(breaks = lab_idx, labels = term_expr(lab), expand = expansion(add = c(0.7, 0.7))) +
    scale_y_continuous(name = ylab, limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = expansion(mult = 0), oob = oob_keep) +
    labs(x = "candidate term") + house_theme() +
    theme(strip.text = element_blank(), strip.background = element_blank(), panel.spacing.y = unit(0.12, "in"),
          axis.text.x = element_text(angle = label_angle, hjust = if (label_angle == 90) 1 else 0.5, vjust = if (label_angle == 90) 0.5 else 1))
  if (!is.null(annotate_df)) {
    annotate_df$index <- match(annotate_df$term, order_terms)
    p <- p + geom_segment(data = annotate_df, aes(x = index + x_off, xend = index, y = y_text - 0.03, yend = y_arrow), linewidth = 0.44,
                          arrow = arrow(length = unit(0.08, "in"), type = "closed"), colour = "black", lineend = "round", linejoin = "mitre") +
      geom_text(data = annotate_df, aes(x = index + x_off, y = y_text, label = label), parse = TRUE, hjust = annotate_df$hjust, vjust = 0,
                size = HOUSE_PT$text / .pt, family = HOUSE_FONT)
  }
  p
}

# Standardised residuals against the fitted value (posterior mean) with the
# binned mean +- 2 sd in 12 equal-count bins (make_figures.R, resid_panel).
resid_panel <- function(e, f, title, tag, ylim = c(-6, 6)) {
  z <- e / sd(e)
  ord <- order(f); bins <- split(seq_along(f)[ord], cut(seq_along(f), 12, labels = FALSE))
  b <- data.frame(f = sapply(bins, function(i) mean(f[i])), m = sapply(bins, function(i) mean(z[i])), s = sapply(bins, function(i) sd(z[i])))
  ggplot() + geom_point(data = data.frame(f, z), aes(f, z), colour = HOUSE_PRIMARY, alpha = 0.25, size = 0.9) +
    geom_ribbon(data = b, aes(f, ymin = m - 2 * s, ymax = m + 2 * s), fill = HOUSE_HIGHLIGHT, alpha = 0.2) +
    geom_line(data = b, aes(f, m), colour = HOUSE_HIGHLIGHT, linewidth = HOUSE_LINE) +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 0.44) +
    scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], length.out = 5), oob = oob_keep) +
    labs(x = "fitted value (posterior mean)", y = "standardised residual", title = title, tag = tag) + house_theme()
}

# copy a figure file to the manuscript folders (imgs or si-plots, and the reply)
copy_to_manuscript <- function(path, subdir = "imgs") {
  mr <- Sys.getenv("ARGOS_MANUSCRIPT_DIR", normalizePath(file.path(dirname(path), "..", "..", "..", "..", "..")))
  for (d in c(file.path(mr, subdir), file.path(mr, "reply-to-reviewers", "Figures"))) {
    if (dir.exists(d)) { file.copy(path, file.path(d, basename(path)), overwrite = TRUE); cat("copied to", file.path(d, basename(path)), "\n") }
  }
}
