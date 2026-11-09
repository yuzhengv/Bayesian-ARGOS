# House style for every data figure of the revised manuscript, its Supplementary
# Information and the reply (plotting standard: reply-to-reviewers/Notes/
# plotting-standard.md, sections 4 and 5). One theme, one palette, one axis
# convention; sourced by the success-rate notebooks and by every figure script
# under ResultsAnalysis/. Never modifies package code.
#
# Authoring route (standard 4.2.3, decision D1): figures are drawn on a canvas
# 2.5 times their final width (18 in for a 183 mm full-width figure, 8.75 in for
# an 89 mm one-column figure) with every size multiplied by 2.5, so that at
# final size tags are 8 pt bold, titles, axis titles and legend text 7 pt, tick
# labels 6.5 pt, ticks 1 mm and spines 0.5 pt.
#
# Usage:  source(file.path(<ResultsAnalysis>, "plot-style", "house_style.R"))
#   house_theme()            the ggplot2 theme (legend off by default)
#   HOUSE_COLOURS, HOUSE_SHAPES, HOUSE_EXTRA, HOUSE_GREY, SUCCESS_BAND
#   scale_x_log_n()          n axis, 10^2 ... 10^5, labels 10^k
#   snr_axis(p, y_break, ...)  SNR axis 10 ... 60 dB plus infinity beyond a break (spine drawn as two segments)
#   snr_x(s)                 position of an SNR value (Inf -> INF_X)
#   house_legend_theme(), house_legend_guides(nrow)   the shared legend row, 9 pt final (2026-11-01)
#   row_label(text)          a rotated row label as a narrow patchwork element
#   success_rate_layers()    the dashed 0.8 line and the band above it
#   scale_y_success()        y axis 0 ... 1
#   scale_y_log10_powers()   y axis with 10^k labels
#   method_scales(methods)   colour, fill and shape scales for the three methods
#   save_figure(p, path, width_mm, aspect)   ggsave on the k = 2.5 canvas
#   save_legend(legend, path)  a legend grob alone, at its natural size on the canvas

suppressPackageStartupMessages({
  library(ggplot2); library(scales); library(latex2exp); library(ggh4x)
})

HOUSE_K <- 2.5                      # canvas = HOUSE_K x final size
HOUSE_PT <- list(tag = 8 * HOUSE_K, text = 7 * HOUSE_K, tick = 6.5 * HOUSE_K)
HOUSE_FONT <- "Helvetica"           # pdf() device base font; Arial is the fallback

# --------------------------------------------------------------- colours, shapes
HOUSE_COLOURS <- c("Bayesian-ARGOS" = "#505ed5", "ARGOS" = "#d74d3d", "SINDy" = "#8eb63d")
HOUSE_SHAPES  <- c("Bayesian-ARGOS" = 16, "ARGOS" = 15, "SINDy" = 17)
HOUSE_EXTRA   <- c("#FA5524", "#EF60AD", "#009E9E", "#8B5A2B")   # fourth to seventh series
HOUSE_PRIMARY <- "#505ed5"; HOUSE_HIGHLIGHT <- "#d74d3d"
HOUSE_GREY    <- "grey75"
SUCCESS_BAND  <- "#BDB6FB"
HOUSE_SERIES  <- c(HOUSE_COLOURS[["Bayesian-ARGOS"]], HOUSE_COLOURS[["ARGOS"]], HOUSE_COLOURS[["SINDy"]], HOUSE_EXTRA)

# ----------------------------------------------------------------------- theme
house_theme <- function(legend = "none", base = HOUSE_PT$text) {
  theme(
    text = element_text(family = HOUSE_FONT, size = base, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.44),      # 0.5 pt at final size
    axis.ticks = element_line(colour = "black", linewidth = 0.44),
    axis.ticks.length = unit(0.25, "cm"),                               # 1 mm at final size
    axis.text = element_text(size = HOUSE_PT$tick, colour = "black"),
    axis.title.x = element_text(size = base, margin = margin(t = 6)),
    axis.title.y = element_text(size = base, angle = 90, vjust = 0.5, margin = margin(r = 6)),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), panel.background = element_blank(),
    plot.background = element_blank(),
    plot.title = element_text(size = base, hjust = 0.5, face = "plain"),
    plot.tag = element_text(size = HOUSE_PT$tag, face = "bold"),
    legend.key = element_blank(), legend.background = element_blank(),
    legend.text = element_text(size = base), legend.title = element_text(size = base, face = "bold"),
    legend.position = legend,
    plot.margin = margin(t = 0.45, r = 0.2, b = 0.2, l = 0.15, unit = "in")
  )
}

# ----------------------------------------------------------------------- axes
LOG_N_BREAKS <- seq(2, 5, 0.5)
LOG_N_LABELS <- c(expression(10^2), expression(10^2.5), expression(10^3), expression(10^3.5),
                  expression(10^4), expression(10^4.5), expression(10^5))
scale_x_log_n <- function(name = expression(italic("n")), limits = NULL, breaks = LOG_N_BREAKS, labels = LOG_N_LABELS, ...) {
  scale_x_continuous(name = name, breaks = breaks, labels = labels, limits = limits, ...)
}

# SNR axis (standard 4.4): finite ticks 10 ... 60 dB; the noiseless case at INF_X
# (the axis starts at 1 dB, the first SNR of the sweep: Yuzheng, 2026-11-09; it started at 6 dB before)
# beyond a break drawn as two slanted separators on an otherwise empty spine.
INF_X <- 73
SNR_BREAK <- c(65.5, 69.5)
SEP_ANGLE <- 75           # slant of the separators
SEP_HALF <- 0.05          # half-height of a separator as a fraction of the (transformed) y range
SEP_ADD_Y <- 0.1 * sin(SEP_ANGLE * pi / 180) / 2   # the half-height the old asymmetric strokes added (compatibility)
snr_x <- function(s) ifelse(is.infinite(s), INF_X, s)

# Two slanted strokes centred on the spine (Yuzheng, 2026-11-01: the spine
# meets each separator at mid-height): from (x - dx, y - half) to
# (x + dx + extra_x, y + half), y and half in the transformed y units.
create_separators <- function(x, extra_x, y, half, angle = SEP_ANGLE, length = .1) {
  add_x <- length * cos(angle * pi / 180)
  list(x = x - add_x, xend = x + add_x + extra_x, y = rep(y - half, length(x)), yend = rep(y + half, length(x)))
}

# p: a ggplot whose x is snr_x(SNR); y_break: the y value at which the spine is
# drawn (the lower limit of the y scale: 0 for a 0-1 axis, 1 for a log axis
# starting at 10^0, ...). The x spine of the theme is switched off and drawn as
# two segments, so that the spine is empty between the two separators whatever
# the ggplot2 version (ggh4x's guide_axis_truncated is deprecated).
# Separator size (2026-11-01): pass the y limits as y_range (and the scale's
# transformation as ytrans, or ylog = TRUE for log10) and each separator spans
# 2 * SEP_HALF = 10 % of the axis height, centred on the spine, whatever the
# quantity and its scale. Without y_range the former arguments apply: extra_y
# and sep_y_scale give the half-height (SEP_ADD_Y * sep_y_scale + extra_y) / 2,
# which is 0.049 for the 0-1 success-rate axis (extra_y = 0.05) and 0.049 R
# for a log axis of R decades (extra_y = 0.05 R, sep_y_scale = R).
snr_axis <- function(p, y_break = 0, extra_y = 0.05, ticks = seq(10, 60, 10), xmin = 1, name = "SNR (dB)",
                     ylog = FALSE, sep_y_scale = 1, ytrans = NULL, y_range = NULL) {
  if (ylog && is.null(ytrans)) ytrans <- scales::transform_log10()
  tf <- if (is.null(ytrans)) identity else ytrans$transform
  inv <- if (is.null(ytrans)) identity else ytrans$inverse
  half <- if (!is.null(y_range)) SEP_HALF * abs(diff(tf(y_range))) else (SEP_ADD_Y * sep_y_scale + extra_y) / 2
  seg <- create_separators(SNR_BREAK, extra_x = 1, y = tf(y_break), half = half)
  xr <- c(xmin, INF_X + 1)
  p + scale_x_continuous(name = name, limits = xr, breaks = c(ticks, INF_X),
                         labels = c(ticks, TeX("$\\infty$")), expand = expansion(mult = 0), oob = oob_keep) +
    theme(axis.line.x = element_blank()) +
    annotate("segment", x = c(xr[1], SNR_BREAK[2] + 0.5), xend = c(SNR_BREAK[1] + 0.5, xr[2]),
             y = y_break, yend = y_break, linewidth = 0.44, colour = "black") +
    annotate("segment", x = seg$x, xend = seg$xend, y = inv(seg$y), yend = inv(seg$yend), linewidth = 0.44, colour = "black") +
    coord_cartesian(clip = "off")
}

success_rate_layers <- function(threshold = 0.8) {
  list(geom_rect(data = data.frame(xmin = -Inf, xmax = Inf, ymin = threshold, ymax = Inf),
                 aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
                 alpha = 0.3, fill = SUCCESS_BAND),
       geom_hline(yintercept = threshold, linetype = 2, linewidth = 0.44))
}
scale_y_success <- function(name = "Success rate", ...) {
  scale_y_continuous(name = name, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02)), oob = oob_keep, ...)
}
scale_y_log10_powers <- function(name, limits = NULL, step = 1, lo = -8, hi = 16) {
  scale_y_continuous(name = name, trans = "log10", breaks = 10^seq(lo, hi, step),
                     labels = trans_format("log10", math_format(10^.x)), limits = limits,
                     expand = expansion(mult = c(0, 0.02)), oob = oob_keep)
}

# ------------------------------------------------------------------- methods
method_scales <- function(methods = names(HOUSE_COLOURS), labels = methods, name = NULL) {
  list(scale_colour_manual(values = HOUSE_COLOURS[methods], breaks = methods, labels = labels, name = name),
       scale_fill_manual(values = HOUSE_COLOURS[methods], breaks = methods, labels = labels, name = name),
       scale_shape_manual(values = HOUSE_SHAPES[methods], breaks = methods, labels = labels, name = name))
}
# point size for 100-trial success-rate sweeps (standard 4.6.1) and summary curves (4.6.2)
HOUSE_POINT <- 3; HOUSE_MARKER <- 4; HOUSE_LINE <- 1.4; HOUSE_CONTEXT_LINE <- 0.7

# ------------------------------------------------------------------- legends
# The shared legend row of a figure (standard 4.7.1) is set larger than the
# panel text (Yuzheng, 2026-11-01: "make this legend larger in all the plots
# using it"): text 9 pt final, keys enlarged to match. house_legend_theme()
# gives the theme elements (position, text, key size, spacing);
# house_legend_guides(nrow) the guides with enlarged point keys for the colour,
# fill and shape mappings of the method and variant legends.
HOUSE_PT$legend <- 9 * HOUSE_K
HOUSE_LEGEND_POINT <- 4.5
house_legend_theme <- function(position = "bottom", key_width = 0.6, key_height = 0.3) {
  theme(legend.position = position, legend.direction = "horizontal", legend.title = element_blank(),
        legend.text = element_text(size = HOUSE_PT$legend),
        legend.key.width = unit(key_width, "in"), legend.key.height = unit(key_height, "in"),
        legend.key.spacing.x = unit(0.25, "in"), legend.margin = margin(2, 2, 2, 2))
}
house_legend_guides <- function(nrow = 1, point_size = HOUSE_LEGEND_POINT, linewidth = NULL) {
  oa <- list(size = point_size); if (!is.null(linewidth)) oa$linewidth <- linewidth
  g <- guide_legend(nrow = nrow, override.aes = oa)
  guides(colour = g, fill = g, shape = g)
}

# A rotated row label (the system of a row, standard 4.8.2) as a narrow
# patchwork element placed at the left of the row with plot_layout(widths =
# c(0.05, 1, 1)), as the Supplementary spurious-terms figure draws it; 7 pt final.
row_label <- function(text, rot = 90, size = HOUSE_PT$text, face = "plain") {
  patchwork::wrap_elements(full = grid::textGrob(text, rot = rot, gp = grid::gpar(fontsize = size, fontfamily = HOUSE_FONT, fontface = face)))
}

# ----------------------------------------------------------------------- save
# width_mm: the final width (183 or 89); aspect = height / width of the figure.
# The canvas is HOUSE_K times the final size. Record width_mm and aspect in the
# header of the calling script (standard 4.12.5).
save_figure <- function(p, path, width_mm = 183, aspect = 0.5, dpi = 300, device = "pdf") {
  # device: "pdf" (default, Helvetica referenced but not embedded, the infinity
  # sign from the Symbol font) or cairo_pdf (embeds the Helvetica substitute
  # Nimbus Sans, infinity from the same face; used by the success-rate rows).
  w <- width_mm / 25.4 * HOUSE_K
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, p, width = w, height = w * aspect, units = "in", dpi = dpi, device = device, family = HOUSE_FONT)
  invisible(path)
}

# A legend grob (cowplot::get_legend(p)) saved on its own at its natural size on
# the k = 2.5 canvas, so that a composite includes it at scale 1 / HOUSE_K
# (pdflatex: \includegraphics[scale=0.4]) and the legend text is 7 pt final.
save_legend <- function(legend, path, device = cairo_pdf, pad_in = 0.05) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  w <- grid::convertWidth(sum(legend$widths), "in", valueOnly = TRUE)
  h <- grid::convertHeight(sum(legend$heights), "in", valueOnly = TRUE)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, cowplot::ggdraw(legend), width = w + 2 * pad_in, height = h + 2 * pad_in, units = "in",
         device = device, family = HOUSE_FONT)
  invisible(path)
}
