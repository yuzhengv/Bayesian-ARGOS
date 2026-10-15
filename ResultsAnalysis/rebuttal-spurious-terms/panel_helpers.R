# Referee 1 Point 3 (R1.3): the sweep panels shared by make_fig6.R (block a,
# means over the seven systems with a min-to-max ribbon) and
# make_si_spurious_by_system.R (one system per panel). House style throughout.
#
# The spurious-term count is drawn on a pseudo-logarithmic axis
# (asinh(y)/ln 10, linear near 0, logarithmic above about 1), because SINDy
# reaches 20 to 50 spurious terms per equation at low SNR or small n while
# ARGOS and Bayesian-ARGOS stay between 0 and 2; on a linear axis the
# comparison of the latter two is lost, and on a logarithmic axis the many
# exact zeros cannot be shown. The axis starts at 0 and is labelled with the
# counts themselves. The other two quantities use linear axes from 0.

source(file.path(root_path, "ResultsAnalysis", "plot-style", "house_style.R"))

SPUR_BREAKS <- c(0, 1, 2, 5, 10, 20, 50)
SPUR_TRANS <- scales::transform_pseudo_log(sigma = 0.5, base = 10)
SNR_XMIN <- 6   # the house SNR axis starts at 6 dB (standard 4.4)

scale_y_quantity <- function(quantity, ymax, name = NULL) {
    if (quantity == "spurious_per_eq") {
        scale_y_continuous(name = name, transform = SPUR_TRANS, breaks = SPUR_BREAKS, limits = c(0, ymax),
                           expand = expansion(mult = c(0, 0.04)), oob = oob_keep)
    } else {
        scale_y_continuous(name = name, limits = c(0, ymax), expand = expansion(mult = c(0, 0.04)), oob = oob_keep)
    }
}
# half-height of the SNR-axis separators: 2.5 % of the axis in transformed units
separator_extra_y <- function(quantity, ymax) {
    if (quantity == "spurious_per_eq") SPUR_TRANS$inverse(0.025 * SPUR_TRANS$transform(ymax)) else 0.025 * ymax
}
QUANTITY_TITLE <- c(spurious_per_eq = "Spurious terms per equation",
                    missed_per_eq = "True terms missed per equation",
                    mean_freq_spurious = "Mean selection frequency of spurious terms")
QUANTITY_YLAB <- c(spurious_per_eq = "Terms per equation", missed_per_eq = "Terms per equation",
                   mean_freq_spurious = "Selection frequency")

# d: columns x, method, y (and optionally ymin, ymax for a ribbon); sweep "n" or "snr".
sweep_panel <- function(d, sweep, quantity, ymax, title = NULL, ylab = NULL, xlab = TRUE, tag = NULL,
                        marker = HOUSE_MARKER, ribbon = "ymin" %in% names(d)) {
    d$method <- factor(d$method, levels = METHODS)
    if (sweep == "snr") { d <- d[d$x >= SNR_XMIN, ]; d$xx <- snr_x(d$x) } else d$xx <- d$x
    d <- d[order(d$method, d$xx), ]
    p <- ggplot(d, aes(x = xx, y = y, colour = method, fill = method, shape = method, group = method))
    if (ribbon) p <- p + geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, colour = NA, show.legend = FALSE)
    p <- p + geom_line(linewidth = HOUSE_LINE) + geom_point(size = marker) +
        method_scales() + house_theme() +
        scale_y_quantity(quantity, ymax, name = ylab) +
        labs(title = title, tag = tag) +
        theme(plot.title = element_text(size = HOUSE_PT$text, hjust = 0.5))
    if (sweep == "n") {
        p <- p + scale_x_log_n(name = if (xlab) expression(italic("n")) else NULL, limits = c(2, 5), expand = expansion(mult = 0.02))
    } else {
        p <- snr_axis(p, y_break = 0, extra_y = separator_extra_y(quantity, ymax), name = if (xlab) "SNR (dB)" else NULL)
    }
    if (!is.null(ylab)) p <- p + theme(axis.title.y = element_text(size = HOUSE_PT$text, angle = 90, margin = margin(r = 6)))
    p
}
row_label <- function(text, rot = 90) {
    patchwork::wrap_elements(full = grid::textGrob(text, rot = rot, gp = grid::gpar(fontsize = HOUSE_PT$text, fontfamily = HOUSE_FONT)))
}
