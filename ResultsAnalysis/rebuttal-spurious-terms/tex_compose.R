# Referee 1 Point 3 (R1.3): composite assembly by a generated standalone TeX
# document compiled with pdflatex (jupyter_tikz is not available here; the
# generated .tex is kept beside the figure so the layout is reproducible).
# Coordinates are in mm with y pointing down; final width 183 mm. Fonts:
# helvet only (tags 8 pt bold, labels 7 pt). The method colours are defined
# for the legend swatches; a tint is the colour mixed with white, which is
# what an alpha fill on white shows.

TEX_PREAMBLE <- c(
    "\\documentclass[border=0pt]{standalone}",
    "\\usepackage{graphicx}",
    "\\usepackage{tikz}",
    "\\usepackage{helvet}",
    "\\renewcommand{\\familydefault}{\\sfdefault}",
    "\\usetikzlibrary{calc}",
    "\\definecolor{ba}{HTML}{505ED5}",
    "\\definecolor{argos}{HTML}{D74D3D}",
    "\\definecolor{sindy}{HTML}{8EB63D}",
    "\\newcommand{\\Tag}[1]{{\\fontsize{8}{9}\\selectfont\\bfseries #1}}",
    "\\newcommand{\\Lab}[1]{{\\fontsize{7}{8}\\selectfont #1}}",
    "\\newcommand{\\Title}[1]{{\\fontsize{9}{10}\\selectfont #1}}",   # figure-level title, above the panel titles (2026-11-01)
    "\\newcommand{\\Leg}[1]{{\\fontsize{9}{10}\\selectfont #1}}",     # legend line, as the shared legend rows of the R figures
    "\\begin{document}",
    "\\begin{tikzpicture}[x=1mm,y=-1mm]"
)
tex_image <- function(x, y, width_mm, file, anchor = "north west") {
    sprintf("\\node[anchor=%s, inner sep=0pt] at (%.3f,%.3f) {\\includegraphics[width=%.3fmm]{%s}};", anchor, x, y, width_mm, file)
}
tex_text <- function(x, y, text, anchor = "north west", style = "Lab") {
    sprintf("\\node[anchor=%s, inner sep=0pt] at (%.3f,%.3f) {\\%s{%s}};", anchor, x, y, style, text)
}
tex_swatch <- function(x, y, colour, size = 2) {
    sprintf("\\fill[%s] (%.3f,%.3f) rectangle ++(%.3f,%.3f);", colour, x, y, size, size)
}
# One legend line, centred on cx with its top edge at y: three dark squares
# (true terms), three light squares (spurious terms), then the reading of the
# bar height; 9 pt text and 3 mm squares (Yuzheng, 2026-11-01: "legend larger").
# The squares are \\rule boxes inside the one text node, so the line is one
# node that TikZ centres (no nested pictures).
tex_fill_legend <- function(cx, y, text_after = "bar height: selection frequency over 100 trials", swatch_mm = 3) {
    sq <- function(col) sprintf("{\\color{%s}\\rule[-0.7mm]{%gmm}{%gmm}}", col, swatch_mm, swatch_mm)
    cols <- c("ba", "argos", "sindy")
    dark <- paste(sapply(cols, sq), collapse = "\\,")
    light <- paste(sapply(paste0(cols, "!40"), sq), collapse = "\\,")
    text <- sprintf("%s\\ true terms\\qquad %s\\ spurious terms\\qquad %s", dark, light, text_after)
    tex_text(cx, y, text, anchor = "north", style = "Leg")
}
# Writes tex_file (in the directory that holds the panel PDFs, referenced by
# relative paths) and compiles it there; returns the PDF path.
compose_tex <- function(tex_file, body) {
    writeLines(c(TEX_PREAMBLE, body, "\\end{tikzpicture}", "\\end{document}"), tex_file)
    dir <- dirname(tex_file)
    log <- system2("sh", c("-c", shQuote(sprintf("cd %s && pdflatex -interaction=nonstopmode -halt-on-error %s", shQuote(dir), shQuote(basename(tex_file))))), stdout = TRUE, stderr = TRUE)
    pdf_file <- sub("\\.tex$", ".pdf", tex_file)
    if (!file.exists(pdf_file) || any(grepl("^!", log))) stop(paste(c("pdflatex failed:", tail(log, 30)), collapse = "\n"))
    for (ext in c(".aux", ".log")) unlink(sub("\\.tex$", ext, tex_file))
    pdf_file
}
pdf_size_mm <- function(file) { # width, height of page 1 in mm
    info <- system2("pdfinfo", shQuote(file), stdout = TRUE)
    s <- grep("^Page size:", info, value = TRUE)
    v <- as.numeric(regmatches(s, gregexpr("[0-9.]+", s))[[1]][1:2])
    v * 25.4 / 72
}
