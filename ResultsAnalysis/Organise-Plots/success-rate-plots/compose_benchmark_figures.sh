#!/bin/bash
# Compose the benchmark figures (revised Figs. 2 and 3), the ablation figure and
# the second-order figure of the Supplementary Information from the success-rate
# rows drawn by the notebooks in ResultsAnalysis/success-rate-plots/ (house
# style, imgs/house/*.pdf), the attractor panels in simulation-plots/ and the
# shared legends. Specification: reply-to-reviewers/Notes/plotting-standard.md
# (blocks M2, M3, S-A, S-B; rules 4.1, 4.2.4, 4.7, 4.8).
#
# Geometry. The page is 183 mm wide plus a 1 mm standalone border, so that
# \includegraphics[width=\textwidth] in the manuscript shows the figure at
# 183 mm. Every R-drawn PDF (rows, legends) was authored on a canvas 2.5 times
# its final size and is included here at scale 0.4, so its text is 7 pt final
# (tags 8 pt bold, tick labels 6.5 pt). Benchmark rows are 138 mm wide
# (canvas 13.58 in) beside a 44 mm attractor column and a 1 mm gap; the SI rows
# are 183 mm wide (canvas 18 in). TeX text: tags \fontsize{8}{9}\bfseries,
# system labels \fontsize{8}{9}\bfseries, equation boxes \fontsize{7}{8}; sfmath
# sets the equations in the sans-serif face as well.
#
# Inputs (relative to this directory):
#   ../../success-rate-plots/imgs/house/<system>_success_rate_row.pdf
#   ../../success-rate-plots/imgs/house/<system>_ablation_success_rate_row[_notitle].pdf
#   ../../success-rate-plots/imgs/house/{duffing,forced-duffing,vdp}_success_rate_row[_notitle].pdf
#   ../../success-rate-plots/imgs/house/legend_methods.pdf, legend_ablation.pdf
#   simulation-plots/<System>_trajectory.pdf, simulation-plots/trajectory_colorbar.pdf
# Outputs: house/<name>.tex (kept beside this script) and house/<name>.pdf for
#   benchmark_main_1 (Lorenz, Thomas, Rossler, Dadras; tags a to d),
#   benchmark_main_2 (Aizawa, Sprott, Halvorsen; tags a to c),
#   ablation_success_rate (Lorenz, Rossler, Thomas; tags a to f),
#   second_order_success_rate (Duffing, forced Duffing, van der Pol; tags a to f).
# With MR_ROOT set (or auto-detected as the manuscript repository five levels
# up), the PDFs are also copied to $MR_ROOT/imgs (benchmark), $MR_ROOT/si-plots
# (SI) and $MR_ROOT/reply-to-reviewers/Figures (all four).
#
# Usage: ./compose_benchmark_figures.sh            (from any directory)
# Runtime: about 10 s (four pdflatex runs).
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
OUT="$DIR/house"; mkdir -p "$OUT"
ROWS="../../success-rate-plots/imgs/house"
MR_ROOT="${MR_ROOT:-$(cd "$DIR/../../../../.." && pwd)}"

# ---- final-size geometry in mm ------------------------------------------------
FIG_W=183          # figure width
ATT_W=44           # attractor column of the benchmark figures
COL_GAP=1          # gap between the attractor column and the row
ROW_W=138          # benchmark row width (= FIG_W - ATT_W - COL_GAP), rows authored at 138 mm x 2.5
ROW_ASPECT=0.30    # row height / row width, as in the notebooks
ROW_H=$(python3 -c "print(round($ROW_W * $ROW_ASPECT, 2))")
ROW_GAP=2          # vertical gap between rows
EQ_GAP=1           # gap between the equation box and the row
LEG_GAP=3          # gap between the last row and the legend line
CBAR_W=57          # width of the time-step colour bar (its 14 pt tick labels become 5.1 pt)
CBAR_TRIM_TOP=42   # pt trimmed off the top of trajectory_colorbar.pdf (its 30 pt title, ink from 6.5 to 35.5 pt); a 7 pt label replaces it
ATT_H=44           # attractor height (centred on the row; the row's own 1.5 mm margins absorb the overlap)
ATT_C=$(python3 -c "print($ATT_W / 2)")
ROW_C=$(python3 -c "print($ATT_W + $COL_GAP + $ROW_W / 2)")
FIG_C=$(python3 -c "print($FIG_W / 2)")
# Attractor PDFs are cropped to their ink (trim = left bottom right top, in pt,
# measured on a 72 dpi raster with 2 pt padding). The square pages (437.76 pt:
# Lorenz, Thomas, Rossler, Sprott) hold the 3-D plot only; the 482.5 x 582.5 pt
# pages (Aizawa, Dadras, Halvorsen) also carry their own "Time steps" colour bar
# below the plot, which the bottom trim removes (the shared bar sits at the foot).
att_trim() {
  case $1 in
    Aizawa|Dadras|Halvorsen) echo "66 123 1 82" ;;
    *) echo "60 16 2 74" ;;
  esac
}

preamble() {
cat <<TEX
\\documentclass[border={1mm 1mm 1mm 1mm}]{standalone}
\\usepackage{graphicx}
\\usepackage{tikz}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{helvet}
\\renewcommand{\\familydefault}{\\sfdefault}
\\usepackage[helvet]{sfmath}
\\usetikzlibrary{calc,positioning}
\\newcommand\\TagFont{\\fontsize{8}{9}\\selectfont\\bfseries}
\\newcommand\\LabelFont{\\fontsize{8}{9}\\selectfont\\bfseries}
\\newcommand\\TextFont{\\fontsize{7}{8}\\selectfont}
\\tikzset{
  img/.style={anchor=north west, inner sep=0pt, outer sep=0pt},
  tag/.style={anchor=west, inner sep=0pt, font=\\TagFont},
  syslabel/.style={inner sep=0pt, font=\\LabelFont},
  eqnbox/.style={anchor=south, inner sep=2.5pt, fill=gray!5, draw=gray!55, rounded corners=1.2pt, font=\\TextFont},
  txt/.style={inner sep=0pt, font=\\TextFont}
}
\\begin{document}
\\begin{tikzpicture}[x=1mm, y=1mm]
TEX
}

# benchmark row: $1 tag, $2 System (attractor file stem), $3 row file stem, $4 label, $5 equations, $6 previous tag or "top"
bench_row() {
  local tag=$1 sys=$2 stem=$3 label=$4 eqs=$5 prev=$6
  if [ "$prev" = top ]; then
    echo "  \\coordinate (eqtop$tag) at ($ROW_C, 0);"
  else
    echo "  \\coordinate (eqtop$tag) at (\$(row$prev.south -| $ROW_C,0) + (0,-$ROW_GAP)\$);"
  fi
cat <<TEX
  \\node[eqnbox, anchor=north] (eq$tag) at (eqtop$tag) {$eqs};
  \\node[img, anchor=north] (row$tag) at (\$(eq$tag.south) + (0,-$EQ_GAP)\$) {\\includegraphics[scale=0.4]{$ROWS/$stem.pdf}};
  \\node[img, anchor=center] at (\$(row$tag.center -| $ATT_C,0)\$) {\\includegraphics[height=${ATT_H}mm, width=${ATT_W}mm, keepaspectratio, trim=$(att_trim "$sys"), clip]{simulation-plots/${sys}_trajectory.pdf}};
  \\node[tag] at (0,0 |- eq$tag.center) {$tag};
  \\node[syslabel] at ($ATT_C,0 |- eq$tag.center) {$label};
TEX
}

# legend line of the benchmark figures: colour bar at the left, method legend centred under the row column
bench_foot() {
  local last=$1
cat <<TEX
  \\coordinate (foot) at (\$(row$last.south -| 0,0) + (0,-$LEG_GAP)\$);
  \\node[txt, anchor=north west] (cbarlab) at (foot) {Time steps};
  \\node[img] (cbar) at (\$(cbarlab.north east) + (1.5,0)\$) {\\includegraphics[width=${CBAR_W}mm, trim=0 0 0 ${CBAR_TRIM_TOP}pt, clip]{simulation-plots/trajectory_colorbar.pdf}};
  \\node[img, anchor=north] at (\$(foot -| $ROW_C,0) + (14,0)\$) {\\includegraphics[scale=0.4]{$ROWS/legend_methods.pdf}};
TEX
}

# SI row: $1 tag of the left panel, $2 tag of the right panel, $3 row file stem, $4 label, $5 previous left tag or "top"
si_row() {
  local ta=$1 tb=$2 stem=$3 label=$4 prev=$5
  if [ "$prev" = top ]; then
    echo "  \\coordinate (hd$ta) at (0, 0);"
  else
    echo "  \\coordinate (hd$ta) at (\$(row$prev.south -| 0,0) + (0,-$ROW_GAP)\$);"
  fi
cat <<TEX
  \\node[tag, anchor=north west] (tag$ta) at (hd$ta) {$ta};
  \\node[syslabel, anchor=north west] at (\$(tag$ta.north east) + (2,0)\$) {$label};
  \\node[tag, anchor=north west] at (\$(hd$ta) + ($FIG_C,0)\$) {$tb};
  \\node[img] (row$ta) at (\$(tag$ta.south west) + (0,-0.5)\$) {\\includegraphics[scale=0.4]{$ROWS/$stem.pdf}};
TEX
}

si_foot() {
  local last=$1 legend=$2
cat <<TEX
  \\node[img, anchor=north] at (\$(row$last.south -| $FIG_C,0) + (0,-$LEG_GAP)\$) {\\includegraphics[scale=0.4]{$ROWS/$legend.pdf}};
TEX
}

closing() { printf '\\end{tikzpicture}\n\\end{document}\n'; }

LorenzEqns='$\dot{x}_1=10(x_2 - x_1),\; \dot{x}_2= x_1(28 - x_3) - x_2,\; \dot{x}_3= x_1 x_2 - \tfrac{8}{3} x_3$'
ThomasEqns='$\dot{x}_1= \sin(x_2) - 0.208186 x_1,\; \dot{x}_2= \sin(x_3) - 0.208186 x_2,\; \dot{x}_3= \sin(x_1) - 0.208186 x_3$'
RosslerEqns='$\dot{x}_1=- x_2 - x_3,\; \dot{x}_2= x_1 + 0.2 x_2,\; \dot{x}_3=0.2 + x_3(x_1 - 5.7)$'
DadrasEqns='$\dot{x}_1= x_2 - 3 x_1 + 2.7 x_2 x_3,\; \dot{x}_2=1.7 x_2 - x_1 x_3 + x_3,\; \dot{x}_3=2 x_1 x_2 - 9 x_3$'
AizawaEqns='$\begin{gathered}\dot{x}_1=-3.5 x_2 + x_1 (-0.7 + x_3),\; \dot{x}_2=3.5 x_1 + x_2 (-0.7 + x_3) \\ \dot{x}_3=0.95x_3 + 0.65 + 0.1x_1^{3}x_3 - \tfrac{1}{3}x_3^{3} - (x_1^{2}+x_2^{2})(0.25x_3+1)\end{gathered}$'
SprottEqns='$\dot{x}_1= x_2 + 2.07 x_1 x_2 + x_1 x_3,\; \dot{x}_2= 1 -1.79 x_1^{2} + x_2 x_3,\; \dot{x}_3= x_1 - x_1^{2} - x_2^{2}$'
HalvorsenEqns='$\dot{x}_1=-1.89 x_1 - 4 x_2 - 4 x_3 - x_2^{2},\; \dot{x}_2=-1.89 x_2 - 4 x_3 - 4 x_1 - x_3^{2},\; \dot{x}_3=-1.89 x_3 - 4 x_1 - 4 x_2 - x_1^{2}$'

# ---- benchmark_main_1: Lorenz, Thomas, Rossler, Dadras -------------------------
{
  preamble
  bench_row a Lorenz  lorenz_success_rate_row  "Lorenz"      "$LorenzEqns"  top
  bench_row b Thomas  thomas_success_rate_row  "Thomas"      "$ThomasEqns"  a
  bench_row c Rossler rossler_success_rate_row "R\\\"ossler" "$RosslerEqns" b
  bench_row d Dadras  dadras_success_rate_row  "Dadras"      "$DadrasEqns"  c
  bench_foot d
  closing
} > "$OUT/benchmark_main_1.tex"

# ---- benchmark_main_2: Aizawa, Sprott, Halvorsen -------------------------------
{
  preamble
  bench_row a Aizawa    aizawa_success_rate_row    "Aizawa"    "$AizawaEqns"    top
  bench_row b Sprott    sprott_success_rate_row    "Sprott"    "$SprottEqns"    a
  bench_row c Halvorsen halvorsen_success_rate_row "Halvorsen" "$HalvorsenEqns" b
  bench_foot c
  closing
} > "$OUT/benchmark_main_2.tex"

# ---- ablation_success_rate: Lorenz, Rossler, Thomas (six screening variants) ----
{
  preamble
  si_row a b lorenz_ablation_success_rate_row          "Lorenz"      top
  si_row c d rossler_ablation_success_rate_row_notitle "R\\\"ossler" a
  si_row e f thomas_ablation_success_rate_row_notitle  "Thomas"      c
  si_foot e legend_ablation
  closing
} > "$OUT/ablation_success_rate.tex"

# ---- second_order_success_rate: Duffing, forced Duffing, van der Pol -----------
{
  preamble
  si_row a b duffing_success_rate_row                "Duffing"        top
  si_row c d forced-duffing_success_rate_row_notitle "Forced Duffing" a
  si_row e f vdp_success_rate_row_notitle            "Van der Pol"    c
  si_foot e legend_methods
  closing
} > "$OUT/second_order_success_rate.tex"

for name in benchmark_main_1 benchmark_main_2 ablation_success_rate second_order_success_rate; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory "$OUT" "$OUT/$name.tex" > "$OUT/$name.console.txt" 2>&1 \
    || { echo "pdflatex failed for $name; see $OUT/$name.console.txt and $OUT/$name.log"; tail -20 "$OUT/$name.console.txt"; exit 1; }
  rm -f "${OUT:?}/${name:?}.aux" "${OUT:?}/${name:?}.console.txt" "${OUT:?}/${name:?}.log"
  echo "$name: $(pdfinfo "$OUT/$name.pdf" | grep 'Page size')"
done

if [ -d "$MR_ROOT/imgs" ] && [ -d "$MR_ROOT/si-plots" ]; then
  cp "$OUT/benchmark_main_1.pdf" "$OUT/benchmark_main_2.pdf" "$MR_ROOT/imgs/"
  cp "$OUT/ablation_success_rate.pdf" "$OUT/second_order_success_rate.pdf" "$MR_ROOT/si-plots/"
  mkdir -p "$MR_ROOT/reply-to-reviewers/Figures"
  cp "$OUT"/{benchmark_main_1,benchmark_main_2,ablation_success_rate,second_order_success_rate}.pdf "$MR_ROOT/reply-to-reviewers/Figures/"
  echo "copied to $MR_ROOT/{imgs,si-plots,reply-to-reviewers/Figures}"
else
  echo "MR_ROOT=$MR_ROOT has no imgs/ and si-plots/; PDFs left in $OUT"
fi
