#!/usr/bin/env Rscript
#                 ############################
#               ##                            ##
#             ##  ####                    ####  ##
#           ##   ######                  ######   ##
#         ##    ########                ########    ##
#       ##      ##########            ##########      ##
#     ##        ####  ####################  ####        ##
#   ##          #####  ##################  #####          ##
# ##            ################################            ##
#   ##            ############################            ##
#     ##              ####################              ##
#       ##             ##  ##########  ##             ##
#         ##           ##   ########   ##           ##
#           ##          #    ######    #          ##
#             ##              ####              ##
#               ##                            ##
#                 ############################
#
# ==========================================================
#  PROJECT : URGG catalog
#  SCRIPT  : 04_catalogue_labels.R  -  Label placement rules for radial figures
#  AUTHOR  : QRQ
#  INPUT   : None (sourced by figure scripts)
#  OUTPUT  : None
# ==========================================================

suppressPackageStartupMessages(library(ggplot2))

# Wedges below this share of the WHOLE catalogue are not labelled on the
# plot itself — they still appear in the legend panel. 3-4% per
# plotting_instructions.md point 5/3; 3.5 is the working default.
LABEL_MIN_PCT <- 3.5

#' should_label(pct_vec) — TRUE/FALSE per wedge, using LABEL_MIN_PCT.
should_label <- function(pct_vec) pct_vec >= LABEL_MIN_PCT

#' italicise(x) — wraps a species name (or vector of names) for use inside
#' plotmath-aware ggplot text (e.g. geom_text with parse = TRUE), and
#' separately provides a plain-text version formatter for non-plotmath
#' contexts (e.g. plotly, which uses HTML tags instead).
italicise_plotmath <- function(x) {
  # Only italicise genuine two-word binomials; leave family/other labels
  # (e.g. "Bovidae", "Other 5 hosts") in upright text.
  is_binomial <- grepl("^[A-Z][a-z]+ [a-z]+$", x)
  out <- x
  out[is_binomial] <- paste0("italic('", x[is_binomial], "')")
  out
}

italicise_html <- function(x) {
  is_binomial <- grepl("^[A-Z][a-z]+ [a-z]+$", x)
  out <- x
  out[is_binomial] <- paste0("<i>", x[is_binomial], "</i>")
  out
}

#' build_center_annotation(total, subtitle_lines, size_total, size_sub)
#'
#' Returns a list of two ggplot annotate() layers: a large bold total and
#' smaller descriptive lines beneath it — for the hole of a donut/sunburst.
#' Add both elements to a plot with `plot + build_center_annotation(...)`.
#'
#' @param total          numeric or pre-formatted string, e.g. "25,451"
#' @param subtitle_lines character vector, one per line, e.g.
#'                       c("species representatives", "95% ANI catalogue", "17 host species")
#' @param x, y           position in the plot's coordinate space (defaults
#'                       assume a standard geom_arc_bar hole at the origin)
build_center_annotation <- function(total, subtitle_lines,
                                    x = 0, y = 0,
                                    size_total = 8, size_sub = 3.1,
                                    gap = 0.5) {
  total_txt <- if (is.numeric(total)) format(total, big.mark = ",") else total
  sub_txt <- paste(subtitle_lines, collapse = "\n")
  list(
    annotate("text", x = x, y = y + gap, label = total_txt,
            size = size_total, fontface = "bold", colour = "grey12"),
    annotate("text", x = x, y = y - gap, label = sub_txt,
            size = size_sub, colour = "grey35", lineheight = 1.05)
  )
}
