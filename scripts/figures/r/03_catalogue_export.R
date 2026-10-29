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
#  SCRIPT  : 03_catalogue_export.R  -  Export helper for catalogue figures (PDF, SVG, PNG)
#  AUTHOR  : QRQ
#  INPUT   : None (sourced by figure scripts)
#  OUTPUT  : Figure files in the output folder given to the helper
# ==========================================================

suppressPackageStartupMessages(library(ggplot2))

#' export_catalogue_plot(plot, out_dir, basename, width, height, formats)
#'
#' @param plot     a ggplot object (or grid/cowplot object accepted by ggsave)
#' @param out_dir  destination directory (created if absent)
#' @param basename filename without extension, e.g. "Main_Fig1_sunburst"
#' @param width, height  in inches
#' @param formats  any of "pdf", "svg", "png" (default: all three)
#' @param png_dpi  resolution for PNG export (default 600, per
#'                 plotting_instructions.md point 8)
export_catalogue_plot <- function(plot, out_dir, basename,
                                  width = 9, height = 9,
                                  formats = c("pdf", "svg", "png"),
                                  png_dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character(0)

  if ("pdf" %in% formats) {
    f <- file.path(out_dir, paste0(basename, ".pdf"))
    ggsave(f, plot, width = width, height = height, device = cairo_pdf)
    written <- c(written, f)
  }
  if ("svg" %in% formats) {
    if (requireNamespace("svglite", quietly = TRUE)) {
      f <- file.path(out_dir, paste0(basename, ".svg"))
      ggsave(f, plot, width = width, height = height, device = svglite::svglite)
      written <- c(written, f)
    } else {
      message("Skipping SVG export for ", basename, " — install.packages('svglite') to enable it.")
    }
  }
  if ("png" %in% formats) {
    f <- file.path(out_dir, paste0(basename, ".png"))
    ggsave(f, plot, width = width, height = height, dpi = png_dpi, bg = "white")
    written <- c(written, f)
  }

  cat("Exported ", basename, ":\n", sep = "")
  for (f in written) cat("  ", f, "\n", sep = "")
  invisible(written)
}
