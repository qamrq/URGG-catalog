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
#  SCRIPT  : catalogue_export.R  -  Plot export used by 8a
#  AUTHOR  : QRQ
#  INPUT   : none
#  OUTPUT  : export_catalogue_plot()
# ==========================================================

suppressPackageStartupMessages(library(ggplot2))

CATALOGUE_PNG_DPI <- 300

extract_legend_catalogue <- function(plot) {
  if (requireNamespace("cowplot", quietly = TRUE) &&
      "get_plot_component" %in% getNamespaceExports("cowplot")) {
    legend <- tryCatch(
      cowplot::get_plot_component(plot, "guide-box", return_all = TRUE),
      error = function(e) NULL
    )
    if (!is.null(legend)) {
      if (is.list(legend) && !inherits(legend, "grob")) {
        legend <- Filter(function(g) !inherits(g, "zeroGrob"), legend)
        legend <- if (length(legend) > 0) legend[[1]] else NULL
      }
      if (!is.null(legend) && !inherits(legend, "zeroGrob")) return(legend)
    }
  }
  if (requireNamespace("cowplot", quietly = TRUE)) {
    legend <- tryCatch(cowplot::get_legend(plot), error = function(e) NULL)
    if (!is.null(legend) && !inherits(legend, "zeroGrob")) return(legend)
  }
  legend <- tryCatch({
    built <- ggplotGrob(plot)
    idx <- which(sapply(built$grobs, function(g) grepl("guide-box", g$name %||% "")))
    if (length(idx) > 0) built$grobs[[idx[1]]] else NULL
  }, error = function(e) NULL)
  if (!is.null(legend) && !inherits(legend, "zeroGrob")) return(legend)

  NULL
}
`%||%` <- function(a, b) if (is.null(a)) b else a

export_catalogue_plot <- function(plot, out_dir, basename,
                                  width = 20, height = 14,
                                  png_dpi = CATALOGUE_PNG_DPI,
                                  save_legend = TRUE,
                                  legend_width = NULL, legend_height = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  f <- file.path(out_dir, paste0(basename, ".png"))
  ggsave(f, plot, width = width, height = height, dpi = png_dpi, bg = "white")

  cat("Exported ", basename, ":\n  ", f, "\n", sep = "")

  if (isTRUE(save_legend)) {
    legend_grob <- extract_legend_catalogue(plot)
    if (!is.null(legend_grob)) {
      PADDING_FACTOR <- 1.25
      MIN_SIZE_IN <- 2

      measured_w <- tryCatch(
        grid::convertWidth(grid::grobWidth(legend_grob), "in", valueOnly = TRUE),
        error = function(e) NA_real_)
      measured_h <- tryCatch(
        grid::convertHeight(grid::grobHeight(legend_grob), "in", valueOnly = TRUE),
        error = function(e) NA_real_)

      use_w <- if (!is.null(legend_width)) legend_width else
        if (!is.na(measured_w)) max(measured_w * PADDING_FACTOR, MIN_SIZE_IN) else 10
      use_h <- if (!is.null(legend_height)) legend_height else
        if (!is.na(measured_h)) max(measured_h * PADDING_FACTOR, MIN_SIZE_IN) else 4

      if (is.null(legend_width) && is.na(measured_w)) {
        message("  Could not measure legend width automatically -- using a wide fallback (10in). ",
                "If it's still clipped or has excess empty space, pass legend_width= explicitly.")
      }

      f_legend <- file.path(out_dir, paste0(basename, "_legend.png"))
      if (requireNamespace("cowplot", quietly = TRUE)) {
        ggsave(f_legend, cowplot::ggdraw(legend_grob),
              width = use_w, height = use_h, dpi = png_dpi, bg = "white")
      } else {
        png(f_legend, width = use_w * png_dpi, height = use_h * png_dpi, res = png_dpi)
        grid::grid.newpage(); grid::grid.draw(legend_grob)
        dev.off()
      }
      cat("  ", f_legend, sprintf(" (legend only, auto-sized to %.1fin x %.1fin)\n", use_w, use_h), sep = "")
    } else {
      message("  (no legend found on this plot -- skipping legend export)")
    }
  }

  invisible(f)
}