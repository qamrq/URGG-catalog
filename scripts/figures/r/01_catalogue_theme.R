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
#  SCRIPT  : 01_catalogue_theme.R  -  Shared ggplot theme for catalogue figures
#  AUTHOR  : QRQ
#  INPUT   : None (sourced by figure scripts)
#  OUTPUT  : None
# ==========================================================

suppressPackageStartupMessages(library(ggplot2))

CATALOGUE_FONT_FAMILY <- "sans"   # overridden automatically if a custom font
                                  # was loaded by 00_install_packages.R's
                                  # fonts/ handling — set explicitly here if
                                  # you want a specific loaded family, e.g.
                                  # CATALOGUE_FONT_FAMILY <- "Arial"

CATALOGUE_COLORS <- list(
  text_primary   = "grey12",
  text_secondary = "grey35",
  grid           = "grey88",
  panel_bg       = "white",
  wedge_border   = "white"
)

#' Shared theme for radial (sunburst/donut) figures — no axes, no grid,
#' generous margins for labels that sit outside the plot area.
theme_catalogue_radial <- function(base_size = 11) {
  theme_void(base_size = base_size, base_family = CATALOGUE_FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", size = base_size * 1.4,
                                colour = CATALOGUE_COLORS$text_primary,
                                hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size * 1.0,
                                   colour = CATALOGUE_COLORS$text_secondary,
                                   hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(14, 14, 14, 14),
      legend.position = "none"   # legends are a separate assembled panel
    )
}

#' Shared theme for standard cartesian figures (bar charts, scatter, etc.)
#' used alongside the radial figures in the same section.
theme_catalogue_cartesian <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = CATALOGUE_FONT_FAMILY) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = CATALOGUE_COLORS$grid, linewidth = .3),
      plot.title = element_text(face = "bold", size = base_size * 1.15,
                                colour = CATALOGUE_COLORS$text_primary),
      plot.subtitle = element_text(size = base_size * 0.95,
                                   colour = CATALOGUE_COLORS$text_secondary),
      axis.text = element_text(colour = CATALOGUE_COLORS$text_primary),
      axis.title = element_text(colour = CATALOGUE_COLORS$text_primary)
    )
}

#' Legend-panel theme: no axes/grid, italic species names, compact.
theme_catalogue_legend <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = CATALOGUE_FONT_FAMILY) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(face = "italic", colour = CATALOGUE_COLORS$text_primary),
      axis.title = element_blank(),
      plot.margin = margin(5, 5, 5, 0)
    )
}
