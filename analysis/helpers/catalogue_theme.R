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
#  SCRIPT  : catalogue_theme.R  -  Plot theme used by 8a
#  AUTHOR  : QRQ
#  INPUT   : none
#  OUTPUT  : theme_catalogue_*() functions
# ==========================================================

suppressPackageStartupMessages(library(ggplot2))

CATALOGUE_FONT_FAMILY <- local({
  requested <- "Times New Roman"
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    available <- systemfonts::system_fonts()$family
    if (requested %in% available) {
      requested
    } else {
      message("'", requested, "' not found among installed system fonts. ",
              "Falling back to 'serif' (will resolve to whatever serif font ",
              "is actually installed -- often Liberation Serif on Linux, a ",
              "close Times New Roman substitute). To get the exact font, ",
              "install it, e.g. on Debian/Ubuntu: ",
              "'sudo apt install fonts-liberation ttf-mscorefonts-installer' ",
              "then re-run this script.")
      "serif"
    }
  } else {
    message("Package 'systemfonts' not installed -- cannot verify whether ",
            "'", requested, "' is available. Using 'serif' as a safe default ",
            "(install.packages('systemfonts') to enable a proper check and ",
            "use the exact font name when available).")
    "serif"
  }
})

CATALOGUE_COLORS <- list(
  text_primary   = "grey12",
  text_secondary = "grey35",
  grid           = "grey88",
  panel_bg       = "white",
  wedge_border   = "white"
)

theme_catalogue_radial <- function(base_size = 16) {
  theme_void(base_size = base_size, base_family = CATALOGUE_FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", size = base_size * 1.4,
                                colour = CATALOGUE_COLORS$text_primary,
                                hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size * 1.0,
                                   colour = CATALOGUE_COLORS$text_secondary,
                                   hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(14, 14, 14, 14),
      legend.position = "none"
    )
}

theme_catalogue_cartesian <- function(base_size = 22) {
  theme_minimal(base_size = base_size, base_family = CATALOGUE_FONT_FAMILY) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = CATALOGUE_COLORS$grid, linewidth = .3),
      plot.title = element_text(face = "bold", size = base_size * 1.15,
                                colour = CATALOGUE_COLORS$text_primary),
      plot.subtitle = element_text(size = base_size * 0.95,
                                   colour = CATALOGUE_COLORS$text_secondary),
      axis.text = element_text(size = rel(1), colour = CATALOGUE_COLORS$text_primary),
      axis.title = element_text(size = rel(1), colour = CATALOGUE_COLORS$text_primary),
      legend.text = element_text(size = rel(0.95)),
      legend.title = element_text(size = rel(1))
    )
}

theme_catalogue_legend <- function(base_size = 13) {
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