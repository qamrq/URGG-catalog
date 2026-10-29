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
#  SCRIPT  : 02_catalogue_palette.R  -  Species and family colour palettes
#  AUTHOR  : QRQ
#  INPUT   : None (sourced by figure scripts)
#  OUTPUT  : None
# ==========================================================

#' species_palette() — hand-curated hex per host species.
#' First 9 entries are exactly as supplied in plotting_instructions.md;
#' the remaining 8 (Cervidae/Camelidae/Giraffidae) were added by hand
#' following the same dark-to-light-within-family logic.
species_palette <- function() {
  c(
    # --- Bovidae (blues -> greens -> olive/tan), as supplied -------------
    "Ovis aries"            = "#355C9A",
    "Pantholops hodgsonii"  = "#4776B3",
    "Bos taurus"            = "#5C8FD1",
    "Bos grunniens"         = "#7DA8E6",
    "Capra hircus"          = "#56B4B0",
    "Bubalus bubalis"       = "#5BAE61",
    "Bos indicus"           = "#7FC97F",
    "Procapra picticaudata" = "#C4C43A",
    "Bos frontalis"         = "#D7D7A8",
    # --- Cervidae (teal-greens, distinct hue from Bovidae's green tail) --
    "Hydropotes inermis"    = "#2E7D5B",
    "Alces alces"           = "#4FA37F",
    "Cervus elaphus"        = "#7FC7A6",
    "Rangifer tarandus"     = "#A9D9C6",
    "Alces americanus"      = "#CFE9DE",
    # --- Camelidae (oranges) ----------------------------------------------
    "Camelus bactrianus"    = "#C97B2E",
    "Camelus dromedarius"   = "#E0A667",
    # --- Giraffidae (brown) -----------------------------------------------
    "Giraffa camelopardalis"= "#8B5A2B"
  )
}

#' family_palette() — one base colour per family, matched to that family's
#' darkest/most abundant species so inner-ring and outer-ring hues cohere.
family_palette <- function() {
  c(
    "Bovidae"    = "#355C9A",
    "Cervidae"   = "#2E7D5B",
    "Camelidae"  = "#C97B2E",
    "Giraffidae" = "#8B5A2B"
  )
}

#' family_of() — species -> family lookup. Kept here (not duplicated per
#' figure script) so every figure agrees on the same taxonomy.
family_of <- function() {
  c(
    "Bos taurus"             = "Bovidae",
    "Bos indicus"            = "Bovidae",
    "Bos grunniens"          = "Bovidae",
    "Bos frontalis"          = "Bovidae",
    "Bubalus bubalis"        = "Bovidae",
    "Capra hircus"           = "Bovidae",
    "Ovis aries"             = "Bovidae",
    "Pantholops hodgsonii"   = "Bovidae",
    "Procapra picticaudata"  = "Bovidae",
    "Cervus elaphus"         = "Cervidae",
    "Alces alces"            = "Cervidae",
    "Alces americanus"       = "Cervidae",
    "Rangifer tarandus"      = "Cervidae",
    "Hydropotes inermis"     = "Cervidae",
    "Camelus bactrianus"     = "Camelidae",
    "Camelus dromedarius"    = "Camelidae",
    "Giraffa camelopardalis" = "Giraffidae"
  )
}

#' generate_family_shades(base_hex, n) — UTILITY for extending the palette
#' to species not yet hand-curated above. Steps dark -> medium -> light ->
#' blue-grey (desaturated), never fading to white, per
#' plotting_instructions.md point 2. Uses HSV space via base grDevices only
#' (no extra dependency).
#'
#' Usage: generate_family_shades("#355C9A", 5)
#' Returns 5 hex colours: 3 lightness steps of the base hue, then
#' progressively desaturated ("blue-grey") tones for any remainder.
generate_family_shades <- function(base_hex, n) {
  stopifnot(n >= 1)
  hsv_base <- rgb2hsv(col2rgb(base_hex))[, 1]
  h <- hsv_base["h"]; s <- hsv_base["s"]

  # First up to 3 steps: dark -> medium -> light, same hue/saturation.
  v_steps <- c(0.55, 0.75, 0.92)
  n_light <- min(n, 3)
  out <- hsv(h, s, v_steps[seq_len(n_light)])

  # Any additional steps: desaturate progressively at moderate lightness
  # ("blue-grey" family) rather than continuing to lighten toward white.
  n_remaining <- n - n_light
  if (n_remaining > 0) {
    s_steps <- s * seq(0.55, 0.15, length.out = n_remaining)
    out <- c(out, hsv(h, pmax(s_steps, 0.08), 0.72))
  }
  unname(out)
}

#' validate_palette_coverage(species_vec) — call this from a figure script
#' after computing the actual species present in the data, to fail loudly
#' if a species has no hand-curated colour rather than silently falling
#' back to a default grey.
validate_palette_coverage <- function(species_vec) {
  pal <- species_palette()
  missing <- setdiff(species_vec, names(pal))
  if (length(missing) > 0) {
    stop("species_palette() has no colour for: ", paste(missing, collapse = ", "),
         "\nAdd them to 02_catalogue_palette.R (hand-pick a hex, or use ",
         "generate_family_shades() to propose one), then re-run.")
  }
  invisible(TRUE)
}
