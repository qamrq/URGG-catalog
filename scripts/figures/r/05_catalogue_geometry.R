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
#  SCRIPT  : 05_catalogue_geometry.R  -  Data loading and arc geometry for radial figures
#  AUTHOR  : QRQ
#  INPUT   : mag_provenance.tsv (PROV_F)
#  OUTPUT  : None (sourced by figure scripts)
# ==========================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })

# ------------------------------------------------------- data loading ------

#' load_catalogue_data(prov_path) — species- and family-level summary
#' tables for the final representative genomes, computed from
#' mag_provenance.tsv. Fails loudly (not silently) if any species in the
#' data has no family assignment.
load_catalogue_data <- function(prov_path, fam_lookup = family_of()) {
  if (!file.exists(prov_path)) {
    stop("mag_provenance.tsv not found at ", prov_path,
         "\nProduced by 33a_mag_provenance.R")
  }

  prov <- read_tsv(prov_path, show_col_types = FALSE) %>%
    filter(!is.na(species)) %>%
    mutate(species_disp = gsub("_", " ", species))

  species_df <- prov %>%
    filter(genome == final_representative) %>%
    distinct(genome, species_disp) %>%
    count(species_disp, name = "n") %>%
    mutate(family = fam_lookup[species_disp])

  if (any(is.na(species_df$family))) {
    missing_spp <- species_df %>% filter(is.na(family)) %>% pull(species_disp)
    stop("family_of() (02_catalogue_palette.R) is missing: ",
         paste(missing_spp, collapse = ", "),
         "\nAdd them there before proceeding.")
  }

  total <- sum(species_df$n)
  species_df <- species_df %>% mutate(pct = 100 * n / total)

  family_df <- species_df %>%
    group_by(family) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(pct = 100 * n / total) %>%
    arrange(desc(n))

  # Order species: grouped by family (in family_df's order), abundance
  # descending within family. This ordering is what makes a species ring
  # nest correctly beneath its parent family wedge in compute_ring_geometry.
  species_df <- species_df %>%
    arrange(match(family, family_df$family), desc(n))

  list(species = species_df, family = family_df, total = total,
      n_species = nrow(species_df), n_family = nrow(family_df))
}

# ------------------------------------------------------ geometry engine ----

#' compute_ring_geometry(df, value_col, r0, r1, start_at_top = TRUE)
#'
#' Adds start/end/mid angle (radians) and inner/outer radius columns to a
#' data frame, ready for ggforce::geom_arc_bar(aes(x0=0, y0=0, r0=r0, r=r1,
#' start=start, end=end, fill=...)).
#'
#' Angle convention matches ggforce::geom_arc_bar: 0 radians = 12 o'clock,
#' increasing CLOCKWISE. Rows are stacked in the data frame's existing row
#' order — callers control ordering (and therefore nesting alignment
#' between rings) by sorting the input df BEFORE calling this function.
#'
#' @param df        data frame with one row per wedge, already in the
#'                  desired stacking order
#' @param value_col name of the numeric column to size wedges by (e.g. "n")
#' @param r0, r1    inner/outer radius of this ring. Use WIDE ranges
#'                  (e.g. r0=1.15, r1=2.3) per plotting_instructions.md
#'                  point 5 — thick rings, not thin coord_polar-style bars.
compute_ring_geometry <- function(df, value_col, r0, r1) {
  vals <- df[[value_col]]
  total <- sum(vals)
  frac <- vals / total
  end_cum <- cumsum(frac) * 2 * pi
  start_cum <- end_cum - frac * 2 * pi

  df$start <- start_cum
  df$end   <- end_cum
  df$mid   <- (start_cum + end_cum) / 2
  df$r0    <- r0
  df$r1    <- r1

  # Label angle in degrees for geom_text angle= (arcs read outward from
  # centre). Flip labels on the left half of the circle by 180 deg so text
  # stays upright / reads left-to-right rather than upside-down.
  mid_deg <- df$mid * 180 / pi
  # ggforce's 0 = top, clockwise; convert to standard math degrees for
  # readability flip logic (measuring from 3 o'clock, counter-clockwise).
  math_deg <- 90 - mid_deg
  df$label_angle <- ifelse(math_deg < -90 | math_deg > 90, math_deg + 180, math_deg)
  df$label_hjust <- ifelse(mid_deg > 180, 1, 0)   # right half vs left half
  df$label_r <- (r0 + r1) / 2

  # Cartesian position for labels placed OUTSIDE the ring (e.g. leader-line
  # style for minor-but-still-shown wedges), at radius r1 + offset.
  df$label_x_out <- (r1 + 0.15) * sin(df$mid)
  df$label_y_out <- (r1 + 0.15) * cos(df$mid)
  # Cartesian position for labels placed INSIDE the ring (on the arc itself).
  df$label_x_in <- df$label_r * sin(df$mid)
  df$label_y_in <- df$label_r * cos(df$mid)

  df
}

#' validate_geometry_totals(family_geo, species_geo) — sanity check that
#' every family's arc [start,end] exactly bounds its member species'
#' arcs, i.e. the two rings genuinely nest rather than merely happening to
#' share a colour scheme. Stops with a specific mismatched family if not.
validate_geometry_totals <- function(family_geo, species_geo) {
  for (i in seq_len(nrow(family_geo))) {
    fam <- family_geo$family[i]
    sp_rows <- species_geo %>% filter(family == fam)
    fam_start <- family_geo$start[i]; fam_end <- family_geo$end[i]
    sp_start <- min(sp_rows$start); sp_end <- max(sp_rows$end)
    if (abs(fam_start - sp_start) > 1e-6 || abs(fam_end - sp_end) > 1e-6) {
      stop(sprintf(
        "Geometry mismatch for family '%s': family arc [%.4f, %.4f] does not ",
        fam, fam_start, fam_end),
        sprintf("match species arc bounds [%.4f, %.4f]. ", sp_start, sp_end),
        "Check that species_df was sorted by family BEFORE compute_ring_geometry() was called.")
    }
  }
  message("Geometry validated: all family arcs exactly bound their species arcs.")
  invisible(TRUE)
}
