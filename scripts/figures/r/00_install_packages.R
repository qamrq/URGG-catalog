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
#  SCRIPT  : 00_install_packages.R  -  Package check and install for the catalogue figures
#  AUTHOR  : QRQ
#  INPUT   : None
#  OUTPUT  : None (loads packages, prints session info)
# ==========================================================

MIN_R_VERSION <- "4.1.0"
if (getRversion() < MIN_R_VERSION) {
  warning(sprintf("R %s detected; this framework was written against R >= %s. Proceed with caution.",
                  getRversion(), MIN_R_VERSION))
}

REQUIRED_PKGS <- c(
  "dplyr", "readr", "tidyr", "stringr",     # data handling
  "ggplot2", "ggforce",                     # geometry engine (geom_arc_bar)
  "cowplot",                                # panel assembly, legend extraction
  "scales",                                 # label formatting
  "treemapify",                             # treemap geometry (31i)
  "plotly", "htmlwidgets"                   # interactive export option
)

OPTIONAL_PKGS <- c("webshot2", "showtext", "sysfonts")  # static export / custom fonts

.catalogue_ensure_packages <- function(pkgs, install_if_missing = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    if (install_if_missing) {
      message("Installing missing packages: ", paste(missing, collapse = ", "))
      install.packages(missing, repos = "https://cloud.r-project.org")
    } else {
      message("Missing (optional) packages, not installed: ", paste(missing, collapse = ", "))
    }
  }
  invisible(NULL)
}

.catalogue_ensure_packages(REQUIRED_PKGS, install_if_missing = TRUE)
.catalogue_ensure_packages(OPTIONAL_PKGS, install_if_missing = FALSE)

suppressPackageStartupMessages({
  for (p in REQUIRED_PKGS) library(p, character.only = TRUE)
})

# Optional custom fonts (fonts/ subdirectory) via showtext, if available and
# the user has dropped .ttf files there. Silent no-op otherwise.
#
# CATALOGUE_R_DIR should be set by the calling figure script (31g/31h/31i)
# BEFORE sourcing this file, e.g.:
#   CATALOGUE_R_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_3_analysis/scripts/results_scripts/R"
#   source(file.path(CATALOGUE_R_DIR, "00_install_packages.R"))
# This avoids relying on sys.frame()$ofile, which does not reliably resolve
# a script's own path when run via `Rscript` rather than `source()`.
if (!exists("CATALOGUE_R_DIR")) {
  CATALOGUE_R_DIR <- getwd()
  message("CATALOGUE_R_DIR not set by caller — defaulting to getwd() (", CATALOGUE_R_DIR, "). ",
         "Set it explicitly before sourcing this file if fonts/ lives elsewhere.")
}
FONTS_DIR <- file.path(CATALOGUE_R_DIR, "fonts")

if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE) &&
    dir.exists(FONTS_DIR)) {
  ttfs <- list.files(FONTS_DIR, pattern = "\\.ttf$", full.names = TRUE, ignore.case = TRUE)
  if (length(ttfs) > 0) {
    for (f in ttfs) {
      fam_name <- tools::file_path_sans_ext(basename(f))
      sysfonts::font_add(fam_name, f)
    }
    showtext::showtext_auto()
    message("Loaded ", length(ttfs), " custom font(s) from ", FONTS_DIR)
  }
}

cat("\n=== Session info (for reproducibility) ===\n")
cat("R version:", as.character(getRversion()), "\n")
cat("Platform:", R.version$platform, "\n")
cat("Required packages loaded:\n")
for (p in REQUIRED_PKGS) {
  cat(sprintf("  %-14s %s\n", p, as.character(packageVersion(p))))
}
cat("===========================================\n\n")
