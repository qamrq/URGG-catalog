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

# Custom fonts in figures/r/fonts (.ttf) are loaded with showtext when the
# packages are available. Nothing happens if the folder is empty or missing.
FONTS_DIR <- file.path(FIG_R_DIR, "fonts")

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
