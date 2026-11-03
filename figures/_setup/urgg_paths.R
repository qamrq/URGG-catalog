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
#  SCRIPT  : urgg_paths.R  -  Path, font and save helpers for R figure scripts
#  AUTHOR  : QRQ
#  INPUT   : none (URGG_PROJECT, URGG_FIGURES or config/local_paths.R)
#  OUTPUT  : urgg_in(), urgg_out(), urgg_save(), URGG_FONT
# ==========================================================

urgg_repo_root <- local({
  d <- normalizePath(if (exists(".urgg_script_dir")) .urgg_script_dir else getwd(), winslash = "/")
  repeat {
    if (file.exists(file.path(d, ".urgg-root"))) break
    parent <- dirname(d)
    if (identical(parent, d)) stop("Repository root not found (no .urgg-root file above ", getwd(), ")")
    d <- parent
  }
  d
})

local({
  cfg <- file.path(urgg_repo_root, "config", "local_paths.R")
  env <- new.env()
  if (file.exists(cfg)) sys.source(cfg, envir = env)
  pick <- function(var, default) {
    v <- Sys.getenv(var, unset = "")
    if (nzchar(v)) return(v)
    if (exists(var, envir = env, inherits = FALSE)) return(get(var, envir = env))
    default
  }
  assign("URGG_PROJECT", normalizePath(pick("URGG_PROJECT", file.path(urgg_repo_root, "data")), winslash = "/", mustWork = FALSE), envir = globalenv())
  assign("URGG_FIGURES", normalizePath(pick("URGG_FIGURES", file.path(urgg_repo_root, "figures_out")), winslash = "/", mustWork = FALSE), envir = globalenv())
})

# input file under the project folder; stops with a clear message if it is missing
urgg_in <- function(rel, must_exist = TRUE) {
  f <- file.path(URGG_PROJECT, rel)
  if (must_exist && !file.exists(f)) stop("Input not found: ", f, "\n  (set URGG_PROJECT or config/local_paths.R)", call. = FALSE)
  f
}
# output folder for one figure, created on demand
urgg_out <- function(fig) { d <- file.path(URGG_FIGURES, fig); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

# one font for all panels: URGG_FONT if set, else the first installed sans family with an italic face
URGG_FONT <- local({
  f <- Sys.getenv("URGG_FONT", unset = "")
  if (!nzchar(f)) {
    f <- "sans"
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      sf <- systemfonts::system_fonts()
      for (fam in c("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"))
        if (any(sf$family == fam & sf$italic)) { f <- fam; break }
    }
  }
  f
})

# PNG (600 dpi) + PDF (+ SVG when svglite is installed)
urgg_save <- function(p, dir, name, w_mm, h_mm) {
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm", dpi = 600, bg = "white",
                  device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm",
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  if (requireNamespace("svglite", quietly = TRUE))
    ggplot2::ggsave(file.path(dir, paste0(name, ".svg")), p, width = w_mm, height = h_mm, units = "mm", device = svglite::svglite)
  message("Written: ", file.path(dir, name), " (.png .pdf", if (requireNamespace("svglite", quietly = TRUE)) " .svg", ")")
}

message("URGG_PROJECT = ", URGG_PROJECT, "\nURGG_FIGURES = ", URGG_FIGURES, "\nfont = ", URGG_FONT)
