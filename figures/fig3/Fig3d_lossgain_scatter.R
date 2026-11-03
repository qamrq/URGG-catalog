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
#  SCRIPT  : Fig3d_lossgain_scatter.R  -  Fig. 3d, observed vs expected loss:gain
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s80_figC_tables/figC_d_castor.csv
#  OUTPUT  : fig3/Fig3d_lossgain_scatter.{png,pdf}
# ==========================================================

.urgg_script_dir <- local({                     # folder of this script (Rscript or RStudio)
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1]))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) dirname(rstudioapi::getSourceEditorContext()$path)
  else getwd()
})
source(file.path(.urgg_script_dir, "..", "_setup", "urgg_paths.R"))

suppressPackageStartupMessages(library(ggplot2))

IN_CASTOR <- urgg_in("annotation/results_phase3/s80_figC_tables/figC_d_castor.csv")
OUT <- urgg_out("fig3")
FONT <- URGG_FONT; FS <- 12
INK  <- "#1F1F1F"; RUST <- "#A3401C"; TEAL <- "#0083A8"; GREYP <- "#9A9A9A"; RED <- "#D0312D"
theme_fig <- function(base_size = FS) theme_classic(base_size = base_size, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), plot.margin = margin(4, 8, 4, 4))

# =============================================================================
# Panel d: observed vs expected loss:gain
# =============================================================================
d <- read.csv(IN_CASTOR, stringsAsFactors = FALSE)
d$stable <- as.logical(d$stable)
fam <- d[d$type == "family", ]
cat_ <- d[d$type == "category", ]
tier_fam <- fam[fam$dir_unadj != "ns", ]

rho <- cor(tier_fam$expected_ratio, tier_fam$observed_ratio, method = "spearman")
message(sprintf("families: %d | categories: %d | unstable: %d",
                nrow(fam), nrow(cat_), sum(!d$stable)))
message(sprintf("loss > gain: %d of %d tier-associated families",
                sum(tier_fam$observed_ratio > 1), nrow(tier_fam)))
message(sprintf("Spearman rho = %.2f", rho))

rng_min <- min(d$expected_ratio, d$observed_ratio)
rng_max <- max(d$expected_ratio, d$observed_ratio)
lo_d <- 10^floor(log10(rng_min)) / 2
hi_d <- 10^ceiling(log10(rng_max)) * 2

fam$col <- ifelse(!fam$stable, RED,
                  ifelse(fam$dir_unadj == "core", RUST,
                         ifelse(fam$dir_unadj == "host", TEAL, GREYP)))

# label placement: multiplicative offsets in log space, as in the Python version
LAB <- data.frame(
  trait = c("CBM98", "GH189", "GT19", "GT117", "GH190", "Siderophores", "Extracellular proteins"),
  text  = c("CBM98", "GH189", "GT19", "GT117", "GH190", "Siderophores\n(3.5×)",
            "Extracellular proteins (unstable)"),
  fx    = c(0.35, 0.15, 0.12, 0.12, 0.12, 0.64, 2.2),
  fy    = c(6.0, 2.6, 0.9, 0.4, 0.12, 1.3, 2.0),
  hj    = c(1, 1, 1, 1, 1, 1, 0),
  stringsAsFactors = FALSE
)
LAB$ex <- sapply(LAB$trait, function(t) d$expected_ratio[d$trait == t][1])
LAB$ob <- sapply(LAB$trait, function(t) d$observed_ratio[d$trait == t][1])
LAB$stable <- sapply(LAB$trait, function(t) d$stable[d$trait == t][1])
LAB$col <- ifelse(LAB$stable, INK, RED)
LAB$lx <- LAB$ex * LAB$fx
LAB$ly <- LAB$ob * LAB$fy

p_scatter <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#666666", linewidth = 0.3) +
  geom_point(data = fam, aes(x = expected_ratio, y = observed_ratio),
             shape = 16, size = 1.8, colour = fam$col, alpha = 0.9) +
  geom_point(data = cat_, aes(x = expected_ratio, y = observed_ratio),
             shape = 18, size = 3.2,
             colour = ifelse(cat_$stable, INK, RED)) +
  geom_segment(data = LAB, aes(x = ex, y = ob, xend = lx, yend = ly), colour = "#777777",
               linewidth = 0.25) +
  geom_text(data = LAB, aes(x = lx, y = ly, label = text, hjust = hj, colour = col),
            family = FONT, size = FS / .pt, vjust = 0.5) +
  annotate("text", x = lo_d * 1.6, y = hi_d / 1.6, hjust = 0, vjust = 1, family = FONT,
           size = FS / .pt, colour = "#555555", lineheight = 1.1,
           label = "More loss-biased than\nprevalence predicts") +
  annotate("text", x = hi_d / 1.6, y = lo_d * 1.6, hjust = 1, vjust = 0, family = FONT,
           size = FS / .pt, colour = "#555555", lineheight = 1.1,
           label = "Less loss-biased than\nprevalence predicts") +
  annotate("text", x = hi_d / 1.6, y = lo_d * 3.6, hjust = 1, vjust = 0, family = FONT,
           size = FS / .pt, colour = INK, label = sprintf("Spearman \u03c1 = %.2f", rho)) +
  scale_colour_identity() +
  scale_x_log10(limits = c(lo_d, hi_d), breaks = 10^seq(log10(lo_d), log10(hi_d)),
                labels = function(x) parse(text = paste0("10^", round(log10(x))))) +
  scale_y_log10(limits = c(lo_d, hi_d), breaks = 10^seq(log10(lo_d), log10(hi_d)),
                labels = function(x) parse(text = paste0("10^", round(log10(x))))) +
  coord_fixed() +
  labs(x = "Expected loss:gain from prevalence, (1 \u2212 p)/p", y = "Observed loss:gain") +
  theme_fig()

urgg_save(p_scatter, OUT, "Fig3d_lossgain_scatter", 120, 120)
