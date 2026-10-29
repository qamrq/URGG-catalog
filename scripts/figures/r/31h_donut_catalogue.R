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
#  SCRIPT  : 31h_donut_catalogue.R  -  Species donut, alternate Fig. 1 layout
#  AUTHOR  : QRQ
#  INPUT   : mag_provenance.tsv (PROV_F)
#  OUTPUT  : Donut figure (PDF, SVG, PNG) in OUT_DIR
# ==========================================================

FIG_R_DIR <- "scripts/figures/r"

source(file.path(FIG_R_DIR, "00_install_packages.R"))
source(file.path(FIG_R_DIR, "01_catalogue_theme.R"))
source(file.path(FIG_R_DIR, "02_catalogue_palette.R"))
source(file.path(FIG_R_DIR, "03_catalogue_export.R"))
source(file.path(FIG_R_DIR, "04_catalogue_labels.R"))
source(file.path(FIG_R_DIR, "05_catalogue_geometry.R"))

suppressPackageStartupMessages({ library(ggforce); library(cowplot) })

PROV_F  <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv"
OUT_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s31_mag_distribution"

data <- load_catalogue_data(PROV_F)
validate_palette_coverage(data$species$species_disp)

# Single WIDE ring (point 5): thick enough to read as a proper donut, not
# a thin coord_polar-style bar.
sp_geo <- compute_ring_geometry(data$species, "n", r0 = 1.3, r1 = 2.6)
sp_pal <- species_palette()
sp_geo$fill_col <- sp_pal[sp_geo$species_disp]

p <- ggplot() +
  geom_arc_bar(data = sp_geo,
              aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end,
                  fill = species_disp),
              color = "white", linewidth = 0.9) +
  scale_fill_manual(values = sp_pal, guide = "none") +
  coord_fixed() +
  theme_catalogue_radial() +
  labs(title = "Distribution of representative genomes across host species",
       subtitle = sprintf("95%% ANI species catalogue; total = %s; %d host species",
                          format(data$total, big.mark = ","), data$n_species))

sp_major <- sp_geo %>% filter(should_label(pct))
p <- p + geom_text(data = sp_major,
                   aes(x = label_x_in, y = label_y_in,
                       label = sprintf("%s\n%s (%.1f%%)",
                                      species_disp, format(n, big.mark = ","), pct)),
                   color = "white", fontface = "bold", size = 3.0, lineheight = .9)

p <- p + build_center_annotation(
  total = data$total,
  subtitle_lines = c("species representatives", "95% ANI catalogue",
                     sprintf("%d host species", data$n_species))
)

legend_df <- data$species %>%
  mutate(species_disp = factor(species_disp, levels = rev(species_disp)),
         txt = sprintf("%s  (%.1f%%)", format(n, big.mark = ","), pct))

p_legend <- ggplot(legend_df, aes(n, species_disp, fill = species_disp)) +
  geom_col(width = .6, show.legend = FALSE) +
  geom_text(aes(label = txt), hjust = -0.08, size = 2.8, colour = "grey20") +
  scale_fill_manual(values = sp_pal) +
  scale_x_continuous(expand = expansion(mult = c(0, .55))) +
  theme_catalogue_legend()

final_plot <- plot_grid(p, p_legend, ncol = 2, rel_widths = c(1.3, 1))
export_catalogue_plot(final_plot, OUT_DIR, "Main_Fig1_donut_catalogue",
                      width = 13, height = 8)

cat("\nDone.\n")
