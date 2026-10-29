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
#  SCRIPT  : 31i_treemap_catalogue.R  -  Family to species treemap, alternate Fig. 1 layout
#  AUTHOR  : QRQ
#  INPUT   : mag_provenance.tsv (PROV_F)
#  OUTPUT  : Treemap figure (PDF, PNG) in OUT_DIR
# ==========================================================

FIG_R_DIR <- "scripts/figures/r"

source(file.path(FIG_R_DIR, "00_install_packages.R"))
source(file.path(FIG_R_DIR, "01_catalogue_theme.R"))
source(file.path(FIG_R_DIR, "02_catalogue_palette.R"))
source(file.path(FIG_R_DIR, "03_catalogue_export.R"))
source(file.path(FIG_R_DIR, "04_catalogue_labels.R"))
source(file.path(FIG_R_DIR, "05_catalogue_geometry.R"))

suppressPackageStartupMessages(library(treemapify))

PROV_F  <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv"
OUT_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s31_mag_distribution"

data <- load_catalogue_data(PROV_F)
validate_palette_coverage(data$species$species_disp)

sp_pal <- species_palette()

plot_df <- data$species %>%
  mutate(pct = 100 * n / data$total,
         label = ifelse(should_label(pct),
                        sprintf("%s\n%s\n(%.1f%%)", species_disp, format(n, big.mark=","), pct),
                        species_disp))

p <- ggplot(plot_df, aes(area = n, fill = species_disp, label = label,
                         subgroup = family)) +
  geom_treemap(color = "white", size = 2) +
  geom_treemap_subgroup_border(color = "grey15", size = 2.2) +
  geom_treemap_subgroup_text(place = "topleft", fontface = "bold",
                             colour = "grey15", alpha = .9, size = 13,
                             padding.x = grid::unit(4, "mm"),
                             padding.y = grid::unit(4, "mm")) +
  geom_treemap_text(colour = "white", place = "center", fontface = "italic",
                    size = 10, min.size = 0, reflow = TRUE) +
  scale_fill_manual(values = sp_pal, guide = "none") +
  labs(title = "Distribution of representative genomes across host species",
       subtitle = sprintf("95%% ANI species catalogue; total = %s; %d host species (grouped by family)",
                          format(data$total, big.mark = ","), data$n_species)) +
  theme_catalogue_radial()   # theme_void-based; still applies cleanly here

export_catalogue_plot(p, OUT_DIR, "Main_Fig1_treemap_catalogue",
                      width = 11, height = 7)

cat("\nDone.\n")
