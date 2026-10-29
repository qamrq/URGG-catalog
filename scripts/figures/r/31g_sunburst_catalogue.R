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
#  SCRIPT  : 31g_sunburst_catalogue.R  -  Family to species sunburst (Fig. 1a)
#  AUTHOR  : QRQ
#  INPUT   : mag_provenance.tsv (PROV_F)
#  OUTPUT  : Sunburst figure (PDF, SVG, PNG) in OUT_DIR
# ==========================================================

# Set this to wherever the R/ framework directory lives on your server.
# Explicit and robust, rather than trying to auto-detect a script's own
# path (which does not resolve reliably when run via `Rscript` — see
# 00_install_packages.R's CATALOGUE_R_DIR note).
CATALOGUE_R_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_3_analysis/scripts/results_scripts/R"

source(file.path(CATALOGUE_R_DIR, "00_install_packages.R"))
source(file.path(CATALOGUE_R_DIR, "01_catalogue_theme.R"))
source(file.path(CATALOGUE_R_DIR, "02_catalogue_palette.R"))
source(file.path(CATALOGUE_R_DIR, "03_catalogue_export.R"))
source(file.path(CATALOGUE_R_DIR, "04_catalogue_labels.R"))
source(file.path(CATALOGUE_R_DIR, "05_catalogue_geometry.R"))

suppressPackageStartupMessages({ library(ggforce); library(cowplot) })

PROV_F  <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv"
OUT_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s31_mag_distribution"

# ------------------------------------------------------------ load data ----
data <- load_catalogue_data(PROV_F)
validate_palette_coverage(data$species$species_disp)

cat("\n=== Family totals ===\n")
print(as.data.frame(data$family %>% mutate(pct = round(pct, 2))))
cat(sprintf("\nTotal representatives: %s | Host species: %d\n",
            format(data$total, big.mark = ","), data$n_species))

# ------------------------------------------------------------- geometry ----
# WIDE ring spacing (plotting_instructions.md point 5): inner (family) ring
# thick, small gap, outer (species) ring thicker still.
fam_geo <- compute_ring_geometry(data$family, "n", r0 = 1.0, r1 = 1.9)
sp_geo  <- compute_ring_geometry(data$species, "n", r0 = 1.95, r1 = 3.3)
validate_geometry_totals(fam_geo, sp_geo)

fam_pal <- family_palette()
sp_pal  <- species_palette()

fam_geo$fill_col <- fam_pal[fam_geo$family]
sp_geo$fill_col  <- sp_pal[sp_geo$species_disp]

# --------------------------------------------------------------- plot -----
p <- ggplot() +
  # Family ring
  geom_arc_bar(data = fam_geo,
              aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end,
                  fill = family),
              color = "white", linewidth = 1.1) +
  # Species ring
  geom_arc_bar(data = sp_geo,
              aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end,
                  fill = species_disp),
              color = "white", linewidth = 0.7) +
  scale_fill_manual(values = c(fam_pal, sp_pal), guide = "none") +
  coord_fixed() +
  theme_catalogue_radial() +
  labs(title = "Distribution of representative genomes across host species",
       subtitle = sprintf("95%% ANI species catalogue; total = %s; %d host species",
                          format(data$total, big.mark = ","), data$n_species))

# Family labels — centred in the inner ring, always shown (only 4 families).
p <- p + geom_text(data = fam_geo,
                   aes(x = label_x_in, y = label_y_in,
                       label = sprintf("%s\n%.1f%%", family, pct)),
                   color = "white", fontface = "bold", size = 3.4, lineheight = .9)

# Species labels — ONLY major wedges on-plot (Nature style: minor taxa go
# to the legend only). Threshold set in 04_catalogue_labels.R.
sp_major <- sp_geo %>% filter(should_label(pct))
p <- p + geom_text(data = sp_major,
                   aes(x = label_x_in, y = label_y_in,
                       label = sprintf("%s\n%s (%.1f%%)",
                                      species_disp, format(n, big.mark = ","), pct)),
                   color = "grey8", fontface = "italic", size = 2.9, lineheight = .9)

# Centre annotation
p <- p + build_center_annotation(
  total = data$total,
  subtitle_lines = c("species representatives", "95% ANI catalogue",
                     sprintf("%d host species", data$n_species))
)

# ------------------------------------------------------- legend panel ------
# Every species listed individually (including those below the on-plot
# label threshold), with proportional bar + count + percentage. Grouped by
# family with family sub-headers so the legend mirrors the ring structure.
legend_df <- data$species %>%
  mutate(species_disp = factor(species_disp, levels = rev(species_disp)),
         txt = sprintf("%s  (%.1f%%)", format(n, big.mark = ","), pct))

p_legend <- ggplot(legend_df, aes(n, species_disp, fill = species_disp)) +
  geom_col(width = .6, show.legend = FALSE) +
  geom_text(aes(label = txt), hjust = -0.08, size = 2.8, colour = "grey20") +
  scale_fill_manual(values = sp_pal) +
  scale_x_continuous(expand = expansion(mult = c(0, .55))) +
  theme_catalogue_legend()

# --------------------------------------------------------- assemble --------
final_plot <- plot_grid(p, p_legend, ncol = 2, rel_widths = c(1.3, 1))

export_catalogue_plot(final_plot, OUT_DIR, "Main_Fig1_sunburst_catalogue",
                      width = 13, height = 8)

cat("\nDone.\n")
