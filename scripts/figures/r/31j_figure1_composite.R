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
#  SCRIPT  : 31j_figure1_composite.R  -  Figure 1 five-panel composite
#  AUTHOR  : QRQ
#  INPUT   : mag_provenance.tsv (PROV_F) and quality and funnel tables
#  OUTPUT  : Figure 1 composite (PDF, PNG) in OUT_DIR
# ==========================================================

CATALOGUE_R_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_3_analysis/scripts/results_scripts/figures_script"

source(file.path(CATALOGUE_R_DIR, "00_install_packages.R"))
source(file.path(CATALOGUE_R_DIR, "01_catalogue_theme.R"))
source(file.path(CATALOGUE_R_DIR, "02_catalogue_palette.R"))
source(file.path(CATALOGUE_R_DIR, "03_catalogue_export.R"))
source(file.path(CATALOGUE_R_DIR, "04_catalogue_labels.R"))
source(file.path(CATALOGUE_R_DIR, "05_catalogue_geometry.R"))

suppressPackageStartupMessages({ library(ggforce); library(cowplot); library(tidyr) })

if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required for panel (c). install.packages('magick')")
}
library(magick)

PROV_F  <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv"
OUT_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s31_mag_distribution"
PANEL_C_IMAGE <- file.path(OUT_DIR, "Main_Fig2_catalogue_flow.png")

if (!file.exists(PANEL_C_IMAGE)) {
  stop("Panel (c) image not found at ", PANEL_C_IMAGE,
       "\nRun 31c_mag_distribution_extra_figures.R first to generate it.")
}

data <- load_catalogue_data(PROV_F)
validate_palette_coverage(data$species$species_disp)
prov <- read_tsv(PROV_F, show_col_types = FALSE)

fam_pal <- family_palette()
sp_pal  <- species_palette()

# =============================================================================
# PANEL (a) — Sunburst (centre, the hero panel)
# =============================================================================
fam_geo <- compute_ring_geometry(data$family, "n", r0 = 1.0, r1 = 1.9)
sp_geo  <- compute_ring_geometry(data$species, "n", r0 = 1.95, r1 = 3.3)
validate_geometry_totals(fam_geo, sp_geo)
fam_geo$fill_col <- fam_pal[fam_geo$family]
sp_geo$fill_col  <- sp_pal[sp_geo$species_disp]

p_a <- ggplot() +
  geom_arc_bar(data = fam_geo, aes(x0=0, y0=0, r0=r0, r=r1, start=start, end=end, fill=family),
              color = "white", linewidth = 1.1) +
  geom_arc_bar(data = sp_geo, aes(x0=0, y0=0, r0=r0, r=r1, start=start, end=end, fill=species_disp),
              color = "white", linewidth = 0.7) +
  scale_fill_manual(values = c(fam_pal, sp_pal), guide = "none") +
  coord_fixed(clip = "off") +
  theme_catalogue_radial() +
  labs(title = "a", subtitle = NULL)

fam_geo$wedge_width <- fam_geo$end - fam_geo$start
fam_inside <- fam_geo %>% filter(wedge_width >= 0.35)
sp_major   <- sp_geo  %>% filter(should_label(pct))

p_a <- p_a +
  geom_text(data = fam_inside, aes(x = label_x_in, y = label_y_in,
                                   label = sprintf("%s\n%.1f%%", family, pct)),
           color = "white", fontface = "bold", size = 3.2, lineheight = .9) +
  geom_text(data = sp_major, aes(x = label_x_in, y = label_y_in,
                                 label = sprintf("%s\n%.1f%%", species_disp, pct)),
           color = "grey8", fontface = "italic", size = 2.6, lineheight = .9) +
  build_center_annotation(data$total,
                          c("representatives", "95% ANI catalogue",
                            sprintf("%d host species", data$n_species)),
                          size_total = 6.5, size_sub = 2.6)

# =============================================================================
# PANEL (b) — Completeness vs. contamination (all QC-passed MAGs)
# =============================================================================
all_mags <- prov %>%
  distinct(genome, completeness, contamination) %>%
  filter(!is.na(completeness), !is.na(contamination)) %>%
  mutate(tier = case_when(
    completeness > 90 & contamination <= 5                       ~ "High-quality",
    completeness >= 70 & completeness <= 90 & contamination <= 5 ~ "Medium-quality",
    TRUE                                                          ~ "Outside QC"
  ))

p_b <- ggplot(all_mags, aes(completeness, contamination, color = tier)) +
  geom_point(size = .3, alpha = .2) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey40", linewidth = .3) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey40", linewidth = .3) +
  scale_color_manual(values = c("High-quality" = "#D62728", "Medium-quality" = "#7F7F7F",
                                "Outside QC" = "#BF8F00"), guide = "none") +
  scale_y_continuous(limits = c(0, 7)) +
  labs(title = "b", x = "Completeness (%)", y = "Contamination (%)") +
  theme_catalogue_cartesian(base_size = 8)

# =============================================================================
# PANEL (c) — Catalogue flow (imported, ggalluvial-based; see header note)
# =============================================================================
img_c <- image_read(PANEL_C_IMAGE)
p_c <- ggdraw() + draw_image(img_c) +
  draw_label("c", x = 0.02, y = 0.97, fontface = "bold", size = 14, hjust = 0)

# =============================================================================
# PANEL (d) — Quality tier breakdown (final representatives)
# =============================================================================
final_reps <- prov %>%
  filter(genome == final_representative) %>%
  distinct(genome, completeness, contamination) %>%
  mutate(quality_tier = case_when(
    completeness > 90 & contamination <= 5                       ~ "High-quality\n(HQ)",
    completeness >= 70 & completeness <= 90 & contamination <= 5 ~ "Medium-quality\n(MQ)",
    TRUE ~ "Other"
  ))

qtier_tbl <- final_reps %>%
  count(quality_tier, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))

p_d <- ggplot(qtier_tbl, aes(reorder(quality_tier, n), n, fill = quality_tier)) +
  geom_col(width = .6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%s (%.1f%%)", format(n, big.mark=","), pct)),
            vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c("#2F5496", "#3C7A3C", "#BF8F00")) +
  scale_y_continuous(expand = expansion(mult = c(0, .22))) +
  labs(title = "d", x = NULL, y = "Representative genomes") +
  theme_catalogue_cartesian(base_size = 8) +
  theme(axis.text.x = element_text(size = 7))

# =============================================================================
# PANEL (e) — Genome distribution by GI sample site
# =============================================================================
if ("sample_site" %in% names(prov)) {
  site_tbl <- prov %>%
    filter(!is.na(sample_site), trimws(sample_site) != "") %>%
    count(sample_site, name = "n") %>%
    mutate(pct = round(100 * n / sum(n), 1))

  p_e <- ggplot(site_tbl, aes(pct, reorder(sample_site, pct))) +
    geom_segment(aes(x = 0, xend = pct, y = reorder(sample_site, pct),
                     yend = reorder(sample_site, pct)),
                color = "grey70", linewidth = .6) +
    geom_point(size = 2.2, color = "#3C7A3C") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.4, size = 2.4) +
    scale_x_continuous(expand = expansion(mult = c(0, .3))) +
    labs(title = "e", x = "% of genomes", y = NULL) +
    theme_catalogue_cartesian(base_size = 8) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.5))
} else {
  p_e <- ggplot() + labs(title = "e", subtitle = "sample_site not available") + theme_void()
}

# =============================================================================
# ASSEMBLE — sunburst centred, four panels at the corners
# =============================================================================
composite <- ggdraw() +
  draw_plot(p_a, x = 0.24, y = 0.22, width = 0.52, height = 0.56) +   # centre
  draw_plot(p_b, x = 0.00, y = 0.55, width = 0.26, height = 0.42) +   # top-left
  draw_plot(p_c, x = 0.74, y = 0.55, width = 0.26, height = 0.42) +   # top-right
  draw_plot(p_d, x = 0.00, y = 0.03, width = 0.26, height = 0.42) +   # bottom-left
  draw_plot(p_e, x = 0.74, y = 0.03, width = 0.26, height = 0.42) +   # bottom-right
  draw_label("Figure 1. Ruminant gastrointestinal MAG catalogue overview",
            x = 0.5, y = 0.99, fontface = "bold", size = 14, hjust = 0.5)

export_catalogue_plot(composite, OUT_DIR, "Main_Fig1_composite",
                      width = 16, height = 12)

cat("\nNOTE: panel layout coordinates (x/y/width/height in draw_plot() above)\n")
cat("are a starting point — inspect the exported PNG and nudge them if any\n")
cat("panel overlaps the sunburst or the outer figure margin.\n")
