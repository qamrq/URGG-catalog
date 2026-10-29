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
#  SCRIPT  : 31j_sunburst_ggplot_relaxed.R  -  Family to species two-ring sunburst (Fig. 1a)
#  AUTHOR  : QRQ
#  INPUT   : Cluster sharing table (CSV)
#  OUTPUT  : Sunburst figure (PDF, PNG) in OUT_DIR
# ==========================================================

FIG_R_DIR <- "scripts/figures/r"
source(file.path(FIG_R_DIR, "00_install_packages.R"))
source(file.path(FIG_R_DIR, "01_catalogue_theme.R"))
source(file.path(FIG_R_DIR, "03_catalogue_export.R"))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(cowplot)
  library(purrr); library(tidyr); library(tibble)
})

PROV_F  <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv"
OUT_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s31_mag_distribution"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. TUNABLE PARAMETERS
# =============================================================================
RELAX_EXPONENT   <- 0.35   # lower = more compression of small wedges. 1 = true proportional (off).
INNER_R0         <- 0      # inner ring inner radius (0 = donut hole starts at center)
INNER_R1         <- 1.0    # inner ring outer radius (family ring)
GAP              <- 0.03   # gap between family ring and species ring
OUTER_R0         <- INNER_R1 + GAP
OUTER_R1         <- OUTER_R0 + 0.55

# =============================================================================
# 1. LOAD DATA — family/species counts from mag_provenance.tsv
#    (falls back to the verified ground-truth table if the source file or
#    expected columns aren't found, with a loud warning)
# =============================================================================
FAMILY_MAP <- c(
  "Ovis aries" = "Bovidae", "Pantholops hodgsonii" = "Bovidae", "Bos taurus" = "Bovidae",
  "Bos grunniens" = "Bovidae", "Capra hircus" = "Bovidae", "Bubalus bubalis" = "Bovidae",
  "Bos indicus" = "Bovidae", "Procapra picticaudata" = "Bovidae", "Bos frontalis" = "Bovidae",
  "Hydropotes inermis" = "Cervidae", "Alces alces" = "Cervidae", "Cervus elaphus" = "Cervidae",
  "Rangifer tarandus" = "Cervidae", "Alces americanus" = "Cervidae",
  "Camelus bactrianus" = "Camelidae", "Camelus dromedarius" = "Camelidae",
  "Giraffa camelopardalis" = "Giraffidae"
)

load_species_counts <- function(f) {
  if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  if (!"host_species" %in% names(d)) return(NULL)
  d %>% count(host_species, name = "n") %>% rename(species = host_species)
}

species_df <- load_species_counts(PROV_F)

if (is.null(species_df)) {
  warning(
    "Could not load/derive species counts from ", PROV_F,
    " (file missing or 'host_species' column not found). Falling back to ",
    "the verified ground-truth table. VERIFY THESE ARE STILL CORRECT."
  )
  species_df <- tibble::tibble(
    species = names(FAMILY_MAP),
    n = c(6476, 4214, 3755, 3605, 2778, 2456, 950, 150, 3,   # Bovidae
          374, 296, 125, 100, 80,                             # Cervidae
          34, 3,                                               # Camelidae
          52)                                                  # Giraffidae
  )
}

species_df <- species_df %>%
  mutate(family = FAMILY_MAP[species]) %>%
  filter(!is.na(family))   # drop any species not in the known 17 -- adjust FAMILY_MAP if this fires

stopifnot("Unmapped species dropped -- check FAMILY_MAP" = nrow(species_df) == length(FAMILY_MAP))

TOTAL <- sum(species_df$n)
N_SPECIES <- nrow(species_df)
N_FAMILIES <- n_distinct(species_df$family)

# Family order fixed to match the reference image's visual weighting
# (largest family first), species ordered by descending count within family
FAMILY_ORDER <- species_df %>% count(family, wt = n, sort = TRUE) %>% pull(family)

species_df <- species_df %>%
  mutate(family = factor(family, levels = FAMILY_ORDER)) %>%
  arrange(family, desc(n)) %>%
  mutate(pct_true = 100 * n / TOTAL)

family_df <- species_df %>%
  group_by(family) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(pct_true = 100 * n / TOTAL) %>%
  arrange(match(family, FAMILY_ORDER))

cat("\n=== Ground truth (verify before trusting the figure) ===\n")
print(as.data.frame(family_df))
cat(sprintf("\nTotal genomes: %s | Families: %d | Species: %d\n",
            format(TOTAL, big.mark = ","), N_FAMILIES, N_SPECIES))
stopifnot(abs(sum(family_df$pct_true) - 100) < 1e-6)

# =============================================================================
# 2. RELAXED ANGLES -- compute visual wedge width from n^RELAX_EXPONENT,
#    renormalized to 360 degrees. Order preserved from species_df (already
#    sorted family-major, count-minor).
# =============================================================================
species_df <- species_df %>%
  mutate(relax_weight = n ^ RELAX_EXPONENT,
         relax_share = relax_weight / sum(relax_weight),
         angle_width = relax_share * 360,
         angle_end   = cumsum(angle_width),
         angle_start = angle_end - angle_width,
         angle_mid   = (angle_start + angle_end) / 2)

family_df <- species_df %>%
  group_by(family) %>%
  summarise(angle_start = min(angle_start), angle_end = max(angle_end),
            n = sum(n), .groups = "drop") %>%
  mutate(angle_mid = (angle_start + angle_end) / 2,
         pct_true = 100 * n / TOTAL) %>%
  arrange(match(family, FAMILY_ORDER))

# =============================================================================
# 3. COLOR PALETTE -- family base color, species = tint of family color
# =============================================================================
FAMILY_COLORS <- c(
  "Bovidae" = "#6FA0D6", "Cervidae" = "#F0A860", "Camelidae" = "#E38B8B", "Giraffidae" = "#8FBF7F"
)

tint <- function(hex, amount) {
  rgb_col <- grDevices::col2rgb(hex) / 255
  rgb_light <- rgb_col + (1 - rgb_col) * amount
  grDevices::rgb(rgb_light[1], rgb_light[2], rgb_light[3])
}

species_df <- species_df %>%
  group_by(family) %>%
  mutate(rank_in_family = row_number(),
         tint_amount = (rank_in_family - 1) / pmax(n() - 1, 1) * 0.55,
         fill_color = mapply(tint, FAMILY_COLORS[as.character(family)], tint_amount)) %>%
  ungroup()

family_df <- family_df %>% mutate(fill_color = FAMILY_COLORS[as.character(family)])

# =============================================================================
# 4. BUILD RING POLYGONS (angles converted to radians, 0 deg = top,
#    clockwise -- matches the reference image orientation)
# =============================================================================
deg2rad_top_cw <- function(deg) (90 - deg) * pi / 180

build_polygons <- function(df, r0, r1, group_offset = 0, n_seg = 80) {
  purrr::map_dfr(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    angs_out <- seq(row$angle_start, row$angle_end, length.out = n_seg)
    angs_in  <- rev(angs_out)
    rad_out <- deg2rad_top_cw(angs_out)
    rad_in  <- deg2rad_top_cw(angs_in)
    tibble(
      group = i + group_offset,
      x = c(r1 * cos(rad_out), r0 * cos(rad_in)),
      y = c(r1 * sin(rad_out), r0 * sin(rad_in)),
      fill_color = row$fill_color
    )
  })
}

inner_poly <- build_polygons(family_df, INNER_R0, INNER_R1, group_offset = 0)
outer_poly <- build_polygons(species_df, OUTER_R0, OUTER_R1, group_offset = 1000)

all_poly <- bind_rows(inner_poly, outer_poly)
fill_lookup <- all_poly %>% distinct(group, fill_color) %>% deframe()

# =============================================================================
# 5. SIDE LEGENDS -- family (left) and species (right). Replaces the
#    leader-line-to-wedge labels entirely.
#
#    WHY THE SWITCH: with 17 species ranging from 6,476 down to 3 genomes,
#    leader lines routed to individual wedges come out visibly UNEVEN --
#    lines to adjacent small wedges bunch together at very different
#    angles while label spacing is forced uniform, so line length/angle
#    varies a lot from label to label even though the labels themselves
#    are evenly spaced. Two flat side-legends avoid this entirely: no
#    line has to reach a specific point on the ring. Each is a clean
#    vertically-stacked list (color swatch + name + count + percentage),
#    sorted largest-to-smallest, same colors as the ring wedges so the
#    mapping to the ring stays obvious without any connecting line.
# =============================================================================
LEGEND_GAP        <- 0.45   # horizontal distance from ring edge to legend column
FAMILY_LEGEND_X   <- -(OUTER_R1 + LEGEND_GAP)
SPECIES_LEGEND_X  <-  (OUTER_R1 + LEGEND_GAP)
LEGEND_Y_TOP      <- 2.5
LEGEND_Y_BOTTOM   <- -2.5

family_legend_df <- family_df %>%
  arrange(desc(n)) %>%
  mutate(legend_x = FAMILY_LEGEND_X,
         legend_y = seq(LEGEND_Y_TOP, LEGEND_Y_BOTTOM, length.out = n())[seq_len(n())],
         label_text = sprintf("%s \u2014 %s (%s%%)", family, format(n, big.mark = ","),
                               formatC(pct_true, digits = 2, format = "f")))

species_legend_df <- species_df %>%
  arrange(desc(n)) %>%
  mutate(legend_x = SPECIES_LEGEND_X,
         legend_y = seq(LEGEND_Y_TOP, LEGEND_Y_BOTTOM, length.out = n())[seq_len(n())],
         count_pct_text = sprintf("%s (%s%%)", format(n, big.mark = ","),
                                   formatC(pct_true, digits = 2, format = "f")))

# =============================================================================
# 6. PLOT
# =============================================================================
p <- ggplot() +
  geom_polygon(data = all_poly, aes(x, y, group = group, fill = factor(group)),
               color = "white", linewidth = 0.6) +
  scale_fill_manual(values = fill_lookup, guide = "none") +
  # center summary box
  annotate("label", x = 0, y = -0.05,
           label = sprintf("Total Genomes: %s\nHost Families: %d\nHost Species: %d",
                            format(TOTAL, big.mark = ","), N_FAMILIES, N_SPECIES),
           size = 3, hjust = 0.5, vjust = 0.5, label.size = 0.3,
           fill = "white", color = "grey15") +
  # -- FAMILY legend (left side): swatch + "Family — count (pct%)" --
  geom_point(data = family_legend_df, aes(x = legend_x, y = legend_y),
             shape = 15, size = 4.2, color = family_legend_df$fill_color) +
  geom_text(data = family_legend_df,
            aes(x = legend_x + 0.16, y = legend_y, label = label_text),
            hjust = 0, size = 3.0, color = "grey15") +
  # -- SPECIES legend (right side): swatch + italic species name, count/pct
  #    on its own line beneath so long names don't crowd the numbers --
  geom_point(data = species_legend_df, aes(x = legend_x, y = legend_y),
             shape = 15, size = 3.4, color = species_legend_df$fill_color) +
  geom_text(data = species_legend_df,
            aes(x = legend_x + 0.16, y = legend_y, label = species),
            hjust = 0, fontface = "italic", size = 2.7, color = "grey15") +
  geom_text(data = species_legend_df,
            aes(x = legend_x + 1.75, y = legend_y, label = count_pct_text),
            hjust = 0, size = 2.5, color = "grey40") +
  coord_fixed(xlim = c(FAMILY_LEGEND_X - 0.3, SPECIES_LEGEND_X + 3.1),
              ylim = c(LEGEND_Y_BOTTOM - 0.2, LEGEND_Y_TOP + 0.2), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10),
        plot.caption = element_text(size = 6.5, color = "grey45", hjust = 0.5, margin = margin(t = 8))) +
  labs(caption = sprintf(
    "Wedge angles are area-relaxed (weight = count^%.2f) so low-abundance host species remain visible; printed counts and percentages are exact.",
    RELAX_EXPONENT))

# =============================================================================
# 7. EXPORT
# =============================================================================
export_catalogue_plot(p, OUT_DIR, "Fig1_sunburst_relaxed", width = 11, height = 9)

cat("\nWritten to:", OUT_DIR, "\n")
cat("Adjust RELAX_EXPONENT (currently", RELAX_EXPONENT, ") if tiny wedges are still",
    "too thin, or too exaggerated, at the top.\n")
