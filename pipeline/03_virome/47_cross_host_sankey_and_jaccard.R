###############################################################################
# cross_host_sankey_and_jaccard_ggplot.R
#
# ggplot2-native version: alluvial/Sankey-style flow diagram (ggalluvial) +
# Jaccard heatmap, combined with patchwork, exported with ggsave().
# No browser/kaleido/orca dependency - fits the same graphics pipeline as
# your other ggforce-based figures (theme_catalogue(), manual palettes, etc.)
#
# Inputs (place alongside this script):
#   sankey_flow.tsv           - species x genus_bucket x n_votus edge list
#   species_pair_jaccard.tsv  - all pairwise species Jaccard values (long format)
###############################################################################

## ---- 0. Package management -------------------------------------------------
required_pkgs <- c("ggplot2", "ggalluvial", "ggnewscale", "patchwork", "dplyr")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, quiet = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))
cat("R version:", R.version.string, "\n")

## =============================================================================
## PART 1: Data handling, manual palettes, factor ordering
## =============================================================================

flow <- read.delim("/disk252/qamar/metagenome_analysis/phase_5_virome/08_results/sankey_flow.tsv", stringsAsFactors = FALSE)
jac  <- read.delim("/disk252/qamar/metagenome_analysis/phase_5_virome/08_results/species_pair_jaccard.tsv", stringsAsFactors = FALSE)
OUT_DIR <- "/disk252/qamar/metagenome_analysis/phase_5_virome/08_results/figures/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

species_palette <- c(
  Bos_taurus             = "#1F4E8C",
  Bos_grunniens          = "#3E72B4",
  Bubalus_bubalis        = "#7FA8DA",
  Ovis_aries             = "#0F6E56",
  Capra_hircus           = "#4FA88C",
  Pantholops_hodgsonii   = "#3B6D11",
  Procapra_picticaudata  = "#8FBE5A",
  Alces_alces            = "#854F0B",
  Hydropotes_inermis     = "#D6923A",
  Giraffa_camelopardalis = "#72243E"
)

species_short <- c(
  Bos_taurus="Bos taurus", Bos_grunniens="Bos grunniens", Bubalus_bubalis="Bubalus bubalis",
  Ovis_aries="Ovis aries", Capra_hircus="Capra hircus",
  Pantholops_hodgsonii="Pantholops hodgsonii", Procapra_picticaudata="Procapra picticaudata",
  Alces_alces="Alces alces", Hydropotes_inermis="Hydropotes inermis",
  Giraffa_camelopardalis="Giraffa camelopardalis"
)
# Utility strata (not taxonomic names) - these stay upright, everything
# else (species + genus names) is drawn italic per convention.
utility_labels <- c("Other genera", "No host predicted")
# Minimum total flow for a stratum to get a text label at all (Nature style:
# label major wedges, leave minor ones to color + legend only).
MIN_LABEL_N <- 15

genus_palette <- c(
  Alistipes           = "#993C1D",
  HGM04593             = "#712B13",
  Prevotella           = "#26215C",
  Akkermansia          = "#3C3489",
  Enterousia           = "#993556",
  Tidjanibacter        = "#72243E",
  `Other genera`       = "#5F5E5A",
  `No host predicted`  = "#B4B2A9"
)

species_order <- c("Ovis_aries","Capra_hircus","Pantholops_hodgsonii","Procapra_picticaudata",
                    "Bos_taurus","Bos_grunniens","Bubalus_bubalis",
                    "Hydropotes_inermis","Alces_alces","Giraffa_camelopardalis")
species_order <- species_order[species_order %in% unique(flow$species)]

genus_order <- c("Alistipes","HGM04593","Prevotella","Akkermansia","Enterousia","Tidjanibacter",
                  "Other genera","No host predicted")
genus_order <- genus_order[genus_order %in% unique(flow$genus_bucket)]

flow$species     <- factor(flow$species, levels = species_order, labels = species_short[species_order])
flow$genus_bucket <- factor(flow$genus_bucket, levels = genus_order)

## =============================================================================
## PART 2: Build the alluvial (Sankey-style) plot + Jaccard heatmap
## =============================================================================

theme_catalogue <- function() {
  theme_minimal(base_family = "Arial", base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "plain", hjust = 0.5),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

# ---- Alluvial diagram: species (axis1) -> bacterial genus (axis2) ----
# Flows colored by target genus (the "hub" story); strata colored on their
# own scale via ggnewscale so species keep their family-shaded palette.
# Labels are placed OUTSIDE each bar (left of axis1, right of axis2), only
# for strata above MIN_LABEL_N, and italicized except for the two utility
# buckets ("Other genera" / "No host predicted").
species_totals <- flow %>% group_by(species) %>% summarise(total = sum(n_votus), .groups = "drop")
genus_totals   <- flow %>% group_by(genus_bucket) %>% summarise(total = sum(n_votus), .groups = "drop")
totals_lookup <- setNames(
  c(species_totals$total, genus_totals$total),
  c(as.character(species_totals$species), as.character(genus_totals$genus_bucket))
)
label_text <- function(stratum_chr) {
  big_enough <- totals_lookup[stratum_chr] >= MIN_LABEL_N
  big_enough[is.na(big_enough)] <- TRUE
  ifelse(big_enough, stratum_chr, "")
}
label_face <- function(stratum_chr) ifelse(stratum_chr %in% utility_labels, "plain", "italic")

fill_values_2node <- c(setNames(species_palette, species_short[species_order]), genus_palette)
fill_values_2node <- fill_values_2node[!is.na(names(fill_values_2node))]

p_alluvial <- ggplot(flow, aes(axis1 = species, axis2 = genus_bucket, y = n_votus)) +
  geom_alluvium(aes(fill = genus_bucket), width = 1/10, alpha = 0.55, color = NA) +
  scale_fill_manual(values = genus_palette, name = "Predicted genus") +
  new_scale_fill() +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/10, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = fill_values_2node, guide = "none") +
  stat_stratum(
    aes(
      label  = label_text(as.character(after_stat(stratum))),
      x      = after_stat(x) + ifelse(after_stat(x) == 1, -0.16, 0.16),
      hjust  = ifelse(after_stat(x) == 1, 1, 0),
      fontface = label_face(as.character(after_stat(stratum)))
    ),
    geom = "text", size = 3.1, color = "#2C2C2A"
  ) +
  scale_x_discrete(limits = c("Host species", "Bacterial genus"),
                    expand = expansion(mult = c(0.32, 0.32))) +
  labs(title = "Cross-host-shared vOTUs linked through predicted bacterial host genus") +
  theme_catalogue()

# ---- Jaccard heatmap, same species order ----
n <- length(species_order)
mat <- matrix(NA_real_, n, n, dimnames = list(species_order, species_order))
for (i in seq_len(nrow(jac))) {
  a <- jac$sp_a[i]; b <- jac$sp_b[i]
  if (a %in% species_order && b %in% species_order) {
    mat[a, b] <- jac$jaccard[i]; mat[b, a] <- jac$jaccard[i]
  }
}
heat_df <- as.data.frame(as.table(mat))
names(heat_df) <- c("sp_a", "sp_b", "jaccard")
heat_df$sp_a <- factor(species_short[as.character(heat_df$sp_a)], levels = rev(species_short[species_order]))
heat_df$sp_b <- factor(species_short[as.character(heat_df$sp_b)], levels = species_short[species_order])

p_heat <- ggplot(heat_df, aes(x = sp_b, y = sp_a, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradientn(colors = c("#F1EFE8", "#F0997B", "#712B13"), na.value = "grey95",
                        name = "Jaccard", limits = c(0, NA)) +
  labs(title = "Pairwise vOTU-sharing Jaccard index (all vOTUs, same species order)") +
  theme_catalogue() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "italic"))

fig <- p_alluvial / p_heat + plot_layout(heights = c(1.4, 1))

## =============================================================================
## PART 3: Export (ggsave - PDF, SVG, 600dpi PNG)
## =============================================================================


ggsave(file.path(OUT_DIR, "cross_host_sankey_and_jaccard.png"), fig, width = 12.5, height = 12, dpi = 600, bg = "white")
ggsave(file.path(OUT_DIR, "cross_host_sankey_and_jaccard.pdf"), fig, width = 12.5, height = 12, bg = "white")
ggsave(file.path(OUT_DIR, "cross_host_sankey_and_jaccard.svg"), fig, width = 12.5, height = 12, bg = "white")

#ggsave(file.path(OUT_DIR, "sankey_jancard_crossspecies.png"), fig, width = 10.5, height = 11.5, dpi = 300)


cat("Done. Files written to figures/ via ggsave (PNG 600dpi, PDF, SVG).\n")
