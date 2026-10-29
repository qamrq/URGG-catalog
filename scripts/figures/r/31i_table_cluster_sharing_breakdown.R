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
#  SCRIPT  : 31i_table_cluster_sharing_breakdown.R  -  Host-sharing breakdown table (Fig. 1 table)
#  AUTHOR  : QRQ
#  INPUT   : Cluster sharing table (DATA_DIR)
#  OUTPUT  : Breakdown table (CSV) in OUT_DIR
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(gt)
})

# ---------------------------------------------------------------------------
# 0. PATHS — match your existing s33_cluster_sharing/ convention
# ---------------------------------------------------------------------------
DATA_DIR <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s33_cluster_sharing"
OUT_DIR  <- "/data/qrq_metagenome_pipeline/metagenome_analysis/phase_3_analysis/scripts/results_scripts/s31_mag_distribution"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SUMMARY_F <- file.path(DATA_DIR, "cluster_host_summary.csv")

# ---------------------------------------------------------------------------
# 1. COMPUTE FROM SOURCE DATA (fallback to hard-coded values below only if
#    the source file / expected columns are not found, with a loud warning
#    so it's never silently wrong)
# ---------------------------------------------------------------------------
compute_from_source <- function(f) {
  if (!file.exists(f)) return(NULL)

  d <- read_csv(f, show_col_types = FALSE)

  # Expected columns (adjust names here if your CSV uses different headers):
  #   cluster_id, n_host_species, n_genomes, sharing_category
  # sharing_category should already classify host_specific / shared / core_broad
  needed <- c("n_host_species", "n_genomes")
  if (!all(needed %in% names(d))) return(NULL)

  total <- nrow(d)

  host_specific <- d %>% filter(n_host_species == 1)
  shared        <- d %>% filter(n_host_species %in% 2:3)
  core_broad    <- d %>% filter(n_host_species > 3)

  singleton      <- host_specific %>% filter(n_genomes == 1)
  host_restricted <- host_specific %>% filter(n_genomes >= 2)

  n_hs  <- nrow(host_specific)
  n_sh  <- nrow(shared)
  n_cb  <- nrow(core_broad)
  n_sg  <- nrow(singleton)
  n_hr  <- nrow(host_restricted)

  list(
    total = total,
    host_specific = n_hs, host_specific_pct = n_hs / total,
    shared = n_sh, shared_pct = n_sh / total,
    core_broad = n_cb, core_broad_pct = n_cb / total,
    singleton = n_sg, singleton_pct_of_hs = n_sg / n_hs,
    host_restricted = n_hr, host_restricted_pct_of_hs = n_hr / n_hs
  )
}

vals <- compute_from_source(SUMMARY_F)

if (is.null(vals)) {
  warning(
    "Could not compute from ", SUMMARY_F,
    " (file missing or columns don't match n_host_species/n_genomes). ",
    "Falling back to the values supplied in chat. VERIFY THESE ARE STILL ",
    "CORRECT before using this table in the manuscript."
  )
  vals <- list(
    total = 25451,
    host_specific = 20186, host_specific_pct = 20186/25451,
    shared = 4811,         shared_pct = 4811/25451,
    core_broad = 454,      core_broad_pct = 454/25451,
    singleton = 12571,          singleton_pct_of_hs = 12571/20186,
    host_restricted = 7615,     host_restricted_pct_of_hs = 7615/20186
  )
}

# Sanity check: the three top-level categories must sum to the total, and
# singleton + host_restricted must sum to host_specific. Fail loudly if not.
stopifnot(
  "Top-level categories don't sum to total" =
    (vals$host_specific + vals$shared + vals$core_broad) == vals$total,
  "Singleton + host-restricted don't sum to host-specific" =
    (vals$singleton + vals$host_restricted) == vals$host_specific
)

pct <- function(x) sprintf("%.1f%%", x * 100)
comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# ---------------------------------------------------------------------------
# 2. BUILD TABLE DATA FRAME
#    Two-space indent prefix on sub-rows renders visually nested once gt
#    left-aligns the Category column with a monospace-safe indent.
# ---------------------------------------------------------------------------
tbl_df <- tibble::tibble(
  Category = c(
    "Host-specific (1 host species)",
    "  of which, singleton (1 genome; insufficient evidence)",
    "  of which, host-restricted (\u2265 2 genomes, evidenced)",
    "Shared (2\u20133 host species)",
    "Core / broad (> 3 host species)",
    "Total species-level clusters"
  ),
  Clusters = c(
    comma(vals$host_specific),
    comma(vals$singleton),
    comma(vals$host_restricted),
    comma(vals$shared),
    comma(vals$core_broad),
    comma(vals$total)
  ),
  `% of total` = c(
    pct(vals$host_specific_pct),
    paste0(pct(vals$singleton_pct_of_hs), " of host-specific"),
    paste0(pct(vals$host_restricted_pct_of_hs), " of host-specific"),
    pct(vals$shared_pct),
    pct(vals$core_broad_pct),
    "100%"
  ),
  is_sub   = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
  is_total = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
)

# ---------------------------------------------------------------------------
# 3. BUILD GT TABLE — Nature-style: no gridlines, thin top/bottom rules
#    only, sub-rows in grey/italic, total row bold with a top border.
#    Matches the "avoid gridlines / decorative elements" convention used
#    across every other figure in this catalogue.
# ---------------------------------------------------------------------------
COL_HEADER_BG   <- "#FFFFFF"
COL_RULE        <- "#1A1A1A"
COL_SUBROW_TEXT <- "#555555"
FONT_FAMILY     <- "Arial"

gt_tbl <- tbl_df %>%
  select(-is_sub, -is_total) %>%
  gt() %>%
  tab_header(
    title = "Table X. Distribution of species-level clusters by host-sharing breadth"
  ) %>%
  cols_label(
    Category = "Category",
    Clusters = "Clusters",
    `% of total` = "% of total"
  ) %>%
  cols_align(align = "left",  columns = Category) %>%
  cols_align(align = "right", columns = c(Clusters, `% of total`)) %>%
  tab_options(
    table.font.names = FONT_FAMILY,
    table.font.size = px(11),
    heading.title.font.size = px(12),
    heading.title.font.weight = "bold",
    column_labels.font.weight = "bold",
    column_labels.border.top.color = COL_RULE,
    column_labels.border.top.width = px(1.5),
    column_labels.border.bottom.color = COL_RULE,
    column_labels.border.bottom.width = px(1),
    table.border.top.style = "hidden",
    table.border.bottom.color = COL_RULE,
    table.border.bottom.width = px(1.5),
    table_body.hlines.style = "hidden",
    table_body.border.bottom.style = "hidden",
    data_row.padding = px(4),
    column_labels.padding = px(4)
  ) %>%
  # Sub-rows: indented text already baked into the string; style grey + italic
  tab_style(
    style = list(cell_text(color = COL_SUBROW_TEXT, style = "italic", size = px(10))),
    locations = cells_body(rows = which(tbl_df$is_sub))
  ) %>%
  # Total row: bold with a rule above it
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_borders(sides = "top", color = COL_RULE, weight = px(1))
    ),
    locations = cells_body(rows = which(tbl_df$is_total))
  ) %>%
  tab_source_note(
    source_note = md(sprintf(
      "*Species-level clusters (95%% ANI), n = %s. Singleton clusters (single-genome, single-host) are shown separately from host-restricted clusters (\u2265 2 genomes, all one host) as the former lack the multi-genome evidence needed to confirm true host restriction.*",
      comma(vals$total)
    ))
  )

# ---------------------------------------------------------------------------
# 4. EXPORT — PNG (300 dpi, matches Nature Microbiology's stated minimum)
#    and PDF (vector, for direct placement without any rasterization).
#    Width kept narrow (~89 mm single-column) since this is a compact
#    summary table meant to sit alongside other Fig 1 panels, not span
#    the full 180 mm page width.
# ---------------------------------------------------------------------------
out_png <- file.path(OUT_DIR, "Fig1_table_cluster_sharing.png")
out_pdf <- file.path(OUT_DIR, "Fig1_table_cluster_sharing.pdf")

gtsave(gt_tbl, out_png, vwidth = 900, vheight = 500, zoom = 3)  # zoom=3 -> ~300dpi equivalent at this vwidth
gtsave(gt_tbl, out_pdf)

message("Table written to:\n  ", out_png, "\n  ", out_pdf)
message("\nComputed values used:")
message(sprintf("  Host-specific: %s (%s)", comma(vals$host_specific), pct(vals$host_specific_pct)))
message(sprintf("    Singleton:        %s (%s of host-specific)", comma(vals$singleton), pct(vals$singleton_pct_of_hs)))
message(sprintf("    Host-restricted:  %s (%s of host-specific)", comma(vals$host_restricted), pct(vals$host_restricted_pct_of_hs)))
message(sprintf("  Shared:        %s (%s)", comma(vals$shared), pct(vals$shared_pct)))
message(sprintf("  Core/broad:    %s (%s)", comma(vals$core_broad), pct(vals$core_broad_pct)))
message(sprintf("  Total:         %s", comma(vals$total)))
