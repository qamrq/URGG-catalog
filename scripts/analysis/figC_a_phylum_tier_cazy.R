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
#  SCRIPT  : figC_a_phylum_tier_cazy.R  -  Phylum x tier x CAZy family, working version (Fig. 3c)
#  AUTHOR  : QRQ
#  INPUT   : figC_a_data.R (carriage_table)
#  OUTPUT  : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/scrap/figC_a_phylum_tier_cazy.png (600 dpi)
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/scrap/figC_a_phylum_tier_cazy.pdf
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/scrap/figC_a_phylum_tier_cazy.csv
# ==========================================================

source("scripts/analysis/figC_a_data.R")

OUT_DIR <- "C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/scrap"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FAMILIES <- list(
  `Strongest core`          = c("GT19", "CE11", "GT3", "GT117"),
  `Strongest host-specific` = c("CE5", "GH112", "GH190"),
  `Within lineages`         = c("CBM98", "CE2", "GH123", "GH1", "GH38", "GT30"),
  `Loss:gain outlier`       = "GH189"
)
TIER_COL <- c("Core/broad" = "#C0392B", Shared = "#E8A33D", "Host-specific" = "#1B7F9E")
CMAP <- colorRampPalette(c("#F5F3EE", "#8FBFC4", "#1B6B76", "#0B3D45"))(100)
FIG_W <- 18
FIG_H <- 14.5

tab <- carriage_table(FAMILIES)
write.csv(tab$res, file.path(OUT_DIR, "figC_a_phylum_tier_cazy.csv"), row.names = FALSE)

# Columns: one per (phylum, tier); NULL marks a gap between phyla
cols <- list()
for (i in seq_along(tab$phy_levels)) {
  if (i > 1) cols <- c(cols, list(NULL))
  for (t in TIERS) cols <- c(cols, list(c(tab$phy_levels[i], t)))
}
n_cols <- length(cols)

# Rows: families, with a gap row between groups
rows <- list()
for (i in seq_along(FAMILIES)) {
  if (i > 1) rows <- c(rows, list(NULL))
  rows <- c(rows, as.list(FAMILIES[[i]]))
}
n_rows <- length(rows)

cm <- cell_matrix(tab, rows, cols)
top_n <- vapply(cols, function(cc) if (is.null(cc)) NA_real_ else
                  as.numeric(tier_count(tab, cc[1], cc[2])), numeric(1))
fmt_n <- function(v) if (v >= 1000) sprintf("%.1fk", v / 1000) else as.character(as.integer(v))

# Positions are 0-based: column c sits at x = c - 1, row r at y = r - 1
draw_figure <- function() {
  # ---- legend (tier colours) ----
  par(fig = c(0.12, 0.84, 0.935, 0.99), mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = TIERS, fill = TIER_COL[TIERS], border = NA, bty = "n",
         ncol = 3, cex = 13 / 12, x.intersp = 0.8,
         title = "Sharing tier (columns in each phylum, left to right)")

  # ---- genome-count row ----
  par(fig = c(0.12, 0.84, 0.845, 0.90), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(-0.5, n_cols - 0.5), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  for (c in seq_along(cols)) {
    if (is.null(cols[[c]])) next
    x <- c - 1
    rect(x - 0.45, 0.05, x + 0.45, 0.95, col = TIER_COL[[cols[[c]][2]]], border = "white", lwd = 1.1)
    text(x, 0.5, fmt_n(top_n[c]), col = "white", font = 2, cex = 9.5 / 12)
  }
  text(-0.7, 0.5, "Genomes (n)", adj = c(1, 0.5), cex = 13 / 12, xpd = NA)

  # ---- heatmap ----
  par(fig = c(0.12, 0.84, 0.10, 0.84), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(-0.5, n_cols - 0.5), ylim = c(n_rows - 0.5, -0.5), xaxs = "i", yaxs = "i")
  for (r in seq_len(n_rows)) {
    if (is.null(rows[[r]])) next
    y <- r - 1
    for (c in seq_len(n_cols)) {
      if (is.null(cols[[c]])) next
      x <- c - 1
      v <- cm$mat[r, c]
      face <- if (is.na(v)) "#FFFFFF" else CMAP[max(1, min(100, round(v) + 1))]
      rect(x - 0.45, y - 0.45, x + 0.45, y + 0.45, col = face, border = "white", lwd = 1.1)
      if (nzchar(cm$lab[r, c])) {
        text(x, y, cm$lab[r, c], cex = 12 / 12,
             col = if (!is.na(v) && v > 55) "white" else "#222222")
      }
    }
  }
  for (r in seq_len(n_rows)) {
    if (is.null(rows[[r]])) next
    text(-0.75, r - 1, rows[[r]], adj = c(1, 0.5), font = 2, cex = 14 / 12,
         col = "#333333", xpd = NA)
  }
  for (grp in names(FAMILIES)) {
    rs <- which(vapply(rows, function(f) !is.null(f) && tab$fam_group[[f]] == grp, logical(1)))
    text(n_cols - 0.3, mean(rs) - 1, grp, adj = c(0, 0.5), cex = 13 / 12,
         col = "#333333", xpd = NA)
  }

  # ---- phylum labels above the count row ----
  for (p in tab$phy_levels) {
    cs <- which(vapply(cols, function(cc) !is.null(cc) && cc[1] == p, logical(1)))
    label <- paste(strwrap(gsub("_", " ", p), width = 20), collapse = "\n")
    text(mean(cs) - 1, 1.2, label, adj = c(0.5, 0), font = 2, cex = 11 / 12, xpd = NA)
  }

  # ---- colour bar ----
  par(fig = c(0.36, 0.64, 0.02, 0.05), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 100), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  for (k in 1:100) rect(k - 1, 0, k, 1, col = CMAP[k], border = NA)
  axis(1, at = c(0, 25, 50, 75, 100), cex.axis = 12 / 12, lwd = 0, tcl = -0.3, line = -0.5)
  mtext("Genomes carrying the CAZy family (%)", side = 3, line = 0.4, cex = 13 / 12)
}

png(file.path(OUT_DIR, "figC_a_phylum_tier_cazy.png"), width = FIG_W, height = FIG_H,
    units = "in", res = 600, bg = "white", pointsize = 12)
draw_figure()
invisible(dev.off())

pdf(file.path(OUT_DIR, "figC_a_phylum_tier_cazy.pdf"), width = FIG_W, height = FIG_H,
    pointsize = 12)
draw_figure()
invisible(dev.off())
