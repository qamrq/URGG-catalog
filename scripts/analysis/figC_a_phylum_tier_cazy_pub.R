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
#  SCRIPT  : figC_a_phylum_tier_cazy_pub.R  -  Phylum x tier x CAZy family, publication layout (Fig. 3c)
#  AUTHOR  : QRQ
#  INPUT   : figC_a_data.R (carriage_table)
#  OUTPUT  : C:/qrq_metagenome_pipeline/Plots and Figures/figC_a_phylum_tier_cazy.png (600 dpi)
#             C:/qrq_metagenome_pipeline/Plots and Figures/figC_a_phylum_tier_cazy.pdf
#             C:/qrq_metagenome_pipeline/Plots and Figures/figC_a_phylum_tier_cazy.csv
# ==========================================================

source("scripts/analysis/figC_a_data.R")

OUT_DIR <- "C:/qrq_metagenome_pipeline/Plots and Figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FAMILIES <- list(
  `Strongest core`          = c("GT19", "CE11", "GT3", "GT117"),
  `Strongest host-specific` = c("CE5", "GH112", "GH190"),
  `Within lineages`         = c("CBM98", "CE2", "GH123", "GH1", "GH38", "GT30")
)
TIER_COL <- c("Core/broad" = "#B06040", Shared = "#C89040", "Host-specific" = "#2898B8")
CMAP <- colorRampPalette(c("#F3F8FC", "#BBD8EA", "#6AA6C9", "#1F5F87"))(100)
FIG_W <- 7.2
FIG_H <- 7.8
LEFT <- 0.165; RIGHT <- 0.815; TOP <- 0.905; BOTTOM <- 0.14
COUNT_H <- 0.164      # height of the top strip (phylum names and genome counts)
GAP_W <- 0.35         # width of the gap between phyla, in column units

tab <- carriage_table(FAMILIES)
write.csv(tab$res, file.path(OUT_DIR, "figC_a_phylum_tier_cazy.csv"), row.names = FALSE)

# Columns: one per (phylum, tier); NULL marks a gap between phyla
cols <- list()
for (i in seq_along(tab$phy_levels)) {
  if (i > 1) cols <- c(cols, list(NULL))
  for (t in TIERS) cols <- c(cols, list(c(tab$phy_levels[i], t)))
}
n_cols <- length(cols)

# Rows: families, a gap row between groups, and a gap row under the last group
rows <- list()
for (i in seq_along(FAMILIES)) {
  if (i > 1) rows <- c(rows, list(NULL))
  rows <- c(rows, as.list(FAMILIES[[i]]))
}
rows <- c(rows, list(NULL))
n_rows <- length(rows)

cm <- cell_matrix(tab, rows, cols)

# x positions: real columns are 1 unit wide, gaps between phyla are narrower
XC <- rep(NA_real_, n_cols)
x_cursor <- 0
for (c in seq_len(n_cols)) {
  if (is.null(cols[[c]])) {
    x_cursor <- x_cursor + GAP_W
  } else {
    XC[c] <- x_cursor + 0.5
    x_cursor <- x_cursor + 1
  }
}
XW <- x_cursor

tier_n <- lapply(tab$phy_levels, function(p) vapply(TIERS, function(t) tier_count(tab, p, t), numeric(1)))
names(tier_n) <- tab$phy_levels

FS_LABEL <- 10
FS_FAM   <- 10
FS_GROUP <- 10
FS_CELL  <- 7
PT <- function(pt) pt / 12

draw_figure <- function() {
  heat_top <- TOP - COUNT_H - 0.012

  # ---- top strip: phylum names and genome-count boxes ----
  par(fig = c(LEFT, RIGHT, TOP - COUNT_H, TOP), mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, XW), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  for (p in tab$phy_levels) {
    cs <- which(vapply(cols, function(cc) !is.null(cc) && cc[1] == p, logical(1)))
    text(mean(XC[cs]) - 0.2, 0.52, gsub("_", " ", p), adj = c(0, 0), srt = 45,
         cex = PT(FS_LABEL), xpd = NA)
  }
  for (c in seq_len(n_cols)) {
    if (is.null(cols[[c]])) next
    n_tier <- tier_n[[cols[[c]][1]]][[cols[[c]][2]]]
    rect(XC[c] - 0.46, 0.02, XC[c] + 0.46, 0.44, col = TIER_COL[[cols[[c]][2]]],
         border = "white", lwd = 0.8)
    text(XC[c], 0.23, formatC(n_tier, big.mark = ",", format = "d"), srt = 90,
         adj = c(0.5, 0.5), cex = PT(FS_CELL), col = "white", font = 2)
  }
  text(-0.1, 0.23, "Genomes (n)", adj = c(1, 0.5), cex = PT(FS_LABEL), xpd = NA)

  # ---- heatmap: rows are families ----
  par(fig = c(LEFT, RIGHT, BOTTOM, heat_top), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, XW), ylim = c(n_rows - 0.5, -0.5), xaxs = "i", yaxs = "i")
  for (r in seq_len(n_rows)) {
    if (is.null(rows[[r]])) next
    y <- r - 1
    for (c in seq_len(n_cols)) {
      if (is.null(cols[[c]])) next
      v <- cm$mat[r, c]
      is_na <- is.na(v)
      face <- if (is_na) "#EFEFEF" else CMAP[max(1, min(100, round(v) + 1))]
      rect(XC[c] - 0.48, y - 0.48, XC[c] + 0.48, y + 0.48, col = face,
           border = if (is_na) "#9A9A9A" else "white", lwd = 0.5)
      if (is_na) {
        text(XC[c], y, "n.a.", cex = PT(FS_CELL), col = "#444444")
      } else if (nzchar(cm$lab[r, c])) {
        text(XC[c], y, cm$lab[r, c], cex = PT(FS_CELL),
             col = if (v > 70) "white" else "#1A1A1A")
      }
    }
  }
  for (r in seq_len(n_rows)) {
    if (is.null(rows[[r]])) next
    text(-0.5, r - 1, rows[[r]], adj = c(1, 0.5), font = 2, cex = PT(FS_FAM), xpd = NA)
  }
  for (grp in names(FAMILIES)) {
    rs <- which(vapply(rows, function(f) !is.null(f) && tab$fam_group[[f]] == grp, logical(1)))
    text(XW / 2, max(rs), grp, adj = c(0.5, 0.5), cex = PT(FS_GROUP), xpd = NA)
  }

  # ---- overlay: tier key and colour bar in figure coordinates ----
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend(RIGHT + 0.01, 0.80, legend = TIERS, fill = TIER_COL[TIERS], border = NA, bty = "n",
         title = "Sharing tier\n(left to right in\neach phylum)", title.adj = 0,
         cex = PT(FS_LABEL), x.intersp = 0.5, y.intersp = 1.1, xpd = NA)

  bar_x0 <- 0.33; bar_x1 <- 0.67; bar_y0 <- 0.075; bar_y1 <- 0.087
  for (k in 1:100) {
    x_a <- bar_x0 + (k - 1) / 100 * (bar_x1 - bar_x0)
    x_b <- bar_x0 + k / 100 * (bar_x1 - bar_x0)
    rect(x_a, bar_y0, x_b, bar_y1, col = CMAP[k], border = NA)
  }
  for (t in c(0, 25, 50, 75, 100)) {
    tx <- bar_x0 + t / 100 * (bar_x1 - bar_x0)
    segments(tx, bar_y0, tx, bar_y0 - 0.004, lwd = 0.6)
    text(tx, bar_y0 - 0.008, t, adj = c(0.5, 1), cex = PT(FS_LABEL - 0.5))
  }
  text((bar_x0 + bar_x1) / 2, bar_y1 + 0.012, "Genomes carrying the CAZy family (%)",
       adj = c(0.5, 0), cex = PT(FS_LABEL))
}

png(file.path(OUT_DIR, "figC_a_phylum_tier_cazy.png"), width = FIG_W, height = FIG_H,
    units = "in", res = 600, bg = "white", pointsize = 12)
draw_figure()
invisible(dev.off())

pdf(file.path(OUT_DIR, "figC_a_phylum_tier_cazy.pdf"), width = FIG_W, height = FIG_H,
    pointsize = 12)
draw_figure()
invisible(dev.off())
