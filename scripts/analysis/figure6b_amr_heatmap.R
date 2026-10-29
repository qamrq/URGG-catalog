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
#  SCRIPT  : figure6b_amr_heatmap.R  -  AMR gene families across host species (Fig. 6b)
#  AUTHOR  : QRQ
#  INPUT   : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig6/result_files/amr_heatmap_matrix_reps.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig6/result_files/amr_family_meta.tsv
#  OUTPUT  : C:/qrq_metagenome_pipeline/Plots and Figures/figure6b_amr_heatmap.png (600 dpi)
# ==========================================================
#
# Rows    : AMR gene families (AMRFinderPlus, collapsed by name), clustered
# Columns : host species with at least MAG_SPLIT MAGs, clustered
# Fill    : log10(hits per 1,000 representative genomes); grey = no hit
# Strips  : drug class | resistance mechanism | gene-family sharing tier

TABLE_DIR <- "C:/qrq_metagenome_pipeline/manuscript/Figures/Fig6/result_files"
OUT_DIR   <- "C:/qrq_metagenome_pipeline/Plots and Figures"
OUT_FILE  <- file.path(OUT_DIR, "figure6b_amr_heatmap.png")
DPI       <- 600
MAG_SPLIT <- 1000   # hosts with at least this many MAGs are shown

CANVAS_W <- 2300    # canvas units, 100 per inch
CANVAS_H <- 1520

# Font sizes in points
FS_ROW   <- 24
FS_COL   <- 24
FS_AXIS  <- 27
FS_STRIP <- 21
FS_LEG_T <- 24
FS_LEG   <- 21
FS_TICK  <- 20

RENAME <- c(
  "AAC(6')-I family aminoglycoside 6'-N-acetyltransferase" = "AAC(6')-I acetyltransferase",
  "aminoglycoside 6'-N-acetyltransferase" = "AAC(6')-N-acetyltransferase",
  "ANT(9) family aminoglycoside nucleotidyltransferase" = "ANT(9) nucleotidyltransferase",
  "CfxA family broad-spectrum class A beta-lactamase" = "CfxA family beta-lactamase",
  "tetracycline resistance ribosomal protection protein" = "tetracycline ribosomal protection protein",
  "23S rRNA (adenine(2058)-N(6))-methyltransferase" = "23S rRNA methyltransferase",
  "trimethoprim-resistant dihydrofolate reductase" = "trimethoprim-resistant DHFR"
)

BINOM <- c(
  Ovis_aries = "Ovis aries", Bos_taurus = "Bos taurus", Bos_grunniens = "Bos grunniens",
  Pantholops_hodgsonii = "Pantholops hodgsonii", Capra_hircus = "Capra hircus",
  Bubalus_bubalis = "Bubalus bubalis", Bos_indicus = "Bos indicus",
  Hydropotes_inermis = "Hydropotes inermis", Alces_alces = "Alces alces",
  Procapra_picticaudata = "Procapra picticaudata", Cervus_elaphus = "Cervus elaphus",
  Giraffa_camelopardalis = "Giraffa camelopardalis", Rangifer_tarandus = "Rangifer tarandus",
  Alces_americanus = "Alces americanus", Camelus_bactrianus = "Camelus bactrianus",
  Camelus_dromedarius = "Camelus dromedarius", Bos_frontalis = "Bos frontalis"
)

CLASS_PAL <- c(Tetracycline = "#D85A30", Aminoglycoside = "#0F6E56", Lincosamide = "#93A8D6",
               `Lincosamide/macrolide/streptogramin` = "#3C3489", Macrolide = "#F0997B",
               `Beta-lactam` = "#854F0B", Nitroimidazole = "#26215C", Trimethoprim = "#D6923A")
MECH_PAL  <- c(`Antibiotic efflux` = "#1F2A5C", `Antibiotic inactivation` = "#D85A30",
               `Antibiotic target alteration` = "#5FA83C", `Antibiotic target protection` = "#4FA8B4")
TIER_PAL  <- c(`Host-specific` = "#2898B8", Shared = "#C89040", `Core/broad` = "#B06040")
TIER_LABEL <- c(`Host-specific` = "Host specific", `Core/broad` = "Core")
ZERO_COL  <- "#B4B2A9"
RAMP      <- colorRampPalette(c("#3B2D96", "#F1EFE8", "#CB519B"))(100)

sentence_case <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

# ---- Data ------------------------------------------------------------------
mat  <- read.delim(file.path(TABLE_DIR, "amr_heatmap_matrix_reps.tsv"), stringsAsFactors = FALSE)
meta <- read.delim(file.path(TABLE_DIR, "amr_family_meta.tsv"), stringsAsFactors = FALSE)

shorten <- function(x) ifelse(x %in% names(RENAME), RENAME[x], x)
mat$fam  <- shorten(mat$family)
meta$fam <- shorten(meta$family)

# Sharing tier from the number of host species with a hit, across all hosts
hit <- mat[mat$n_hits > 0, ]
n_species <- tapply(hit$host, hit$fam, function(x) length(unique(x)))
meta$n_species <- as.integer(ifelse(is.na(n_species[meta$fam]), 0, n_species[meta$fam]))
meta$tier <- ifelse(meta$n_species <= 1, "Host-specific",
              ifelse(meta$n_species <= 4, "Shared", "Core/broad"))
meta$class_lab <- sentence_case(meta$class_)
meta$mech_lab  <- sentence_case(meta$mechanism)

total_mags <- unique(mat[, c("host", "total_mags")])
hosts_sel  <- total_mags$host[total_mags$total_mags >= MAG_SPLIT]
sub <- mat[mat$host %in% hosts_sel, ]

# ---- Clustering on log rate (zero cells = log10(0.01)) ----------------------
fam_lv <- sort(unique(sub$fam))
W <- matrix(log10(0.01), nrow = length(fam_lv), ncol = length(hosts_sel),
            dimnames = list(fam_lv, hosts_sel))
W[cbind(match(sub$fam, fam_lv), match(sub$host, hosts_sel))] <- log10(sub$rate_per_1000 + 0.01)

hc_row <- hclust(dist(W), method = "average")
hc_col <- hclust(dist(t(W)), method = "average")
row_order <- fam_lv[hc_row$order]
col_order <- hosts_sel[hc_col$order]
nR <- length(row_order)
nC <- length(col_order)

# Cell matrices (family x host), zero where absent
pivot <- function(value) {
  p <- tapply(sub[[value]], list(sub$fam, sub$host), sum)
  p <- p[row_order, col_order, drop = FALSE]
  p[is.na(p)] <- 0
  p
}
N_HITS <- pivot("n_hits")
RATE   <- pivot("rate_per_1000")
has    <- N_HITS > 0
logr   <- log10(RATE[has])
r0 <- min(logr)
r1 <- max(logr)
if (r1 == r0) { r0 <- r0 - 0.1; r1 <- r1 + 0.1 }
idx <- pmin(pmax(round(1 + 99 * (log10(pmax(RATE, 1e-9)) - r0) / (r1 - r0)), 1), 100)
CELL_COL <- ifelse(has, RAMP[idx], ZERO_COL)

# ---- Geometry (canvas units, y axis points down) ----------------------------
x0 <- 760; cw <- 112
y0 <- 215; rh <- 46
x1 <- x0 + nC * cw
y1 <- y0 + nR * rh
strip_x <- c(24, 60, 96); strip_w <- 30
d_top <- 60; d_bot <- y0 - 22
lx <- 1740

# Column dendrogram: leaf i sits at the centre of column i
leaf_x <- x0 + (seq_len(nC) - 0.5) * cw
height_y <- function(h, hmax) d_bot - h / hmax * (d_bot - d_top)
hmax <- max(hc_col$height)

# Position of a merge child: a leaf (v < 0) or an earlier merge step (v > 0)
dendro_child <- function(v, merge_x, merge_h) {
  if (v < 0) list(x = leaf_x[match(-v, hc_col$order)], h = 0)
  else list(x = merge_x[v], h = merge_h[v])
}

# Stem segments for every merge step, and the two arms of the root
dendro_geometry <- function() {
  n <- nrow(hc_col$merge)
  merge_h <- hc_col$height
  merge_x <- numeric(n)
  stems <- vector("list", n)
  for (k in seq_len(n)) {
    a <- dendro_child(hc_col$merge[k, 1], merge_x, merge_h)
    b <- dendro_child(hc_col$merge[k, 2], merge_x, merge_h)
    stems[[k]] <- list(xa = a$x, ya = height_y(a$h, hmax),
                       xb = b$x, yb = height_y(b$h, hmax),
                       yk = height_y(merge_h[k], hmax))
    merge_x[k] <- (a$x + b$x) / 2
  }
  arms <- lapply(hc_col$merge[n, ], dendro_child, merge_x = merge_x, merge_h = merge_h)
  list(stems = stems, arms = arms)
}
dendro <- dendro_geometry()

# ---- Plot canvas ------------------------------------------------------------
png(OUT_FILE, width = CANVAS_W / 100, height = CANVAS_H / 100, units = "in",
    res = DPI, bg = "white", pointsize = 12)
op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, CANVAS_W), ylim = c(CANVAS_H, 0))

# Column dendrogram
for (s in dendro$stems) {
  segments(s$xa, s$ya, s$xa, s$yk, lwd = 3.3, lend = 2)
  segments(s$xb, s$yb, s$xb, s$yk, lwd = 3.3, lend = 2)
  segments(s$xa, s$yk, s$xb, s$yk, lwd = 3.3, lend = 2)
}

# Label the two main clades under the root of the dendrogram
for (arm in 1:2) {
  side <- if (arm == 1) -1 else 1
  label <- if (arm == 1) "Clade I" else "Clade II"
  X <- dendro$arms[[arm]]$x
  Y <- height_y(dendro$arms[[arm]]$h, hmax)
  text(X + side * 8, Y - 14, label, adj = c(if (side < 0) 1 else 0, 0),
       font = 2, cex = FS_STRIP / 12)
}

# Heat cells
for (i in seq_len(nR)) for (j in seq_len(nC)) {
  rect(x0 + (j - 1) * cw, y0 + (i - 1) * rh, x0 + j * cw, y0 + i * rh,
       col = CELL_COL[i, j], border = "white", lwd = 1.9)
}

# Strips: drug class | mechanism | sharing tier
mrow <- meta[match(row_order, meta$fam), ]
strip_cols <- list(CLASS_PAL[mrow$class_lab], MECH_PAL[mrow$mech_lab], TIER_PAL[mrow$tier])
for (k in 1:3) for (i in seq_len(nR)) {
  rect(strip_x[k], y0 + (i - 1) * rh, strip_x[k] + strip_w, y0 + i * rh,
       col = strip_cols[[k]][i], border = "white", lwd = 1.3)
}
for (k in 1:3) {
  label <- c("Drug class", "Mechanism", "Sharing tier")[k]
  text(strip_x[k] + strip_w / 2, y0 - 10, label, srt = 90, adj = c(0, 0.5),
       cex = FS_STRIP / 12)
}

# Row labels
text(x0 - 14, y0 + (seq_len(nR) - 0.5) * rh, row_order, adj = c(1, 0.5), cex = FS_ROW / 12)

# Column labels, italic, 45 degrees
col_labels <- BINOM[col_order]
text(x0 + (seq_len(nC) - 0.5) * cw + 14, y1 + 16, col_labels, srt = 45, adj = c(1, 1),
     font = 3, cex = FS_COL / 12)
lowest <- max(y1 + 16 + strwidth(col_labels, cex = FS_COL / 12) * sin(pi / 4))
text((x0 + x1) / 2, lowest + 45, "Host species", adj = c(0.5, 0.5), cex = FS_AXIS / 12)

# Legends
legend_block <- function(y, title, labels, colours, step = 44) {
  text(lx, y, title, adj = c(0, 0.5), font = 2, cex = FS_LEG_T / 12)
  for (i in seq_along(labels)) {
    yy <- y + 46 + (i - 1) * step
    rect(lx, yy - 17, lx + 34, yy + 17, col = colours[i], border = "grey", lwd = 0.6)
    text(lx + 50, yy, labels[i], adj = c(0, 0.5), cex = FS_LEG / 12)
  }
  y + 46 + (length(labels) - 1) * step + 17
}

cls_used  <- sort(unique(meta$class_lab))
mech_used <- sort(unique(meta$mech_lab))
tier_used <- names(TIER_PAL)[names(TIER_PAL) %in% meta$tier]

bottom <- legend_block(y0 + 8, "Drug class", cls_used, CLASS_PAL[cls_used])
bottom <- legend_block(bottom + 50, "Resistance mechanism", mech_used, MECH_PAL[mech_used])
bottom <- legend_block(bottom + 50, "Sharing tier",
                       ifelse(tier_used %in% names(TIER_LABEL), TIER_LABEL[tier_used], tier_used),
                       TIER_PAL[tier_used])

# Colour bar
bar_w <- 340; bar_h <- 26
y_bar <- bottom + 110
bar_top <- y_bar + 70
bar_bottom <- bar_top + bar_h
text(lx, y_bar, "log10(hits per 1,000\nrepresentative genomes)", adj = c(0, 0),
     font = 2, cex = FS_LEG_T / 12)
seg_w <- bar_w / 100
for (k in 1:100) {
  rect(lx + (k - 1) * seg_w, bar_top, lx + k * seg_w + 0.5, bar_bottom,
       col = RAMP[k], border = NA)
}
rect(lx, bar_top, lx + bar_w, bar_bottom, border = "grey", lwd = 0.6)
ticks <- pretty(c(r0, r1), n = 3)
ticks <- ticks[ticks >= r0 & ticks <= r1]
for (t in ticks) {
  tx <- lx + (t - r0) / (r1 - r0) * bar_w
  segments(tx, bar_bottom, tx, bar_bottom + 7, lwd = 0.9)
  text(tx, bar_bottom + 26, sprintf("%.1f", t), cex = FS_TICK / 12)
}

# No-hit swatch
nh_y <- bar_bottom + 66
rect(lx, nh_y - 17, lx + 34, nh_y + 17, col = ZERO_COL, border = "grey", lwd = 0.6)
text(lx + 50, nh_y, "No hits", adj = c(0, 0.5), cex = FS_LEG / 12)

par(op)
invisible(dev.off())

message("Hosts shown (", nC, "): ", paste(BINOM[col_order], collapse = ", "))
message("Figure written to ", OUT_FILE, " at ", DPI, " dpi")
