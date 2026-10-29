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
#  SCRIPT  : figC_b_subpathway_phylum.R  -  Sub-pathway phylum composition and tier prevalence (Fig. 3b)
#  AUTHOR  : QRQ
#  INPUT   : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig2/supp/representative_host_range.csv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/subpathway_presence_propagated_50.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/gtdbtk.bac120.summary.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/gtdbtk.ar53.summary.tsv
#  OUTPUT  : C:/qrq_metagenome_pipeline/Plots and Figures/figC_b_subpathway_phylum.png (600 dpi)
#             C:/qrq_metagenome_pipeline/Plots and Figures/figC_b_subpathway_phylum.pdf
#             C:/qrq_metagenome_pipeline/Plots and Figures/figC_b_subpathway_phylum_composition.csv
# ==========================================================
#
# Left : sub-pathway labels, a bar per sub-pathway with the phylum composition of
#        the genomes that carry it, and a pie for all genomes as the reference.
# Right: prevalence by sharing tier (%), unadjusted (open) and within-lineage
#        (filled) odds ratios on a log axis, and the >= 75% / 100% marker tiles.

BASE       <- "C:/qrq_metagenome_pipeline/manuscript/Figures"
SHARING_F  <- file.path(BASE, "Fig2/supp/representative_host_range.csv")
PATH_F     <- file.path(BASE, "Fig3/supp/subpathway_presence_propagated_50.tsv")
GTDB_BAC_F <- file.path(BASE, "Fig3/supp/gtdbtk.bac120.summary.tsv")
GTDB_AR_F  <- file.path(BASE, "Fig3/supp/gtdbtk.ar53.summary.tsv")
OUT_DIR    <- "C:/qrq_metagenome_pipeline/Plots and Figures"

N_PHYLA <- 6
X_MAX   <- 10
TIER_COL <- c(Core = "#B06040", Shared = "#C89040", `Host-specific` = "#2898B8")
INK     <- "#1A1A1A"
NS      <- "#8C8C8C"
NS_FILL <- "#E9E6E1"
PHYLUM_COL <- c("#1F5F87", "#6AA6C9", "#3B3B98", "#4FA36B", "#8E5BA8", "#E0B97A")
OTHER_COL  <- "#D0D0D0"
DARK_PHYLA <- c("Bacillota", "Bacillota_I", "Verrucomicrobiota")

MAP <- c(butyrate_kinase_route = "Butyrate (kinase)",
         propionate_succinate = "Propionate (succinate)",
         butyrate_transferase_route = "Butyrate (CoA-transferase)",
         acetate_production = "Acetate",
         cobalamin_biosynthesis = "Cobalamin",
         urea_utilization = "Urease",
         ammonia_assimilation = "Ammonia assimilation",
         nitrate_reduction = "Nitrate reduction",
         bcaa_biosynthesis = "Branched-chain amino acids",
         lysine_biosynthesis = "Lysine")

# Forest values from the template script (panel b)
sp <- data.frame(
  label = c("Butyrate (kinase)", "Propionate (succinate)", "Butyrate (CoA-transferase)",
            "Acetate", "Cobalamin", "Urease", "Ammonia assimilation", "Nitrate reduction",
            "Branched-chain amino acids", "Lysine"),
  markers = c(3, 5, 3, 2, 4, 3, 3, 2, 2, 2),
  Core = c(32.4, 1.32, 9.25, 89.9, 17.4, 6.61, 74.2, 1.76, 57.7, 88.6),
  Shared = c(21.9, 0.69, 6.80, 91.8, 14.8, 2.60, 69.5, 2.85, 50.7, 88.7),
  `Host-specific` = c(15.8, 0.54, 5.68, 91.4, 13.9, 1.40, 66.9, 2.81, 46.4, 88.1),
  or_u = c(2.55, 2.47, 1.69, 0.83, 1.30, 4.99, 1.43, 0.62, 1.58, 1.05),
  lo_u = c(2.08, 0.88, 1.20, 0.61, 1.01, 3.27, 1.15, 0.27, 1.30, 0.78),
  hi_u = c(3.13, 5.59, 2.34, 1.16, 1.67, 7.39, 1.78, 1.24, 1.91, 1.43),
  p_u  = c(1.1e-16, 0.079, 6.1e-3, 0.35, 0.079, 7.8e-11, 1.9e-3, 0.35, 8.7e-6, 1.0),
  or_w = c(1.24, 5.65, 1.86, 0.82, 1.10, 1.26, 0.73, 0.69, 1.15, 0.39),
  lo_w = c(0.85, 0.56, 0.90, 0.46, 0.73, 0.70, 0.53, 0.25, 0.82, 0.21),
  hi_w = c(1.81, 57.2, 3.83, 1.49, 1.66, 2.28, 1.01, 1.90, 1.63, 0.69),
  p_w  = c(0.55, 0.55, 0.55, 0.73, 0.85, 0.73, 0.017, 0.55, 0.73, 0.040),
  or_75 = c(1.32, NA, 2.09, 1.41, 1.34, 5.87, 1.51, 3.71, 1.76, 1.13),
  call_75 = c("ns", "nt", "core", "core", "ns", "core", "core", "ns", "core", "ns"),
  stringsAsFactors = FALSE, check.names = FALSE
)
sp$or_100   <- ifelse(sp$markers <= 3, sp$or_75, NA)
sp$call_100 <- ifelse(sp$markers <= 3, sp$call_75, "nt")

# Significant = BH p < 0.05 and the CI excludes 1
call_rule <- function(lo, hi, p) p < 0.05 & (lo > 1 | hi < 1)
sp$sig_u <- call_rule(sp$lo_u, sp$hi_u, sp$p_u)
sp$sig_w <- call_rule(sp$lo_w, sp$hi_w, sp$p_w)
stopifnot(sum(sp$sig_u) == 5, identical(sp$label[sp$sig_w], "Lysine"))

# Row positions: groups are separated by gaps (SCFA | cofactor | nitrogen | amino acids)
GAP <- c(rep(1.5, 4), 1.0, rep(0.5, 3), 0, 0)
sp$y <- (nrow(sp):1) + GAP
TOP_Y <- max(sp$y)
REF_Y <- TOP_Y + 2.0

# ---- Phylum of every representative ---------------------------------------
read_gtdb <- function(path) {
  x <- read.delim(path, colClasses = "character", stringsAsFactors = FALSE)
  x[, c("user_genome", "classification")]
}

sh  <- read.csv(SHARING_F, colClasses = "character", stringsAsFactors = FALSE)["representative"]
tax <- rbind(read_gtdb(GTDB_BAC_F), read_gtdb(GTDB_AR_F))
phylum_of <- function(cls) {
  ph <- ifelse(grepl("p__", cls), sub(".*p__([^;]*).*", "\\1", cls), "")
  ifelse(ph == "", "Unclassified", ph)
}
tax$phylum <- phylum_of(ifelse(is.na(tax$classification), "", tax$classification))

g <- merge(sh, tax[, c("user_genome", "phylum")], by.x = "representative",
           by.y = "user_genome", all.x = TRUE, sort = FALSE)
g$phylum[is.na(g$phylum)] <- "Unclassified"

counts <- sort(table(g$phylum[g$phylum != "Unclassified"]), decreasing = TRUE)
top <- names(counts)[seq_len(min(N_PHYLA, length(counts)))]
g$ph <- ifelse(g$phylum %in% top, g$phylum, "Other phyla")
LV  <- c(top, "Other phyla")
PHC <- setNames(c(PHYLUM_COL[seq_along(top)], OTHER_COL), LV)
message("Phyla shown: ", paste(top, collapse = ", "))

# ---- Carriers of each sub-pathway, by phylum -------------------------------
pr <- read.delim(PATH_F, colClasses = "character", stringsAsFactors = FALSE)[, 1:2]
names(pr) <- c("genome", "subpathway")
pr <- unique(pr[pr$subpathway %in% names(MAP), ])
pr$label <- MAP[pr$subpathway]
m <- merge(pr, g[, c("representative", "ph")], by.x = "genome", by.y = "representative",
           all.x = TRUE, sort = FALSE)
matched <- mean(!is.na(m$ph))
message(sprintf("Carrier rows: %d | matched to a phylum: %.1f%%", nrow(m), 100 * matched))
if (matched < 0.95) stop("Under 95% of carrier genomes matched the representatives: check genome ID format.")
m <- m[!is.na(m$ph), ]

comp <- as.data.frame(table(label = m$label, ph = m$ph), stringsAsFactors = FALSE)
comp <- comp[comp$Freq > 0, ]
names(comp)[3] <- "N"
comp$share <- 100 * comp$N / ave(comp$N, comp$label, FUN = sum)
ref <- as.data.frame(table(ph = g$ph), stringsAsFactors = FALSE)
names(ref)[2] <- "N"
ref$label <- "All genomes"
ref$share <- 100 * ref$N / sum(ref$N)
allc <- rbind(comp[, c("label", "ph", "N", "share")], ref[, c("label", "ph", "N", "share")])
write.csv(allc, file.path(OUT_DIR, "figC_b_subpathway_phylum_composition.csv"), row.names = FALSE)

shares_of <- function(label) {
  d <- allc$share[allc$label == label][match(LV, allc$ph[allc$label == label])]
  d[is.na(d)] <- 0
  d
}

# Break a long name onto two lines at the space that balances the halves
break_name <- function(s, limit = 14) {
  if (nchar(s) <= limit || !grepl(" ", s)) return(s)
  words <- strsplit(s, " ")[[1]]
  best <- NULL
  for (i in seq_len(length(words) - 1)) {
    a <- paste(words[1:i], collapse = " ")
    b <- paste(words[(i + 1):length(words)], collapse = " ")
    w <- max(nchar(a), nchar(b))
    if (is.null(best) || w < best$w) best <- list(w = w, s = paste0(a, "\n", b))
  }
  best$s
}

# Polygon points of a circular sector (angles in degrees, counter-clockwise from +x)
sector <- function(cx, cy, r, a1, a2, n = 60) {
  th <- seq(a1, a2, length.out = n) * pi / 180
  list(x = c(cx, cx + r * cos(th)), y = c(cy, cy + r * sin(th)))
}

# ---- Figure ---------------------------------------------------------------
FIG_W <- 7.2
FIG_H <- 7.0
AX_BOTTOM <- 0.17
AX_H <- 0.70
PT <- function(pt) pt / 12        # font size in points, for png pointsize 12
FS <- 10; FS_SM <- 8.5; FS_TILE <- 8

Y_LO <- min(sp$y) - 0.7
Y_HI <- 17.0

draw_figure <- function() {

# ---- Left panel: phylum composition ----
unit_in <- AX_H * FIG_H / (Y_HI - Y_LO)
W_UNITS <- 0.46 * FIG_W / unit_in
BAR_W   <- 1.45 / unit_in
BAR_X0  <- W_UNITS - BAR_W - 0.1
BAR_H   <- 0.62
REF_R   <- 1.1

par(fig = c(0.01, 0.47, AX_BOTTOM, AX_BOTTOM + AX_H), mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, W_UNITS), ylim = c(Y_LO, Y_HI), xaxs = "i", yaxs = "i")

rows_l <- data.frame(label = c(sp$label, "All genomes"), y = c(sp$y, REF_Y),
                     is_ref = c(rep(FALSE, nrow(sp)), TRUE), stringsAsFactors = FALSE)
for (k in seq_len(nrow(rows_l))) {
  lab <- rows_l$label[k]; y <- rows_l$y[k]; is_ref <- rows_l$is_ref[k]
  txt <- break_name(if (lab == "Lysine") paste(lab, "\u2020") else lab)
  text(BAR_X0 - 0.35, y, txt, adj = c(1, 0.5), cex = PT(FS),
       col = if (is_ref) "#555555" else INK)

  vals <- shares_of(lab)
  if (is_ref) {
    # the all-genome reference is a pie
    cx <- BAR_X0 + BAR_W / 2
    start <- 90
    for (i in seq_along(vals)) {
      if (vals[i] <= 0) next
      sweep <- 360 * vals[i] / 100
      s <- sector(cx, y, REF_R, start - sweep, start)
      polygon(s$x, s$y, col = PHC[LV[i]], border = "white", lwd = 0.5)
      if (vals[i] >= 8) {
        mid <- (start - sweep / 2) * pi / 180
        text(cx + 0.6 * REF_R * cos(mid), y + 0.6 * REF_R * sin(mid),
             sprintf("%.0f", vals[i]), cex = PT(FS_TILE),
             col = if (LV[i] %in% DARK_PHYLA) "white" else INK)
      }
      start <- start - sweep
    }
  } else {
    # sub-pathways are 100% bars
    x_cur <- BAR_X0
    for (i in seq_along(vals)) {
      if (vals[i] <= 0) next
      w <- BAR_W * vals[i] / 100
      rect(x_cur, y - BAR_H / 2, x_cur + w, y + BAR_H / 2, col = PHC[LV[i]],
           border = "white", lwd = 0.5)
      if (vals[i] >= 8) {
        text(x_cur + w / 2, y, sprintf("%.0f", vals[i]), cex = PT(FS_TILE),
             col = if (LV[i] %in% DARK_PHYLA) "white" else INK)
      }
      x_cur <- x_cur + w
    }
  }
}

# ---- Right panel: forest (log10 odds-ratio scale) ----
X_LAB <- c(Core = -2.55, Shared = -1.92, `Host-specific` = -1.30)   # log10 positions
TILE_X <- c(`75` = 1.36, `100` = 1.90)                                # log10 positions
par(fig = c(0.52, 0.98, AX_BOTTOM, AX_BOTTOM + AX_H), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
plot.window(xlim = c(-2.66, 2.12), ylim = c(Y_LO, Y_HI), xaxs = "i", yaxs = "i")

abline(v = 0, col = "#8A8A8A", lty = 2, lwd = 0.8)
abline(v = log10(c(0.2, 0.5, 2, 5)), col = "#E6E6E6", lwd = 0.6)
axis(1, at = log10(c(0.2, 0.5, 1, 2, 5, 10)), labels = c("0.2", "0.5", "1", "2", "5", "10"),
     cex.axis = PT(FS - 1), col.axis = INK, lwd = 0.6, tcl = -0.3)
mtext("Odds ratio, core vs host-specific (log scale)", side = 1, line = 2.6,
      cex = PT(FS - 1), col = INK)

for (tier in names(X_LAB)) {
  text(X_LAB[tier], Y_HI - 0.6, tier, adj = c(0, 0), srt = 45, font = 2,
       cex = PT(FS - 2), col = TIER_COL[tier])
}
for (r in seq_len(nrow(sp))) {
  for (tier in names(X_LAB)) {
    text(X_LAB[tier], sp$y[r], sprintf("%.1f", sp[[tier]][r]), cex = PT(FS_SM), col = INK)
  }
}

# Odds ratios: open circle = unadjusted, filled = within lineages
for (r in seq_len(nrow(sp))) {
  for (model in c("u", "w")) {
    filled <- model == "w"
    y_off <- if (model == "u") 0.17 else -0.17
    orv <- sp[[paste0("or_", model)]][r]
    lo  <- sp[[paste0("lo_", model)]][r]
    hi  <- sp[[paste0("hi_", model)]][r]
    sig <- sp[[paste0("sig_", model)]][r]
    col <- if (sig) INK else NS
    y   <- sp$y[r] + y_off
    hi_d <- min(hi, X_MAX)
    segments(log10(lo), y, log10(hi_d), y, col = col, lwd = 1.2)
    if (hi > X_MAX) {
      arrows(log10(X_MAX / 1.25), y, log10(X_MAX), y, col = col, length = 0.05, angle = 25, lwd = 0.8)
      text(log10(X_MAX), y + (if (model == "u") 0.33 else -0.33), sprintf("to %g", hi),
           adj = c(1, 0.5), cex = PT(FS_TILE - 0.5), col = NS)
    }
    points(log10(orv), y, pch = if (filled) 16 else 1, cex = if (filled) 1.0 else 0.95,
           col = col, bg = if (filled) col else "white", lwd = 0.9)
  }
}

# Threshold tiles (>= 75% and 100% markers)
for (r in seq_len(nrow(sp))) {
  for (key in c("75", "100")) {
    call <- sp[[paste0("call_", key)]][r]
    or_v <- sp[[paste0("or_", key)]][r]
    cx <- TILE_X[key]
    if (call == "core") { fc <- TIER_COL["Core"]; tc <- "white" }
    else if (call == "ns") { fc <- NS_FILL; tc <- INK }
    else { fc <- "white"; tc <- NS }
    rect(cx - 0.15, sp$y[r] - 0.41, cx + 0.15, sp$y[r] + 0.41, col = fc, border = "#C9C9C9", lwd = 0.5)
    txt <- if (is.na(or_v)) "n.t." else sprintf("%.2f", or_v)
    text(cx, sp$y[r], txt, cex = PT(FS_TILE), col = tc)
  }
}
for (key in names(TILE_X)) {
  text(TILE_X[key], Y_HI - 0.6, if (key == "75") "\u226575%" else "100%", font = 2,
       cex = PT(FS - 2), col = INK)
}

# ---- Full-figure overlay: titles and legend ----
par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.01, 0.955, "Phylum composition of carrying genomes", adj = c(0, 0.5), cex = PT(FS - 0.5), col = INK)
text(0.52, 0.955, "Prevalence (%) and odds ratios", adj = c(0, 0.5), cex = PT(FS - 0.5), col = INK)

legend_fill  <- c(PHC[LV], NA, NA, NA, NA)
legend_pch   <- c(rep(NA, length(LV)), 1, 16, NA, NA)
legend_lty   <- c(rep(NA, length(LV)), NA, NA, 1, 1)
legend_col   <- c(rep(NA, length(LV)), INK, INK, INK, NS)
legend_text  <- c(LV, "Unadjusted", "Within lineages", "BH p < 0.05, CI excludes 1", "Not significant")
legend("bottom", inset = 0.005, legend = legend_text, fill = legend_fill, border = NA,
       pch = legend_pch, lty = legend_lty, col = legend_col, lwd = 1,
       ncol = 4, bty = "n", cex = PT(9.5), x.intersp = 0.6)
}

png(file.path(OUT_DIR, "figC_b_subpathway_phylum.png"), width = FIG_W, height = FIG_H,
    units = "in", res = 600, bg = "white", pointsize = 12)
draw_figure()
invisible(dev.off())

pdf(file.path(OUT_DIR, "figC_b_subpathway_phylum.pdf"), width = FIG_W, height = FIG_H,
    pointsize = 12)
draw_figure()
invisible(dev.off())
