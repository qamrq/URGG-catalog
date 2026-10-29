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
#  SCRIPT  : figA_thiamine_guilds.R  -  Thiamine producer and salvager guilds (Fig. A)
#  AUTHOR  : QRQ
#  INPUT   : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig4/supp/63_genome_thiamine_guilds.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig4/supp/table_s49_thiamine_guild_presence_models.csv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig2/supp/representative_host_range.csv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/gtdbtk.bac120.summary.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/gtdbtk.ar53.summary.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig5/supp/lineage_full_classification_table.csv
#  OUTPUT  : C:/qrq_metagenome_pipeline/Plots and Figures/Fig_thiamine_guilds.png (600 dpi)
#             C:/qrq_metagenome_pipeline/Plots and Figures/Fig_thiamine_guilds.pdf
#             C:/qrq_metagenome_pipeline/Plots and Figures/Fig_thiamine_guilds.svg
#             C:/qrq_metagenome_pipeline/Plots and Figures/panel_c_plot_table.csv
#             C:/qrq_metagenome_pipeline/Plots and Figures/figA_genus_thiamine_summary.csv
# ==========================================================
#
# a  pathway schematic: de novo route (KEGG M00127, ThiH type) | membrane ring around
#    the thiamine pool | salvage routes, with representative genera
# b  guild composition by host-sharing tier
# c  raw (Fisher) vs within-lineage (GLMM) odds ratios, core/broad vs host-specific
#
# The whole figure is drawn on one canvas in millimetres (183 x 170 mm).

ROOT       <- "C:/qrq_metagenome_pipeline/manuscript/Figures"
G63_F      <- file.path(ROOT, "Fig4/supp/63_genome_thiamine_guilds.tsv")
GLMM_F     <- file.path(ROOT, "Fig4/supp/table_s49_thiamine_guild_presence_models.csv")
SHARING_F  <- file.path(ROOT, "Fig2/supp/representative_host_range.csv")
GTDB_BAC_F <- file.path(ROOT, "Fig3/supp/gtdbtk.bac120.summary.tsv")
GTDB_AR_F  <- file.path(ROOT, "Fig3/supp/gtdbtk.ar53.summary.tsv")
LINEAGE_F  <- file.path(ROOT, "Fig5/supp/lineage_full_classification_table.csv")
OUT_DIR    <- "C:/qrq_metagenome_pipeline/Plots and Figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXCLUDE <- character(0)   # genera left out of panel c

TIER_HOST   <- "Host-specific (1 host species)"
TIER_SHARED <- "Shared (2-3 host species)"
TIER_CORE   <- "Core / broad (>3 host species)"
TIERS       <- c(TIER_CORE, TIER_SHARED, TIER_HOST)
TIER_LAB    <- c(setNames("Core/broad", TIER_CORE), setNames("Shared", TIER_SHARED),
                 setNames("Host-specific", TIER_HOST))

PRODUCER <- list(dark = "#A9423F", mid = "#CF6A5F", light = "#EDA99E", wash = "#FCF1EE", line = "#E2A097")
SALVAGER <- list(dark = "#1B6B70", mid = "#3A9CA0", light = "#93CDC9", wash = "#EBF6F5", line = "#86C4C0")
NEUTRAL  <- list(ink = "#222222", grey = "#7F7D76", light = "#BDBAB1", wash = "#F4F3EF", pool = "#3D4652")
SHADOW   <- "#3A3A3A"
GUILD_COL <- c("1_de_novo" = PRODUCER$dark, "2_partial_de_novo" = PRODUCER$light,
               "3_precursor_salvage" = SALVAGER$light, "4_thiamine_salvage" = SALVAGER$dark,
               "5_none_detected" = NEUTRAL$light)
GUILD_LAB <- c("1_de_novo" = "Complete de novo", "2_partial_de_novo" = "Partial de novo (thiC+)",
               "3_precursor_salvage" = "Precursor salvage", "4_thiamine_salvage" = "Thiamine salvage",
               "5_none_detected" = "No route detected")
GUILDS <- names(GUILD_COL)

FIG_W <- 183
FIG_H <- 170
FIG_BOTTOM <- 22                                # lowest drawn height (mm); canvas is cut here
FIG_OUT_H <- FIG_H - FIG_BOTTOM                 # output height (mm)
PT_PER_MM <- 72.27 / 25.4                       # mm -> points
pt  <- function(p) p / 12                       # points -> cex (png and pdf use pointsize 12)
lwd_of <- function(lw) max(lw / 0.75, 0.01)     # points -> lwd, never zero
mm  <- function(v) pt(v * PT_PER_MM)           # millimetres of text height -> cex

host_count <- function(s) {
  s <- s[!is.na(s)]
  if (!length(s)) return(NA_integer_)
  names <- trimws(unlist(strsplit(as.character(s), "[;,|]")))
  length(unique(names[nzchar(names)]))
}
fmt_int <- function(x) formatC(round(x), format = "d", big.mark = ",")

# ---- Data ------------------------------------------------------------------
for (f in c(G63_F, GLMM_F, SHARING_F)) if (!file.exists(f)) stop("Missing input: ", f)

g <- read.delim(G63_F, stringsAsFactors = FALSE)
sh <- read.csv(SHARING_F, colClasses = "character", stringsAsFactors = FALSE)
sh <- sh[, intersect(c("representative", "host_species_list"), names(sh)), drop = FALSE]
names(sh)[names(sh) == "representative"] <- "genome"
g <- merge(g, sh, by = "genome", all.x = TRUE, sort = FALSE)
if (!"host_species_list" %in% names(g)) g$host_species_list <- NA_character_

tax_parts <- lapply(Filter(file.exists, c(GTDB_BAC_F, GTDB_AR_F)), function(f) {
  x <- read.delim(f, colClasses = "character", quote = "", stringsAsFactors = FALSE)
  data.frame(genome = x$user_genome, classification = x$classification, stringsAsFactors = FALSE)
})
if (length(tax_parts)) {
  tax <- do.call(rbind, tax_parts)
  genus <- sub(".*g__([^;]*).*", "\\1", ifelse(is.na(tax$classification), "", tax$classification))
  genus[!grepl("g__", ifelse(is.na(tax$classification), "", tax$classification))] <- ""
  tax$genus  <- ifelse(nzchar(genus), genus, "(no genus)")
  tax$domain <- ifelse(grepl("d__Archaea", tax$classification), "Archaea", "Bacteria")
  g <- merge(g, tax[, c("genome", "genus", "domain")], by = "genome", all.x = TRUE, sort = FALSE)
}
if (!"genus" %in% names(g)) g$genus <- "(no genus)"
if (!"domain" %in% names(g)) g$domain <- "Bacteria"
g$genus[is.na(g$genus)] <- "(no genus)"
g$domain[is.na(g$domain)] <- "Bacteria"

N_ALL <- nrow(g)
message("Genomes: ", fmt_int(N_ALL))
gpct <- table(g$guild)
gp <- function(x) if (x %in% names(gpct)) 100 * gpct[[x]] / N_ALL else 0
pct_thiC <- 100 * mean(g$hmp_denovo)

# Representative genera per guild
g$is_complete <- g$guild == "1_de_novo"
g$is_salv     <- g$guild == "4_thiamine_salvage"
g$is_prec     <- g$guild == "3_precursor_salvage"
gg <- g[g$genus != "(no genus)", ]
gen_parts <- lapply(split(gg, list(gg$genus, gg$domain), drop = TRUE), function(s) {
  data.frame(genus = s$genus[1], domain = s$domain[1], n = nrow(s),
             n_thiC = sum(s$hmp_denovo), n_complete = sum(s$is_complete),
             n_salv = sum(s$is_salv), n_prec = sum(s$is_prec), n_tr = sum(s$transport),
             n_hosts = host_count(s$host_species_list), stringsAsFactors = FALSE)
})
gen <- do.call(rbind, gen_parts)
gen$f_thiC     <- gen$n_thiC / gen$n
gen$f_complete <- gen$n_complete / gen$n
gen$f_salv     <- gen$n_salv / gen$n
gen$f_tr       <- gen$n_tr / gen$n
gen$named      <- grepl("^[A-Z][a-z]+(_[A-Z])?$", gen$genus)

pick_distinct <- function(d, k) {
  d <- d[!(d$genus %in% EXCLUDE), ]
  base <- sub("_[A-Z]+$", "", d$genus)     # skip sister placeholders (Phocaeicola_A after Phocaeicola)
  d <- d[!duplicated(base), ]
  d[seq_len(min(k, nrow(d))), ]
}
s_prod <- gen[gen$n >= 20 & gen$f_thiC >= 0.5 & gen$domain == "Bacteria", ]
s_prod <- s_prod[order(-s_prod$n_complete, -s_prod$n_thiC), ]
rep_prod <- pick_distinct(s_prod, 4)
s_salv <- gen[gen$n >= 20 & gen$f_salv >= 0.6 & gen$f_thiC < 0.05 & gen$f_tr >= 0.5, ]
s_salv <- s_salv[order(-s_salv$n_salv), ]
rep_salv <- pick_distinct(s_salv, 3)
write.csv(gen[order(-gen$n_thiC), ], file.path(OUT_DIR, "figA_genus_thiamine_summary.csv"),
          row.names = FALSE)
message("Representative producer genera: ", paste(rep_prod$genus, collapse = ", "))
message("Representative salvager genera: ", paste(rep_salv$genus, collapse = ", "))

# RUC lineage
ruc <- NULL
if (file.exists(LINEAGE_F)) {
  lin <- read.csv(LINEAGE_F, stringsAsFactors = FALSE)$genome
  r <- g[g$genome %in% lin, ]
  if (nrow(r)) {
    ruc <- list(n = nrow(r), tr = 100 * mean(r$transport), tpk = 100 * mean(r$tpk_thiN),
                thiL = 100 * mean(r$tmp_kinase_thiL), salv = 100 * mean(r$is_salv),
                n_hosts = host_count(r$host_species_list))
  }
}

# GLMM results (panel c)
trait_lab <- c(hmp_denovo = "thiC (de novo, committed step)",
               tpk_thiN = "TPK (activates imported thiamine)",
               hmp_kinase = "HMP kinase (ThiD/DE/DN)",
               denovo_module_ge_thresh = "De novo module (at least 50% complete)",
               tmp_kinase_thiL = "ThiL (TMP kinase)",
               transport = "Thiamine transporter",
               `guild:1_de_novo` = "Complete de novo guild",
               `guild:2_partial_de_novo` = "Partial de novo guild (thiC+)",
               `guild:4_thiamine_salvage` = "Thiamine-salvage guild",
               `guild:5_none_detected` = "No route detected")
gm <- read.csv(GLMM_F, stringsAsFactors = FALSE)
full <- gm[gm$trait %in% names(trait_lab) &
           gm$adjustment == "adjusted for log genome size + MAG completeness", ]
full <- full[match(names(trait_lab), full$trait), ]
full <- full[!is.na(full$trait), ]
full$label <- unname(trait_lab[full$trait])
# BH-adjusted LRT p < 0.05 and the within-lineage CI excludes 1.
# Non-converged fits have LRT p = NA: drawn open with a dagger, never filled.
full$nonconv <- grepl("check_convergence", full$status)
full$sig <- !full$nonconv & !is.na(full$LRT_p_BH) & full$LRT_p_BH < 0.05 &
            (full$OR_glmm_low > 1 | full$OR_glmm_high < 1)
full$plot_label <- paste0(full$label, ifelse(full$nonconv, " \u2020", ""))
write.csv(data.frame(trait = full$trait, label = full$label, OR_fisher = full$OR_fisher,
                     OR_glmm = full$OR_glmm, CI_low = full$OR_glmm_low, CI_high = full$OR_glmm_high,
                     LRT_p_BH = full$LRT_p_BH, status = full$status, nonconv = full$nonconv,
                     filled = full$sig),
          file.path(OUT_DIR, "panel_c_plot_table.csv"), row.names = FALSE)

# ---- Coordinate maps ----------------------------------------------------------
# Schematic: data units 200 x 114, drawn at 0.75 mm per unit, centred in its box
SCH_K  <- 0.75
SCH_X  <- function(x) 14 + SCH_K * x
SCH_Y  <- function(y) 79 + SCH_K * (y + 9)
SZ     <- function(v) SCH_K * v

# Panels b and c sit below the schematic; BC_DY moves both up together
BC_DY <- 18

# Panel b: x 0-100 % in 32-92 mm; y tiers in 34-74 mm (plus BC_DY)
PB_X <- function(v) 32 + v * (60 / 100)
PB_Y <- function(y) 34 + BC_DY + (y + 0.6) * (40 / 3.3)

# Panel c: log10 odds ratio 0.2-7 in 139-178 mm
PC_LO <- log10(0.2); PC_HI <- log10(7)
PC_X  <- function(v) 139 + (log10(v) - PC_LO) * (39 / (PC_HI - PC_LO))

# Genus dot plots: one per guild panel (x in mm of the schematic box)
S_X <- 158 / 200
S_Y <- 90 / 114
LABEL_MM <- 14
PLOT_MM  <- 34

# ---- Drawing helpers (millimetre canvas) --------------------------------------
rounded_poly <- function(x, y, w, h, r, n = 6) {
  r <- min(r, w / 2, h / 2)
  corner <- function(cx, cy, a0) {
    th <- seq(a0, a0 + pi / 2, length.out = n)
    cbind(cx + r * cos(th), cy + r * sin(th))
  }
  pts <- rbind(corner(x + w - r, y + r, -pi / 2), corner(x + w - r, y + h - r, 0),
               corner(x + r, y + h - r, pi / 2), corner(x + r, y + r, pi))
  list(x = pts[, 1], y = pts[, 2])
}
rr_mm <- function(x, y, w, h, r, fc, ec = NA, lw = 0.5, alpha = 1) {
  p <- rounded_poly(x, y, w, h, r)
  polygon(p$x, p$y, col = adjustcolor(fc, alpha), border = ec, lwd = lwd_of(lw))
}
# Schematic shapes in data units
rr <- function(x, y, w, h, r, fc, ec = NA, lw = 0.5, alpha = 1) {
  rr_mm(SCH_X(x), SCH_Y(y), SZ(w), SZ(h), SZ(r), fc, ec, lw, alpha)
}
circle_d <- function(cx, cy, r, fc, ec = NA, lw = 0.5, alpha = 1) {
  th <- seq(0, 2 * pi, length.out = 90)
  polygon(SCH_X(cx) + SZ(r) * cos(th), SCH_Y(cy) + SZ(r) * sin(th),
          col = adjustcolor(fc, alpha), border = ec, lwd = lwd_of(lw))
}
ellipse_d <- function(cx, cy, w, h, angle, fc, ec = NA, lw = 0.5) {
  th <- seq(0, 2 * pi, length.out = 90)
  a <- angle * pi / 180
  ex <- w / 2 * cos(th); ey <- h / 2 * sin(th)
  px <- ex * cos(a) - ey * sin(a); py <- ex * sin(a) + ey * cos(a)
  polygon(SCH_X(cx) + SZ(px), SCH_Y(cy) + SZ(py), col = fc, border = ec, lwd = lwd_of(lw))
}
arc_d <- function(cx, cy, d, t1, t2, col, lw) {
  th <- seq(t1, t2, length.out = 60) * pi / 180
  lines(SCH_X(cx) + SZ(d / 2) * cos(th), SCH_Y(cy) + SZ(d / 2) * sin(th),
        col = col, lwd = lwd_of(lw))
}
arrow_d <- function(x1, y1, x2, y2, col, lw = 1.3) {
  arrows(SCH_X(x1), SCH_Y(y1), SCH_X(x2), SCH_Y(y2), col = col, length = 0.045,
         angle = 22, lwd = lwd_of(lw))
}
txt_d <- function(x, y, s, mm_size, col = NEUTRAL$ink, ...) {
  text(SCH_X(x), SCH_Y(y), s, cex = mm(mm_size), col = col, ...)
}

node <- function(x, y, label, side, w = 15, h = 6.4, size = 2.5) {
  fc <- if (side == "p") PRODUCER$wash else SALVAGER$wash
  ec <- if (side == "p") PRODUCER$line else SALVAGER$line
  rr(x - w / 2 + 0.7, y - h / 2 - 0.8, w, h, 1.4, SHADOW, alpha = 0.13)
  rr(x - w / 2, y - h / 2, w, h, 1.4, fc, ec, 1.0)
  txt_d(x, y, label, size)
}

label_h <- function(x, y, enzyme, ko, col) {
  txt_d(x, y + 2.2, enzyme, 2.25, col = col, font = 4)
  if (nzchar(ko)) txt_d(x, y - 2.3, ko, 1.85, col = NEUTRAL$grey)
}
label_v <- function(x, y, enzyme, ko, col, adj = 0) {
  txt_d(x, y, enzyme, 2.25, col = col, adj = c(adj, 0.5))
  if (nzchar(ko)) {
    n_e <- lengths(regmatches(enzyme, gregexpr("\n", enzyme))) + 1
    n_k <- lengths(regmatches(ko, gregexpr("\n", ko))) + 1
    half_e <- (n_e - 1) * 2.3 + 1.2            # half height of the enzyme block
    txt_d(x, y - half_e - 1.6 - (n_k - 1) * 1.0 - 0.9, ko, 1.85, col = NEUTRAL$grey,
          adj = c(adj, 0.5))
  }
}

# ---- Canvas ---------------------------------------------------------------------
draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, FIG_W), ylim = c(FIG_BOTTOM, FIG_H))

  # ---- a: schematic ----
  rr(2, 22, 88, 87.5, 3.0, PRODUCER$wash, PRODUCER$line, 0.7, alpha = 0.8)
  rr(112, 22, 86, 87.5, 3.0, SALVAGER$wash, SALVAGER$line, 0.7, alpha = 0.8)

  PD <- PRODUCER$dark; SD <- SALVAGER$dark
  txt_d(5, 104.6, "De novo producer guild", 3.3, col = PD, font = 2, adj = c(0, 0.5))
  txt_d(5, 100.4, sprintf("Complete route %.1f%%  |  thiC+ (complete or partial) %.1f%% of %s genomes",
                          gp("1_de_novo"), pct_thiC, fmt_int(N_ALL)), 1.9, adj = c(0, 0.5))
  node(11, 93, "AIR", "p", w = 13)
  node(44, 93, "HMP-P", "p")
  node(70, 93, "HMP-PP", "p")
  node(11, 77, "DXP +\ntyrosine", "p", w = 14, h = 9.6, size = 2.3)
  node(44, 77, "HET-P", "p")
  node(70, 69, "TMP", "p")
  node(70, 54, "TPP", "p")
  arrow_d(17.5, 93, 36.5, 93, PD); label_h(27, 93, "ThiC", "K03147", PD)
  arrow_d(51.5, 93, 62.5, 93, PD); label_h(57, 93, "ThiD", "K00941", PD)
  arrow_d(17.5, 77, 36.5, 77, PD); label_h(27, 77, "ThiH + ThiG", "K03150, K03149", PD)
  arrow_d(70, 89.9, 70, 72.1, PD)
  label_v(70.6, 84.6, "ThiE /\nThiDE /\nThiDN", "K00788, K14153,\nK21219", PD)
  arrow_d(51.5, 77, 62.5, 70.6, PD)
  arrow_d(70, 65.9, 70, 57.1, PD)
  label_v(72.3, 64.4, "ThiL", "K00946", PD)
  txt_d(30, 85.2, "de novo synthesis (KEGG M00127, ThiH type)", 1.95, col = PRODUCER$mid, font = 3)

  txt_d(115, 104.6, "Salvager guilds", 3.3, col = SD, font = 2, adj = c(0, 0.5))
  txt_d(115, 100.4, sprintf("Thiamine salvage %.1f%%  |  precursor salvage %.1f%%  |  no route %.1f%%",
                            gp("4_thiamine_salvage"), gp("3_precursor_salvage"), gp("5_none_detected")),
        1.8, adj = c(0, 0.5))
  node(122, 93, "Thiamine\n(external)", "s", w = 17.5, h = 8.4)
  node(160, 93, "Thiamine", "s")
  node(190, 93, "TPP", "s", w = 12)
  node(124, 77, "HMP", "s", w = 13)
  node(124, 63, "HET", "s", w = 13)
  node(159, 77, "HMP-PP", "s")
  node(159, 63, "HET-P", "s")
  node(189, 70, "TMP", "s", w = 13)

  ellipse_d(141.3, 93, 7.2, 4.4, 0, SALVAGER$light, SALVAGER$dark, 0.9)
  ellipse_d(141.3, 93, 3.4, 1.9, 0, SALVAGER$dark, NA, 0)
  arrow_d(130.75, 93, 151.5, 93, SD)
  txt_d(141.3, 97.6, "ThiT / ThiBPQ", 2.25, col = SD, font = 4)
  txt_d(141.3, 89.4, "K16789; K02062-4", 1.85, col = NEUTRAL$grey)
  arrow_d(167.5, 93, 182, 93, SD); label_h(174.8, 93, "TPK", "K00949", SD)
  arrow_d(130.5, 77, 151.5, 77, SD); label_h(141, 77, "ThiD", "K00941", SD)
  arrow_d(130.5, 63, 151.5, 63, SD); label_h(141, 63, "ThiM", "K00878", SD)
  arrow_d(166.5, 77, 182.5, 71.6, SD)
  arrow_d(166.5, 63, 182.5, 68.4, SD)
  txt_d(172.5, 70.6, "ThiE", 2.25, col = SD, font = 4)
  txt_d(172.5, 67.6, "K00788", 1.85, col = NEUTRAL$grey)
  arrow_d(189, 73.1, 189, 89.6, SD)
  label_v(187.2, 84.5, "ThiL", "K00946", SD, adj = 1)
  txt_d(141, 84.4, "thiamine salvage", 1.95, col = SALVAGER$mid, font = 3)
  txt_d(141, 56.4, "precursor salvage (KEGG M00899)", 1.95, col = SALVAGER$mid, font = 3)

  # Thiamine pool with membrane ring: producer side (left), salvage side (right)
  circle_d(102.6, 77.6, 9.2, SHADOW, alpha = 0.13)
  circle_d(102, 79, 9.2, NEUTRAL$pool)
  arc_d(102, 79, 21.4, 90, 270, PRODUCER$mid, 2.4)
  arc_d(102, 79, 21.4, -90, 90, SALVAGER$mid, 2.4)
  circle_d(102, 79, 9.9, NA, "white", 0.5, 0)
  txt_d(102, 79, "Thiamine\npool", 2.4, col = "white", font = 2)
  ellipse_d(107.0, 88.4, 4.0, 2.4, -28, SALVAGER$light, SALVAGER$dark, 0.9)
  arrow_d(77.5, 54, 94.5, 72.5, NEUTRAL$pool)
  txt_d(91.5, 60, "release", 1.9, col = NEUTRAL$grey, font = 3)
  arrow_d(108.0, 86.0, 113.5, 89.6, NEUTRAL$pool)
  arrow_d(111.2, 77.6, 117, 77, NEUTRAL$pool)
  arrow_d(109, 73.2, 117, 65, NEUTRAL$pool)

  # Representative genera as dot plots
  txt_d(5, 40, "Representative producer genera (this catalogue)", 2.35, col = PD, font = 2,
        adj = c(0, 0.5))
  if (nrow(rep_prod)) genus_plot(rep_prod, "p", 7, 35)
  txt_d(115, 40, "Representative salvager genera (this catalogue)", 2.35, col = SD, font = 2,
        adj = c(0, 0.5))
  salv_tbl <- head(rep_salv, if (is.null(ruc)) 4 else 3)
  if (!is.null(ruc)) {
    ruc_row <- data.frame(genus = "RUC lineage", n = ruc$n, n_hosts = ruc$n_hosts,
                          f_salv = ruc$salv / 100, f_tr = ruc$tr / 100, named = FALSE,
                          f_complete = NA, f_thiC = NA, stringsAsFactors = FALSE)
    salv_tbl <- rbind(salv_tbl[, names(ruc_row)], ruc_row)
  }
  if (nrow(salv_tbl)) genus_plot(salv_tbl, "s", 117, 35)

  # ---- b: guild composition by tier ----
  counts <- table(factor(g$sharing_category, levels = TIERS), factor(g$guild, levels = GUILDS))
  pct <- counts / rowSums(counts) * 100
  ntier <- rowSums(counts)
  ypos <- setNames(2 - (seq_along(TIERS) - 1), TIERS)
  for (t in TIERS) {
    left <- 0
    for (guild in GUILDS) {
      w <- pct[t, guild]
      if (w <= 0) next
      rect(PB_X(left), PB_Y(ypos[t]) - 0.66 / 2 * (40 / 3.3),
           PB_X(left + w), PB_Y(ypos[t]) + 0.66 / 2 * (40 / 3.3),
           col = GUILD_COL[guild], border = "white", lwd = 0.4 / 0.75)
      txt_col <- if (guild %in% c("1_de_novo", "4_thiamine_salvage")) "white" else NEUTRAL$ink
      if (w >= 6) {
        text(PB_X(left + w / 2), PB_Y(ypos[t]), sprintf("%.0f", w), cex = pt(6.6), col = txt_col)
      } else {
        text(PB_X(left + w / 2), PB_Y(ypos[t] + 0.36), sprintf("%.1f", w), cex = pt(5.2),
             col = NEUTRAL$ink, adj = c(0.5, 0))
      }
      left <- left + w
    }
    text(PB_X(-1), PB_Y(ypos[t]), sprintf("%s\n(n = %s)", TIER_LAB[t], fmt_int(ntier[t])),
         adj = c(1, 0.5), cex = pt(6.8))
  }
  for (v in c(0, 25, 50, 75, 100)) {
    text(PB_X(v), PB_Y(-0.6) - 2.0, paste0(v, "%"), cex = pt(6.8), adj = c(0.5, 1))
  }
  text(PB_X(50), PB_Y(-0.6) - 6.4, "Genomes in each guild (%)", cex = pt(7.4), adj = c(0.5, 1))
  for (k in seq_along(GUILDS)) {
    col_i <- (k - 1) %% 2; row_i <- (k - 1) %/% 2
    xk <- 40 + col_i * 34; yk <- 22 + BC_DY - row_i * 3.6
    rect(xk, yk - 1.2, xk + 2.4, yk + 1.2, col = GUILD_COL[GUILDS[k]], border = NA)
    text(xk + 3.4, yk, GUILD_LAB[GUILDS[k]], cex = pt(6.6), adj = c(0, 0.5))
  }

  # ---- c: raw vs within-lineage odds ratios ----
  ylab <- full$plot_label
  n_lab <- length(ylab)
  yy <- setNames(n_lab - seq_len(n_lab), ylab)  # first row at the top
  PC_Y <- function(y) 34 + BC_DY + (y + 0.6) * (40 / n_lab)
  abline_mm <- function(xv, col, lty = 1, lwd = 0.5) {
    segments(xv, 34 + BC_DY, xv, 74 + BC_DY, col = col, lty = lty, lwd = lwd / 0.75)
  }
  abline_mm(PC_X(1), NEUTRAL$grey, lty = 2, lwd = 0.5)
  for (i in seq_len(n_lab)) {
    r <- full[i, ]
    y <- PC_Y(yy[r$plot_label])
    lo <- PC_X(r$OR_glmm_low); hi <- PC_X(r$OR_glmm_high); m <- PC_X(r$OR_glmm)
    segments(lo, y, hi, y, col = NEUTRAL$ink, lwd = 0.6 / 0.75)
    segments(lo, y - 0.8, lo, y + 0.8, col = NEUTRAL$ink, lwd = 0.6 / 0.75)
    segments(hi, y - 0.8, hi, y + 0.8, col = NEUTRAL$ink, lwd = 0.6 / 0.75)
    points(PC_X(r$OR_fisher), y, pch = 5, cex = 0.9, col = NEUTRAL$ink, lwd = 0.6 / 0.75)
    face <- if (r$nonconv) "white" else if (r$sig) PRODUCER$dark else NEUTRAL$light
    points(m, y, pch = 21, cex = 1.0, col = NEUTRAL$ink, bg = face, lwd = 0.5 / 0.75)
  }
  for (v in c(0.25, 0.5, 1, 2, 4)) {
    text(PC_X(v), PC_Y(-0.6) - 2.0, format(v), cex = pt(6.8), adj = c(0.5, 1))
  }
  text(PC_X(1), 34 + BC_DY - 6.6, "Odds ratio, core/broad vs host-specific\n(log scale)",
       cex = pt(7.0), adj = c(0.5, 1))
  for (i in seq_len(n_lab)) {
    text(PC_X(0.2) - 1.0, PC_Y(yy[ylab[i]]), ylab[i], adj = c(1, 0.5), cex = pt(6.8), xpd = NA)
  }
  legend(139, 34 + BC_DY - 17,
         legend = c("Unadjusted (Fisher)", "Within lineages (GLMM)"),
         pch = c(5, 21), pt.bg = c(NA, NEUTRAL$light), col = NEUTRAL$ink, bty = "n",
         cex = pt(6.6), pt.cex = 0.9, y.intersp = 1.2, x.intersp = 0.6, xpd = NA)

  # Panel tags
  text(1, 169, "a", cex = pt(11), font = 2, adj = c(0, 1))
  text(1, 76 + BC_DY, "b", cex = pt(11), font = 2, adj = c(0, 1))
  text(97, 76 + BC_DY, "c", cex = pt(11), font = 2, adj = c(0, 1))
}

# Genus dot plot: circle = first share, diamond = second share, counts at right
genus_plot <- function(tbl, kind, x_data, ytop, step = 7.0) {
  col   <- if (kind == "p") PRODUCER$dark else SALVAGER$dark
  a_key <- if (kind == "p") "f_complete" else "f_salv"
  b_key <- if (kind == "p") "f_thiC" else "f_tr"
  n <- nrow(tbl)
  ylo <- ytop - step * (n - 1) - 3.5
  yhi <- ytop + 3.5
  x_mm <- 10 + x_data * S_X
  y_mm <- 79 + (ylo + 3) * S_Y
  h_mm <- (yhi - ylo) * S_Y
  gx <- function(v) x_mm + LABEL_MM + v / 100 * PLOT_MM
  gy <- function(y) y_mm + (y - ylo) / (yhi - ylo) * h_mm
  x0 <- x_mm + LABEL_MM
  rect(x0, y_mm, x0 + PLOT_MM, y_mm + h_mm,
       col = if (kind == "p") "#F6DFD9" else "#D9EEEC", border = NA)
  for (v in c(0, 50, 100)) segments(gx(v), y_mm, gx(v), y_mm + h_mm, col = "#E6E4DD", lwd = 0.4 / 0.75)
  ys <- ytop - step * (seq_len(n) - 1)
  for (i in seq_len(n)) {
    r <- tbl[i, ]
    a <- 100 * r[[a_key]]; b <- 100 * r[[b_key]]
    y <- gy(ys[i])
    if (!is.na(a) && !is.na(b)) segments(gx(b), y, gx(a), y, col = NEUTRAL$light, lwd = 0.8 / 0.75)
    if (!is.na(a)) points(gx(a), y, pch = 21, cex = 0.8, col = col,
                          bg = if (isTRUE(r$named)) col else "white", lwd = 0.7 / 0.75)
    if (!is.na(b)) points(gx(b), y, pch = 23, cex = 0.7, col = col, bg = "white", lwd = 0.7 / 0.75)
    text(gx(104), y, sprintf("n = %s; %s hosts", fmt_int(r$n), r$n_hosts), cex = pt(5),
         col = NEUTRAL$grey, adj = c(0, 0.5), xpd = NA)
    text(x0 - 1, y, r$genus, cex = pt(5.8), font = 4, col = NEUTRAL$ink, adj = c(1, 0.5), xpd = NA)
  }
  for (v in c(0, 50, 100)) {
    text(gx(v), y_mm - 1.5, paste0(v, "%"), cex = pt(5), adj = c(0.5, 1))
  }
  text(gx(50), y_mm - 4.2,
       sprintf("%% of genomes (circle: %s; diamond: %s)",
               if (kind == "p") "complete route" else "salvage",
               if (kind == "p") "thiC+" else "transporter"),
       cex = pt(5), adj = c(0.5, 1))
}

# ---- Export -----------------------------------------------------------------------
png(file.path(OUT_DIR, "Fig_thiamine_guilds.png"), width = FIG_W, height = FIG_OUT_H, units = "mm",
    res = 600, bg = "white", pointsize = 12)
draw_figure()
invisible(dev.off())

pdf(file.path(OUT_DIR, "Fig_thiamine_guilds.pdf"), width = FIG_W / 25.4, height = FIG_OUT_H / 25.4,
    pointsize = 12)
draw_figure()
invisible(dev.off())

svg(file.path(OUT_DIR, "Fig_thiamine_guilds.svg"), width = FIG_W / 25.4, height = FIG_OUT_H / 25.4,
    pointsize = 12)
draw_figure()
invisible(dev.off())
