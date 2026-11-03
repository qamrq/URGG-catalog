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
#  SCRIPT  : fig4_common_setup.R  -  Shared data and styles for Fig. 4d-f
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s63_*/63_genome_thiamine_guilds.tsv
#            annotation/results_phase3/s66_thiamine_by_host/*.csv
#  OUTPUT  : objects used by Fig4d/e/f
# ==========================================================

if (!exists("URGG_PROJECT")) {
  .urgg_script_dir <- local({                     # folder of this script (Rscript or RStudio)
    a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f[1]))
    else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) dirname(rstudioapi::getSourceEditorContext()$path)
    else getwd()
  })
  source(file.path(.urgg_script_dir, "_setup", "urgg_paths.R"))
}

# =============================================================================
# PART 1  Packages, inputs, palettes, theme, derived tables
# =============================================================================
need <- c("data.table", "ggplot2", "ggforce", "patchwork", "svglite", "ggrepel")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing R packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggforce); library(patchwork) })
HAVE_REPEL <- requireNamespace("ggrepel", quietly = TRUE)
cat(R.version.string, "| ggplot2", as.character(packageVersion("ggplot2")), "| ggforce", as.character(packageVersion("ggforce")), "\n")

# ---- inputs (relative to URGG_PROJECT; set once per machine, see config/local_paths.example.R) ----
R3 <- "annotation/results_phase3"
G63_F      <- urgg_in(file.path(R3, "s63_thiamine_guild_propagated/63_genome_thiamine_guilds.tsv"))
GLMM_F     <- urgg_in(file.path(R3, "s49_trait_presence_taxonomy_glmm/table_s49_thiamine_guild_presence_models.csv"), must_exist = FALSE)
S66_DIR    <- dirname(urgg_in(file.path(R3, "s66_thiamine_by_host/66_host_guilds_all.csv")))
SHARING_F  <- urgg_in(file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
GTDB_BAC_F <- urgg_in("annotation/gtdbtk_output_v232/gtdbtk.bac120.summary.tsv", must_exist = FALSE)
GTDB_AR_F  <- urgg_in("annotation/gtdbtk_output_v232/gtdbtk.ar53.summary.tsv", must_exist = FALSE)
LINEAGE_F  <- urgg_in(file.path(R3, "s32_taxonomy_novelty/lineage_full_classification_table.csv"), must_exist = FALSE)
GENUS_PREV_F <- Sys.getenv("GENUS_PREV_F", "")   # optional: genus, pct_samples (read-mapping prevalence)
OUT_DIR    <- urgg_out("fig4")
FONT       <- URGG_FONT
MIN_HOST_N <- as.integer(Sys.getenv("MIN_HOST_N", "100"))   # hosts shown in panels d-f
MIN_EDGE   <- as.integer(Sys.getenv("MIN_EDGE", "5"))       # min shared thiC+ species to draw a network edge
invisible(urgg_in(file.path(R3, "s66_thiamine_by_host/66_hostpair_shared_producers.csv")))

TIER_HOST <- "Host-specific (1 host species)"; TIER_SHARED <- "Shared (2-3 host species)"; TIER_CORE <- "Core / broad (>3 host species)"
TIER_LAB  <- c(setNames("Core", TIER_CORE), setNames("Shared", TIER_SHARED), setNames("Host-specific", TIER_HOST))

# ---- palettes (curated; shades within a family, no fade-to-white) ----------
PRODUCER <- c(dark = "#B5413A", mid = "#D9675D", light = "#EFA59C", wash = "#FBEDEA", line = "#E08A82")
SALVAGER <- c(dark = "#1F6F78", mid = "#2E959C", light = "#86C7C6", wash = "#E7F4F3", line = "#6FB5B5")
NEUTRAL  <- c(ink = "#2B2B2B", grey = "#8C8A84", light = "#B5B3AA", wash = "#F3F2EE", pool = "#4A4F57")
guild_palette <- c("1_de_novo" = PRODUCER[["dark"]], "2_partial_de_novo" = PRODUCER[["light"]],
                   "3_precursor_salvage" = SALVAGER[["light"]], "4_thiamine_salvage" = SALVAGER[["dark"]],
                   "5_none_detected" = NEUTRAL[["light"]])
guild_labels <- c("1_de_novo" = "Complete de novo", "2_partial_de_novo" = "Partial de novo (thiC+)",
                  "3_precursor_salvage" = "Precursor salvage", "4_thiamine_salvage" = "Thiamine salvage",
                  "5_none_detected" = "No route detected")
# Host palette: Bovidae values from the catalogue palette (ploting_instructions.md);
# Cervidae / Giraffidae / Camelidae shades added here - replace with the master palette if it differs.
species_palette <- c(
  "Ovis aries" = "#355C9A", "Pantholops hodgsonii" = "#4776B3", "Bos taurus" = "#5C8FD1", "Bos grunniens" = "#7DA8E6",
  "Capra hircus" = "#56B4B0", "Bubalus bubalis" = "#5BAE61", "Bos indicus" = "#7FC97F", "Procapra picticaudata" = "#C4C43A",
  "Bos frontalis" = "#D7D7A8",
  "Hydropotes inermis" = "#B5523B", "Cervus elaphus" = "#C9693F", "Alces alces" = "#D9824B", "Alces americanus" = "#E59D66",
  "Rangifer tarandus" = "#EDBE8E",
  "Giraffa camelopardalis" = "#8E6BB0",
  "Camelus bactrianus" = "#9C7B5B", "Camelus dromedarius" = "#C2A383")
host_family <- c(setNames(rep("Bovidae", 9), names(species_palette)[1:9]),
                 setNames(rep("Cervidae", 5), names(species_palette)[10:14]),
                 "Giraffa camelopardalis" = "Giraffidae", "Camelus bactrianus" = "Camelidae", "Camelus dromedarius" = "Camelidae")

theme_catalogue <- function(base = 7) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(text = element_text(colour = NEUTRAL[["ink"]]), axis.text = element_text(colour = NEUTRAL[["ink"]]),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E6E4DD", linewidth = 0.25),
          plot.title = element_text(face = "bold", size = base + 1), plot.tag = element_text(face = "bold", size = base + 3),
          legend.key.size = unit(3, "mm"), legend.text = element_text(size = base - 0.5),
          legend.title = element_text(size = base - 0.5, face = "bold"), plot.margin = margin(4, 6, 4, 4))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
fmtp <- function(x, d = 1) formatC(x, format = "f", digits = d)
fmtn <- function(x) formatC(x, format = "d", big.mark = ",")

# ---- inputs -----------------------------------------------------------------
g <- fread(G63_F)
sh <- fread(SHARING_F); setnames(sh, "representative", "genome")
tax <- rbindlist(lapply(c(GTDB_BAC_F, GTDB_AR_F), function(f) if (file.exists(f)) fread(f, sep = "\t", quote = "", select = c("user_genome", "classification")) else NULL))
if (nrow(tax)) {
  setnames(tax, c("genome", "classification"))
  tax[, genus := fifelse(grepl("g__", classification), sub("^.*g__([^;]*).*$", "\\1", classification), NA_character_)]
  tax[, domain := fifelse(grepl("d__Archaea", classification), "Archaea", "Bacteria")]
  tax[is.na(genus) | genus == "", genus := "(no genus)"]
  g <- merge(g, tax[, .(genome, genus, domain)], by = "genome", all.x = TRUE)
} else g[, `:=`(genus = "(no genus)", domain = "Bacteria")]
if ("host_species_list" %in% names(sh)) g <- merge(g, sh[, .(genome, host_species_list)], by = "genome", all.x = TRUE)
N_ALL <- nrow(g)
cat(sprintf("Genomes: %s\n", fmtn(N_ALL)))

# ---- catalogue-level numbers used in the schematic --------------------------
gpct <- g[, .N, by = guild][, pct := 100 * N / N_ALL]
gp <- function(x) gpct[guild == x, pct]
pct_thiC <- 100 * mean(g$hmp_denovo)
tr_by_guild <- g[, .(pct_transport = 100 * mean(transport)), by = guild]
tp <- function(x) tr_by_guild[guild == x, pct_transport]

# ---- representative genera per guild -----------------------------------------
host_count <- function(lists) {
  if (all(is.na(lists))) return(NA_integer_)
  length(unique(trimws(unlist(strsplit(lists[!is.na(lists)], "[;,|]")))))
}
gen <- g[genus != "(no genus)", .(n = .N, n_thiC = sum(hmp_denovo), n_complete = sum(guild == "1_de_novo"),
                                  n_salv = sum(guild == "4_thiamine_salvage"), n_prec = sum(guild == "3_precursor_salvage"),
                                  n_tr = sum(transport), n_hosts = host_count(host_species_list)), by = .(genus, domain)]
gen[, `:=`(f_thiC = n_thiC / n, f_complete = n_complete / n, f_salv = n_salv / n, f_tr = n_tr / n)]
gen[, named := grepl("^[A-Z][a-z]+(_[A-Z])?$", genus)]             # named genus vs placeholder (e.g. UBA1784, CAG-791)
if (nzchar(GENUS_PREV_F) && file.exists(GENUS_PREV_F)) {
  pv <- fread(GENUS_PREV_F); setnames(pv, tolower(names(pv)))
  gen <- merge(gen, pv[, .(genus, pct_samples)], by = "genus", all.x = TRUE)
} else gen[, pct_samples := NA_real_]
# Producers: bacteria only (the schematic is the bacterial ThiH-type route; archaea use Thi4), ranked by genomes
# with the COMPLETE route, then by thiC+ genomes. Salvagers: must also carry a ThiT/ThiBPQ transporter in most
# genomes, so the depicted salvage route is supported; Bacteroidota import thiamine through TonB-dependent systems
# that the KEGG transporter set does not capture, so their salvage genera are listed separately in the report.
EXCLUDE <- trimws(strsplit(Sys.getenv("EXCLUDE_GENERA", ""), ",")[[1]])         # optional manual override
pick_distinct <- function(d, k) {                                       # skip sister placeholders (Phocaeicola_A after Phocaeicola)
  d <- d[!genus %in% EXCLUDE][, base := sub("_[A-Z]+$", "", genus)]
  d[!duplicated(base)][1:min(k, .N)]
}
rep_prod <- pick_distinct(gen[n >= 20 & f_thiC >= 0.5 & domain == "Bacteria"][order(-n_complete, -n_thiC)], 4)
rep_salv <- pick_distinct(gen[n >= 20 & f_salv >= 0.6 & f_thiC < 0.05 & f_tr >= 0.5][order(-n_salv)], 3)
cat("Archaeal thiC+ genera (not shown in the bacterial schematic):\n")
print(gen[domain == "Archaea" & n >= 10 & f_thiC >= 0.5, .(genus, n, f_thiC = round(f_thiC, 2), n_hosts)])
cat("Salvage-dominated genera WITHOUT a detected ThiT/ThiBPQ transporter (likely TonB-dependent uptake; not shown):\n")
print(gen[n >= 20 & f_salv >= 0.6 & f_thiC < 0.05 & f_tr < 0.5][order(-n_salv)][1:min(5, .N), .(genus, n, f_salv = round(f_salv, 2), f_tr = round(f_tr, 2), n_hosts)])
fwrite(gen[order(-n_thiC)], file.path(OUT_DIR, "67_genus_thiamine_summary.csv"))
cat("Representative producer genera:\n"); print(rep_prod[, .(genus, n, f_complete = round(f_complete, 2), f_thiC = round(f_thiC, 2), n_hosts, named)])
cat("Representative salvager genera:\n"); print(rep_salv[, .(genus, n, f_salv = round(f_salv, 2), f_tr = round(f_tr, 2), n_hosts, named)])

ruc <- NULL
if (file.exists(LINEAGE_F)) {
  lin <- unique(fread(LINEAGE_F)$genome)
  r <- g[genome %chin% lin]
  if (nrow(r)) ruc <- list(n = nrow(r), tr = 100 * mean(r$transport), tpk = 100 * mean(r$tpk_thiN),
                           thiL = 100 * mean(r$tmp_kinase_thiL), salv = 100 * mean(r$guild == "4_thiamine_salvage"),
                           n_hosts = host_count(r$host_species_list))
}

# ---- host tables shared by panels d, e and f (previously defined only in panel d) ----
hg   <- fread(file.path(S66_DIR, "66_host_guilds_all.csv"))
pp   <- fread(file.path(S66_DIR, "66_hostpair_shared_producers.csv"))
keep <- hg[n_genomes >= MIN_HOST_N, host]

# ---- export: PNG 600 dpi + PDF + SVG ----
save_all <- function(p, name, w_mm, h_mm) urgg_save(p, OUT_DIR, name, w_mm, h_mm)
