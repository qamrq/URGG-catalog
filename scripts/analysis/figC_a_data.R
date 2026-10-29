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
#  SCRIPT  : figC_a_data.R  -  Phylum x tier x CAZy carriage table (shared by Fig. 3c)
#  AUTHOR  : QRQ
#  INPUT   : C:/qrq_metagenome_pipeline/manuscript/Figures/Fig2/supp/representative_host_range.csv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/54_cazy_family_presence_propagated.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/gtdbtk.bac120.summary.tsv
#             C:/qrq_metagenome_pipeline/manuscript/Figures/Fig3/supp/gtdbtk.ar53.summary.tsv
#  OUTPUT  : none (function carriage_table() is used by the figure scripts)
# ==========================================================

FIG3_BASE <- "C:/qrq_metagenome_pipeline/manuscript/Figures"
SHARING_F <- file.path(FIG3_BASE, "Fig2/supp/representative_host_range.csv")
PRES_F    <- file.path(FIG3_BASE, "Fig3/supp/54_cazy_family_presence_propagated.tsv")
GTDB_BAC_F <- file.path(FIG3_BASE, "Fig3/supp/gtdbtk.bac120.summary.tsv")
GTDB_AR_F  <- file.path(FIG3_BASE, "Fig3/supp/gtdbtk.ar53.summary.tsv")

TIERS <- c("Core/broad", "Shared", "Host-specific")
MIN_PHYLUM_N <- 300   # phyla with fewer genomes are pooled as "Other phyla"
MIN_CELL_N   <- 10    # phylum x tier cells with fewer genomes are left blank

tier_of <- function(cat) {
  tier <- ifelse(startsWith(cat, "Core"), "Core/broad",
          ifelse(startsWith(cat, "Shared"), "Shared",
          ifelse(startsWith(cat, "Host"), "Host-specific", NA)))
  if (anyNA(tier)) stop("Unrecognised sharing_category values present")
  tier
}

read_gtdb <- function(path) {
  x <- read.delim(path, colClasses = "character", stringsAsFactors = FALSE)
  x[, c("user_genome", "classification")]
}

phylum_of <- function(cls) {
  ph <- ifelse(grepl("p__", cls), sub(".*p__([^;]*).*", "\\1", cls), "")
  ifelse(ph == "", "Unclassified", ph)
}

# One row per family x phylum x tier with carriers, genomes, percentage and group.
# families: named list of family names per group, e.g. list(`Strongest core` = c("GT19", ...))
carriage_table <- function(families) {
  fam_order <- unlist(families, use.names = FALSE)
  fam_group <- setNames(rep(names(families), lengths(families)), fam_order)

  sh <- read.csv(SHARING_F, colClasses = "character", stringsAsFactors = FALSE)
  sh <- data.frame(representative = sh$representative,
                   tier = tier_of(sh$sharing_category), stringsAsFactors = FALSE)

  tax <- rbind(read_gtdb(GTDB_BAC_F), read_gtdb(GTDB_AR_F))
  tax$phylum <- phylum_of(ifelse(is.na(tax$classification), "", tax$classification))
  g <- merge(sh, tax[, c("user_genome", "phylum")], by.x = "representative",
             by.y = "user_genome", all.x = TRUE, sort = FALSE)
  g$phylum[is.na(g$phylum)] <- "Unclassified"
  message(sprintf("Representatives: %d | without GTDB-Tk phylum: %d",
                  nrow(g), sum(g$phylum == "Unclassified")))

  counts <- sort(table(g$phylum[g$phylum != "Unclassified"]), decreasing = TRUE)
  big <- names(counts)[counts >= MIN_PHYLUM_N]
  g$phylum_show <- ifelse(g$phylum %in% big, g$phylum, "Other phyla")
  phy_levels <- c(big, "Other phyla")
  message("Phyla shown: ", paste(big, collapse = ", "))

  pres <- read.delim(PRES_F, colClasses = "character", stringsAsFactors = FALSE)
  pres <- unique(pres[pres$cazy_family %in% fam_order, c("genome", "cazy_family")])
  missing_fam <- setdiff(fam_order, pres$cazy_family)
  if (length(missing_fam)) message("No carriers found for: ", paste(missing_fam, collapse = ", "))

  den <- as.data.frame(table(phylum_show = g$phylum_show, tier = g$tier), stringsAsFactors = FALSE)
  names(den)[3] <- "n"

  res <- do.call(rbind, lapply(fam_order, function(fam) {
    carriers <- unique(pres$genome[pres$cazy_family == fam])
    carrier <- g$representative %in% carriers
    agg <- aggregate(carrier, by = list(phylum_show = g$phylum_show, tier = g$tier),
                     FUN = function(x) c(sum(x), length(x)))
    data.frame(family = fam, phylum_show = agg$phylum_show, tier = agg$tier,
               carriers = agg$x[, 1], n = agg$x[, 2], stringsAsFactors = FALSE)
  }))
  res$pct   <- 100 * res$carriers / res$n
  res$group <- unname(fam_group[res$family])

  list(res = res, den = den, phy_levels = phy_levels, fam_order = fam_order,
       fam_group = fam_group)
}

# Percentage matrix (rows = families in `rows`, columns = list of c(phylum, tier) or NULL)
# with the label text and the number of genomes per cell
cell_matrix <- function(tab, rows, cols) {
  mat <- matrix(NA_real_, nrow = length(rows), ncol = length(cols))
  lab <- matrix("", nrow = length(rows), ncol = length(cols))
  n_mat <- matrix(NA_integer_, nrow = length(rows), ncol = length(cols))
  for (r in seq_along(rows)) {
    if (is.null(rows[[r]])) next
    for (c in seq_along(cols)) {
      if (is.null(cols[[c]])) next
      hit <- tab$res$family == rows[[r]] & tab$res$phylum_show == cols[[c]][1] &
             tab$res$tier == cols[[c]][2]
      if (!any(hit)) next
      cell <- tab$res[hit, ][1, ]
      n_mat[r, c] <- cell$n
      if (cell$n >= MIN_CELL_N) {
        mat[r, c] <- cell$pct
        lab[r, c] <- sprintf("%.0f", cell$pct)
      }
    }
  }
  list(mat = mat, lab = lab, n = n_mat)
}

# Number of genomes in each phylum x tier cell
tier_count <- function(tab, phylum, tier) {
  hit <- tab$den$phylum_show == phylum & tab$den$tier == tier
  if (any(hit)) tab$den$n[hit][1] else 0L
}
