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
#  SCRIPT  : 35a_lineage_cazyme_fisher.R  -  CAZy family enrichment in the RUC lineage (Fisher exact)
#  AUTHOR  : QRQ
#  INPUT   : Per-genome CAZy presence table, lineage and tier table
#  OUTPUT  : Per-family Fisher test results (table)
# ==========================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })

OVERVIEW_F        <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/dbcan_output/overview.tsv"
LINEAGE_MEMBERS_F <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s32_taxonomy_novelty/lineage_full_classification_table.csv"
SHARING_F         <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv"
OUT_F             <- "/data/qrq_metagenome_pipeline/metagenome_analysis/annotation/results_phase3/s34_cazyme_core_specific/table_lineage_family_fisher.csv"

TARGET_FAMILIES <- c("GT121", "GT26", "GH190")

for (f in c(OVERVIEW_F, LINEAGE_MEMBERS_F, SHARING_F)) if (!file.exists(f)) stop("Required input missing: ", f)

# =============================================================================
# STEP 1 — genome x family presence table (tilde/no-tilde fix applied)
# =============================================================================
message("Reading dbCAN overview.tsv and extracting CAZy families ...")
overview <- read_tsv(OVERVIEW_F, show_col_types = FALSE) %>%
  rename(gene_id = `Gene ID`) %>%
  mutate(
    n_tildes_check = lengths(regmatches(gene_id, gregexpr("~", gene_id, fixed = TRUE))),
    genome = ifelse(n_tildes_check > 0, sub("~.*$", "", gene_id), sub("_[0-9]+$", "", gene_id))
  )

extract_family <- function(x) {
  x <- ifelse(is.na(x) | x == "-" | x == "", NA_character_, x)
  x <- sub("\\+.*$", "", x)
  x <- sub("\\(.*$", "", x)
  x <- sub("_e[0-9]+.*$", "", x)
  str_extract(x, "^[A-Z]{2,3}[0-9]+")
}

overview <- overview %>%
  mutate(cazy_family = coalesce(
    extract_family(dbCAN_hmm), extract_family(dbCAN_sub), extract_family(DIAMOND)
  ))

presence <- overview %>%
  filter(!is.na(cazy_family)) %>%
  distinct(genome, cazy_family)

# =============================================================================
# STEP 2 — genome universe + lineage membership
# =============================================================================
sharing <- read_csv(SHARING_F, show_col_types = FALSE) %>%
  select(genome = representative)

lineage_ids <- read_csv(LINEAGE_MEMBERS_F, show_col_types = FALSE) %>%
  distinct(genome) %>%
  pull(genome)

cat(sprintf("Lineage genomes: %d | Catalogue universe: %d\n",
            length(lineage_ids), nrow(sharing)))
if (length(lineage_ids) != 76) {
  warning(sprintf("Expected 76 lineage genomes, found %d. Confirm the lineage classification table is still correctly scoped.", length(lineage_ids)))
}
if (length(lineage_ids) == 0) stop("No lineage genome IDs loaded — check LINEAGE_MEMBERS_F.")

genome_membership <- sharing %>%
  mutate(in_lineage = genome %in% lineage_ids)

# =============================================================================
# STEP 3 — 2x2 Fisher's exact test per target family
# =============================================================================
test_family <- function(fam) {
  fam_genomes <- presence %>% filter(cazy_family == fam) %>% pull(genome)
  tbl <- genome_membership %>% mutate(has_family = genome %in% fam_genomes)

  mat <- matrix(
    c(sum(tbl$in_lineage & tbl$has_family),
      sum(tbl$in_lineage & !tbl$has_family),
      sum(!tbl$in_lineage & tbl$has_family),
      sum(!tbl$in_lineage & !tbl$has_family)),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("lineage", "rest_of_catalogue"), c("has_family", "no_family"))
  )

  ft <- fisher.test(mat)
  tibble(
    cazy_family = fam,
    n_lineage_with_family = mat["lineage", "has_family"],
    n_lineage_total = sum(tbl$in_lineage),
    pct_lineage = round(100 * mat["lineage", "has_family"] / sum(tbl$in_lineage), 1),
    n_catalogue_with_family = mat["rest_of_catalogue", "has_family"],
    n_catalogue_total = sum(!tbl$in_lineage),
    pct_catalogue = round(100 * mat["rest_of_catalogue", "has_family"] / sum(!tbl$in_lineage), 1),
    OR = unname(ft$estimate),
    OR_ci_low = ft$conf.int[1],
    OR_ci_high = ft$conf.int[2],
    p_value = ft$p.value
  )
}

results <- bind_rows(lapply(TARGET_FAMILIES, test_family)) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH)

write_csv(results, OUT_F)
cat("\n=== Lineage vs. catalogue Fisher's exact results (post genome-ID-fix) ===\n")
print(as.data.frame(results))
cat("\nWritten to: ", OUT_F, "\n", sep = "")
