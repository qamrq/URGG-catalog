#!/usr/bin/env Rscript
# Regenerate lineage_full_classification_table.csv from raw GTDB-Tk summaries.
#
# Authoritative definition (per handoff doc, item 8): genomes whose GTDB-Tk
# `classification` string contains "UMGS1840" or "HGM11507". Matched against
# bac120 (expected) and ar53 (safety check only — this lineage is Clostridia,
# so archaeal hits would indicate a problem, not a real result).

suppressPackageStartupMessages(library(tidyverse))

bac120_path <- "/disk252/qamar/metagenome_analysis/annotation/gtdbtk_output_v232/classify/gtdbtk.bac120.summary.tsv"
ar53_path   <- "/disk252/qamar/metagenome_analysis/annotation/gtdbtk_output_v232/classify/gtdbtk.ar53.summary.tsv"

out_dir  <- "/disk252/qamar/metagenome_analysis/annotation/results_phase3/s32_taxonomy_novelty"
out_csv  <- file.path(out_dir, "lineage_full_classification_table.csv")
out_list <- file.path(out_dir, "lineage_genome_ids.txt")

pattern <- "UMGS1840|HGM11507"

read_summary <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    select(genome = user_genome, classification)
}

bac120 <- read_summary(bac120_path)
ar53   <- read_summary(ar53_path)

lineage_bac <- bac120 %>% filter(str_detect(classification, pattern)) %>% mutate(domain = "bac120")
lineage_ar  <- ar53   %>% filter(str_detect(classification, pattern)) %>% mutate(domain = "ar53")

if (nrow(lineage_ar) > 0) {
  warning(sprintf(
    "%d archaeal genome(s) matched UMGS1840/HGM11507 — unexpected for a Clostridia lineage. Inspect manually.",
    nrow(lineage_ar)
  ))
}

lineage_full <- bind_rows(lineage_bac, lineage_ar)

cat(sprintf(
  "Matched %d genomes total (bac120: %d, ar53: %d)\n",
  nrow(lineage_full), nrow(lineage_bac), nrow(lineage_ar)
))

if (nrow(lineage_full) != 76) {
  warning(sprintf(
    "Expected 76 lineage genomes (per manuscript/handoff doc), got %d. Inspect before proceeding to the tree build.",
    nrow(lineage_full)
  ))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(lineage_full, out_csv)
write_lines(lineage_full$genome, out_list)

cat("Wrote:\n -", out_csv, "\n -", out_list, "(plain genome ID list, feeds the tree pipeline)\n")
