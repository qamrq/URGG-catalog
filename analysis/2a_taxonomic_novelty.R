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
#  SCRIPT  : 2a_taxonomic_novelty.R  -  Taxonomic resolution and novelty
#  AUTHOR  : QRQ
#  INPUT   : annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv
#            annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv
#  OUTPUT  : annotation/results_phase3/s32_taxonomy_novelty/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})

GTDB_F   <- paste0(URGG_PROJECT, "/annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv")
PROV_F   <- paste0(URGG_PROJECT, "/annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv")
NOHIT_F  <- paste0(URGG_PROJECT, "/annotation/diamond_nr_output/no_nr_hit_gene_ids.txt")
OUT_DIR  <- paste0(URGG_PROJECT, "/annotation/results_phase3/s32_taxonomy_novelty")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(GTDB_F, PROV_F)) if (!file.exists(f)) stop("Required input missing: ", f)

strip_ext <- function(x) sub("\\.fa$", "", x)

gtdb <- read_tsv(GTDB_F, show_col_types = FALSE) %>%
  mutate(genome = strip_ext(user_genome))

RANKS <- c("domain","phylum","class","order","family","genus","species")

gtdb <- gtdb %>%
  separate(classification, into = RANKS, sep = ";", fill = "right", remove = FALSE) %>%
  mutate(across(all_of(RANKS), ~ sub("^[a-z]__", "", .x)),
         across(all_of(RANKS), ~ ifelse(is.na(.x) | .x == "", NA_character_, .x)),
         red_value = suppressWarnings(as.numeric(red_value)))

N_TOTAL <- nrow(gtdb)
cat(sprintf("\nGenomes classified: %s\n", format(N_TOTAL, big.mark = ",")))

rank_res <- tibble(
  rank = RANKS,
  n_assigned = vapply(RANKS, function(r) sum(!is.na(gtdb[[r]])), integer(1))
) %>% mutate(pct = round(100 * n_assigned / N_TOTAL, 2))

cat("\n=== Rank resolution ===\n"); print(as.data.frame(rank_res))
write_csv(rank_res, file.path(OUT_DIR, "table_rank_resolution.csv"))

gtdb <- gtdb %>%
  mutate(deepest_rank = case_when(
    !is.na(species) ~ "species", !is.na(genus) ~ "genus", !is.na(family) ~ "family",
    !is.na(order) ~ "order", !is.na(class) ~ "class", !is.na(phylum) ~ "phylum",
    !is.na(domain) ~ "domain", TRUE ~ "unclassified"))

deepest_tbl <- gtdb %>% count(deepest_rank, name = "n") %>%
  mutate(pct = round(100 * n / N_TOTAL, 2)) %>%
  arrange(match(deepest_rank, c(RANKS, "unclassified")))
cat("\n=== Deepest rank resolved per genome ===\n"); print(as.data.frame(deepest_tbl))
write_csv(deepest_tbl, file.path(OUT_DIR, "table_deepest_rank.csv"))

gtdb <- gtdb %>%
  mutate(novelty_tier = case_when(
    is.na(phylum)  ~ "Novel phylum candidate",
    is.na(class)   ~ "Novel class candidate",
    is.na(order)   ~ "Novel order candidate",
    is.na(family)  ~ "Novel family candidate",
    is.na(genus)   ~ "Novel genus candidate",
    is.na(species) ~ "Novel species candidate",
    TRUE           ~ "Matches existing GTDB species"),
    novelty_tier = factor(novelty_tier, levels = c(
      "Matches existing GTDB species", "Novel species candidate",
      "Novel genus candidate", "Novel family candidate",
      "Novel order candidate", "Novel class candidate", "Novel phylum candidate")))

novelty_tbl <- gtdb %>% count(novelty_tier, name = "n") %>%
  mutate(pct = round(100 * n / N_TOTAL, 2))
cat("\n=== Novelty tiers ===\n"); print(as.data.frame(novelty_tbl))
write_csv(novelty_tbl, file.path(OUT_DIR, "table_novelty_tiers.csv"))

cat(sprintf("\nHEADLINE for 3.2:\n  no species-level match: %s (%.1f%%)\n  no genus-level match:   %s (%.1f%%)\n",
    format(sum(is.na(gtdb$species)), big.mark=","), 100*mean(is.na(gtdb$species)),
    format(sum(is.na(gtdb$genus)),   big.mark=","), 100*mean(is.na(gtdb$genus))))

method_tbl <- gtdb %>% count(classification_method, name = "n") %>%
  mutate(pct = round(100 * n / N_TOTAL, 2)) %>% arrange(desc(n))
cat("\n=== Classification method ===\n"); print(as.data.frame(method_tbl))
write_csv(method_tbl, file.path(OUT_DIR, "table_classification_method.csv"))

red_summary <- gtdb %>% filter(!is.na(red_value)) %>%
  group_by(novelty_tier) %>%
  summarise(n = n(), red_median = round(median(red_value), 4),
            red_min = round(min(red_value), 4), red_max = round(max(red_value), 4),
            .groups = "drop")
cat(sprintf("\n=== RED values (n = %s with a RED value; %s assigned by ANI) ===\n",
            format(sum(!is.na(gtdb$red_value)), big.mark=","),
            format(sum(is.na(gtdb$red_value)), big.mark=",")))
print(as.data.frame(red_summary))
write_csv(red_summary, file.path(OUT_DIR, "table_red_by_novelty.csv"))

for (r in c("phylum", "class", "genus")) {
  tbl <- gtdb %>% filter(!is.na(.data[[r]])) %>% count(.data[[r]], name = "n") %>%
    rename(taxon = 1) %>% mutate(pct = round(100 * n / N_TOTAL, 2)) %>% arrange(desc(n))
  write_csv(tbl, file.path(OUT_DIR, sprintf("table_composition_%s.csv", r)))
  cat(sprintf("\n=== Top 10 %s ===\n", r)); print(as.data.frame(head(tbl, 10)))
}

dom_tbl <- gtdb %>% count(domain, name = "n") %>% mutate(pct = round(100*n/N_TOTAL, 2))
cat("\n=== Domain ===\n"); print(as.data.frame(dom_tbl))

prov <- read_tsv(PROV_F, show_col_types = FALSE)

rep_host <- prov %>% filter(genome == final_representative) %>%
  distinct(genome, species) %>%
  transmute(genome, host_species = gsub("_", " ", species))

gtdb_host <- gtdb %>% inner_join(rep_host, by = "genome")
cat(sprintf("\nGenomes joined to a host: %s of %s\n",
            format(nrow(gtdb_host), big.mark=","), format(N_TOTAL, big.mark=",")))

NOVEL_GENUS_OR_DEEPER <- c("Novel genus candidate", "Novel family candidate",
                           "Novel order candidate", "Novel class candidate",
                           "Novel phylum candidate")

novelty_by_host <- gtdb_host %>%
  group_by(host_species) %>%
  summarise(n_genomes = n(),
            n_novel_species = sum(novelty_tier != "Matches existing GTDB species"),
            n_novel_genus   = sum(novelty_tier %in% NOVEL_GENUS_OR_DEEPER),
            .groups = "drop") %>%
  mutate(pct_novel_species = round(100 * n_novel_species / n_genomes, 1),
         pct_novel_genus   = round(100 * n_novel_genus   / n_genomes, 1)) %>%
  arrange(desc(n_genomes))

cat("\n=== Novelty by host species ===\n"); print(as.data.frame(novelty_by_host))
write_csv(novelty_by_host, file.path(OUT_DIR, "table_novelty_by_host.csv"))

novelty_by_phylum <- gtdb %>% filter(!is.na(phylum)) %>%
  group_by(phylum) %>%
  summarise(n_genomes = n(),
            n_novel = sum(novelty_tier != "Matches existing GTDB species"),
            .groups = "drop") %>%
  mutate(pct_novel = round(100 * n_novel / n_genomes, 1)) %>%
  filter(n_genomes >= 20) %>% arrange(desc(pct_novel))
print(as.data.frame(head(novelty_by_phylum, 15)))
write_csv(novelty_by_phylum, file.path(OUT_DIR, "table_novelty_by_phylum.csv"))

if (file.exists(NOHIT_F)) {
  nohit <- read_lines(NOHIT_F)
  nohit_genome <- sub("~.*$", "", nohit)
  novel_gene_counts <- as.data.frame(table(nohit_genome), stringsAsFactors = FALSE)
  names(novel_gene_counts) <- c("genome", "n_novel_genes")

  gtdb_genes <- gtdb %>%
    left_join(novel_gene_counts, by = "genome") %>%
    mutate(n_novel_genes = replace_na(n_novel_genes, 0L),
           taxonomically_novel = novelty_tier != "Matches existing GTDB species")

  gene_summary <- gtdb_genes %>%
    group_by(taxonomically_novel) %>%
    summarise(n_genomes = n(),
              median_novel_genes = median(n_novel_genes),
              mean_novel_genes = round(mean(n_novel_genes), 1),
              .groups = "drop")
  print(as.data.frame(gene_summary))

  if (n_distinct(gtdb_genes$taxonomically_novel) == 2) {
    wt <- wilcox.test(n_novel_genes ~ taxonomically_novel, data = gtdb_genes)
    cat(sprintf("\nWilcoxon rank-sum: W = %s, p = %s\n",
                format(wt$statistic), format.pval(wt$p.value)))
  }

  by_tier <- gtdb_genes %>% group_by(novelty_tier) %>%
    summarise(n_genomes = n(), median_novel_genes = median(n_novel_genes),
              mean_novel_genes = round(mean(n_novel_genes), 1), .groups = "drop")
  cat("\n--- by novelty tier ---\n"); print(as.data.frame(by_tier))
  write_csv(by_tier, file.path(OUT_DIR, "table_novelgenes_by_noveltytier.csv"))
  write_csv(gtdb_genes %>% select(genome, novelty_tier, n_novel_genes),
            file.path(OUT_DIR, "genome_novelty_and_novelgenes.csv"))
} else {
  message("\nSkipping DIAMOND-nr cross-reference - not found at ", NOHIT_F)
}

write_csv(
  gtdb %>% select(genome, classification, all_of(RANKS), deepest_rank,
                  novelty_tier, classification_method, red_value,
                  closest_genome_ani, closest_genome_af, msa_percent),
  file.path(OUT_DIR, "genome_taxonomy_novelty.csv"))

cat("\nAll outputs written to: ", OUT_DIR, "\n", sep = "")
