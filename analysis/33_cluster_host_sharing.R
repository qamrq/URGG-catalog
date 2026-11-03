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
#  SCRIPT  : 33_cluster_host_sharing.R  -  Host-sharing tiers of species clusters
#  AUTHOR  : QRQ
#  INPUT   : drep_analysis Cdb.csv, metadata.tsv
#            annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv
#  OUTPUT  : annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#            annotation/results_phase3/s33_cluster_sharing/mag_provenance.tsv
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(stringr); library(ggplot2)
})

GTDB_F   <- paste0(URGG_PROJECT, "/annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv")
OUT_DIR  <- paste0(URGG_PROJECT, "/annotation/results_phase3/s33_cluster_sharing")

PROV_F <- file.path(OUT_DIR, "mag_provenance.tsv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(PROV_F)) {
  stop("mag_provenance.tsv not found at ", PROV_F,
       "\nRun 33a_mag_provenance.R first.")
}

strip_ext <- function(x) sub("\\.fa$", "", x)

message("Reading mag_provenance.tsv (full stage-1 + stage-2 chain) ...")
prov <- read_tsv(PROV_F, show_col_types = FALSE)

stopifnot(all(c("genome", "final_cluster", "final_representative",
               "species", "bioproject") %in% names(prov)))

prov <- prov %>%
  mutate(species = str_trim(species),
         species_disp = ifelse(!is.na(species), gsub("_", " ", species), NA_character_))

n_raw   <- n_distinct(prov$genome)
n_final <- n_distinct(prov$final_cluster)
n_meta_ok <- sum(!is.na(prov$species))

message(sprintf("Raw genomes in provenance chain: %d", n_raw))
message(sprintf("Final (95%% ANI) clusters:        %d", n_final))
message(sprintf("Rows with species metadata:      %d (%.1f%%)",
                n_meta_ok, 100 * n_meta_ok / n_raw))
if (n_meta_ok < 0.99 * n_raw) {
  warning(n_raw - n_meta_ok, " raw genomes have no species metadata and will ",
          "be dropped from host-range counts.")
}

message("Building cluster -> host species membership (TRUE membership) ...")
memb <- prov %>%
  filter(!is.na(species)) %>%
  select(genome, cluster = final_cluster, species, species_disp,
         bioproject, completeness = any_of("completeness"))

cluster_summary <- memb %>%
  group_by(cluster) %>%
  summarise(
    n_member_genomes  = n_distinct(genome),
    n_host_species    = n_distinct(species),
    n_bioprojects     = n_distinct(bioproject),
    host_species_list = paste(sort(unique(species_disp)), collapse = "; "),
    .groups = "drop"
  )

cluster_summary <- cluster_summary %>%
  mutate(sharing_category = case_when(
    n_host_species == 1 ~ "Host-specific (1 host species)",
    n_host_species <= 3 ~ "Shared (2-3 host species)",
    TRUE                ~ "Core / broad (>3 host species)"
  ),
  sharing_category = factor(sharing_category,
    levels = c("Host-specific (1 host species)",
               "Shared (2-3 host species)",
               "Core / broad (>3 host species)")),
  specificity_evidence = case_when(
    n_host_species >  1 ~ "shared",
    n_member_genomes > 1 ~ "host-restricted (informative: multiple genomes)",
    TRUE                 ~ "singleton (uninformative: 1 genome)"
  ),
  sharing_multi_study = n_host_species > 1 & n_bioprojects > 1)

reps <- prov %>%
  filter(!is.na(final_representative)) %>%
  distinct(cluster = final_cluster, representative = final_representative)

rep_meta <- prov %>%
  filter(genome == final_representative) %>%
  distinct(representative = genome, rep_species = species_disp,
           rep_completeness = completeness, rep_bioproject = bioproject)

cluster_summary <- cluster_summary %>%
  left_join(reps, by = "cluster") %>%
  left_join(rep_meta, by = "representative")

have_gtdb <- file.exists(GTDB_F)
if (have_gtdb) {
  message("Joining GTDB-Tk taxonomy ...")
  gtdb <- read_tsv(GTDB_F, show_col_types = FALSE) %>%
    select(representative = user_genome, classification) %>%
    mutate(representative = strip_ext(representative)) %>%
    separate(classification,
             into = c("domain","phylum","class","order","family","genus","species_gtdb"),
             sep = ";", fill = "right", remove = FALSE) %>%
    mutate(across(domain:species_gtdb, ~ sub("^[a-z]__", "", .x)),
           across(domain:species_gtdb, ~ ifelse(.x == "", NA_character_, .x)))
  cluster_summary <- cluster_summary %>% left_join(gtdb, by = "representative")
} else {
  message("GTDB-Tk summary not found — taxonomy columns will be omitted.")
}

write_csv(cluster_summary, file.path(OUT_DIR, "cluster_host_summary.csv"))

rep_host_range <- cluster_summary %>%
  filter(!is.na(representative)) %>%
  select(representative, cluster, n_member_genomes, n_host_species,
         host_species_list, sharing_category, specificity_evidence,
         any_of(c("phylum","class","order","family","genus","species_gtdb")))
write_csv(rep_host_range, file.path(OUT_DIR, "representative_host_range.csv"))

table2 <- cluster_summary %>%
  count(sharing_category, name = "n_clusters") %>%
  mutate(pct = round(100 * n_clusters / sum(n_clusters), 1))
write_csv(table2, file.path(OUT_DIR, "table2_sharing_categories.csv"))

evidence_tbl <- cluster_summary %>%
  count(specificity_evidence, name = "n_clusters") %>%
  mutate(pct = round(100 * n_clusters / sum(n_clusters), 1))
write_csv(evidence_tbl, file.path(OUT_DIR, "table2b_specificity_evidence.csv"))

n_species_dist <- cluster_summary %>%
  count(n_host_species, name = "n_clusters") %>% arrange(n_host_species)
write_csv(n_species_dist, file.path(OUT_DIR, "cluster_host_species_distribution.csv"))

if (have_gtdb) {
  table3 <- cluster_summary %>%
    filter(n_host_species > 3) %>%
    arrange(desc(n_host_species), desc(n_member_genomes)) %>%
    select(cluster, n_host_species, n_member_genomes,
           phylum, class, order, family, genus, species_gtdb,
           host_species_list, representative)
  write_csv(table3, file.path(OUT_DIR, "table3_core_taxa.csv"))

  core_genera <- cluster_summary %>%
    filter(n_host_species > 3, !is.na(genus)) %>%
    count(genus, sort = TRUE, name = "n_core_clusters")
  write_csv(core_genera, file.path(OUT_DIR, "table3b_core_genera.csv"))
}

comp_df <- cluster_summary %>%
  filter(!is.na(rep_completeness)) %>%
  mutate(grp = ifelse(n_host_species == 1, "host-specific", "shared"))

comp_txt <- file.path(OUT_DIR, "section31_completeness_comparison.txt")
if (n_distinct(comp_df$grp) == 2) {
  wt <- wilcox.test(rep_completeness ~ grp, data = comp_df)
  med <- comp_df %>% group_by(grp) %>%
    summarise(n = n(),
              median_completeness = round(median(rep_completeness), 2),
              mean_completeness   = round(mean(rep_completeness), 2),
              .groups = "drop")
  writeLines(c(
    "Section 3.1 — representative completeness by cluster sharing status",
    "Grouping variable: species (Latin binomial), NOT common-name host.",
    "",
    capture.output(print(as.data.frame(med))),
    "",
    sprintf("Wilcoxon rank-sum: W = %s, p = %s",
            format(wt$statistic), format.pval(wt$p.value)),
    "",
    "Interpretation: lower completeness among single-host clusters is",
    "consistent with Mash-based ANI pre-clustering under-merging partial",
    "genomes, i.e. an artefact inflating apparent host specificity rather",
    "than a biological result."
  ), comp_txt)
  print(med)
} else {
  writeLines("Could not compute: only one group present.", comp_txt)
}

theme_pub <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

p3 <- ggplot(table2, aes(sharing_category, n_clusters, fill = sharing_category)) +
  geom_col(width = .65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(n_clusters, big.mark = ","), pct)),
            vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("#8C6D31", "#2F5496", "#3C7A3C")) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "Species-level clusters by host-species breadth",
       subtitle = "95% ANI catalogue; grouped by Latin binomial",
       x = NULL, y = "Number of clusters") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(OUT_DIR, "Fig3_sharing_categories.pdf"), p3, width = 6.5, height = 4.2)
ggsave(file.path(OUT_DIR, "Fig3_sharing_categories.png"), p3, width = 6.5, height = 4.2, dpi = 300)

core_ids <- cluster_summary %>% filter(n_host_species > 3) %>% pull(cluster)
if (length(core_ids) > 0) {
  pres <- memb %>%
    filter(cluster %in% core_ids) %>%
    distinct(cluster, species_disp) %>%
    mutate(present = 1L)

  ord <- cluster_summary %>% filter(cluster %in% core_ids) %>%
    arrange(desc(n_host_species)) %>% pull(cluster)
  sp_ord <- pres %>% count(species_disp, sort = TRUE) %>% pull(species_disp)

  p4 <- pres %>%
    mutate(cluster = factor(cluster, levels = ord),
           species_disp = factor(species_disp, levels = sp_ord)) %>%
    ggplot(aes(species_disp, cluster)) +
    geom_tile(fill = "#2F5496") +
    labs(title = sprintf("Host-species occurrence of core clusters (n = %d)",
                         length(core_ids)),
         subtitle = "clusters present in >3 host species",
         x = NULL, y = "Species-level cluster") +
    theme_pub +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

  h <- max(4, min(11, length(core_ids) * 0.012))
  ggsave(file.path(OUT_DIR, "Fig4_core_cluster_host_heatmap.pdf"), p4, width = 7, height = h)
  ggsave(file.path(OUT_DIR, "Fig4_core_cluster_host_heatmap.png"), p4, width = 7, height = h, dpi = 300)
}

cat("\n=====================================================\n")
cat("SECTION 3.3 — numbers for Table 2 (species-based)\n")
cat("=====================================================\n")
print(as.data.frame(table2))
cat("\nTotal clusters: ", format(sum(table2$n_clusters), big.mark = ","), "\n", sep = "")
cat("\nHost species represented: ", n_distinct(memb$species), "\n", sep = "")
cat("\n--- specificity evidence (true membership; singletons separated) ---\n")
print(as.data.frame(evidence_tbl))

multi_study <- cluster_summary %>% filter(n_host_species > 1) %>%
  count(sharing_multi_study, name = "n_clusters") %>%
  mutate(pct = round(100 * n_clusters / sum(n_clusters), 1))
cat("\n--- of shared clusters, corroborated by >1 BioProject? ---\n")
print(as.data.frame(multi_study))
cat("\n--- clusters by number of host species ---\n")
print(as.data.frame(n_species_dist))
if (have_gtdb && exists("core_genera")) {
  cat("\n--- most frequent genera among core (>3 host species) clusters ---\n")
  print(as.data.frame(head(core_genera, 15)))
}
cat("\nOutputs written to: ", OUT_DIR, "\n", sep = "")
cat("\nCROSS-CHECK vs draft Table 2 (20,180 / 4,792 / 501, sum 25,473):\n")
cat("  If these differ, the draft numbers were computed on the 21-value\n")
cat("  common-name `host` column. Use the figures above and update Table 2,\n")
cat("  the 3.3 percentages, and the Discussion sentence together.\n")
