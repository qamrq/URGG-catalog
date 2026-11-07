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
#  SCRIPT  : 10a_cooperative_genes_and_pathways.R  -  Cooperative-gene shares and categories, gene-level KO join, pathway and sub-pathway completeness
#  AUTHOR  : QRQ
#  INPUT   : $URGG_DB/evo_ko_list/ko_dataset_Dataset S4.xlsx
#            annotation/eggnog_output/nr_catalog.emapper.annotations
#            annotation/gene_catalog/nr_gene_catalog_cluster.tsv
#            annotation/genes/per_genome_gene_counts.tsv
#            annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv
#  OUTPUT  : annotation/results_phase3/s36_cooperative_ko_core_specific/ (incl. gene-level join checkpoint)
#            phase_3_analysis/results_phase3/s36_methanogenesis/
#  USAGE   : Rscript 10a_cooperative_genes_and_pathways.R [coop_share coop_categories coop_presence pathways subpathways]   (no argument: all parts in this order)
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")
URGG_DB <- Sys.getenv("URGG_DB", "/path/to/databases")

.parts <- commandArgs(trailingOnly = TRUE)
run_part <- function(p) {
  if (length(.parts) == 0 || p %in% .parts) { message("==> ", p); TRUE } else FALSE
}

if (run_part("coop_share")) local({
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(readxl); library(stringr)
})

KO_LIST_F      <- paste0(URGG_DB, "/evo_ko_list/ko_dataset_Dataset S4.xlsx")
EGGNOG_F       <- paste0(URGG_PROJECT, "/annotation/eggnog_output/nr_catalog.emapper.annotations")
CLUSTER_F      <- paste0(URGG_PROJECT, "/annotation/gene_catalog/nr_gene_catalog_cluster.tsv")
GENE_COUNTS_F  <- paste0(URGG_PROJECT, "/annotation/genes/per_genome_gene_counts.tsv")
SHARING_F      <- paste0(URGG_PROJECT, "/annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv")
OUT_DIR        <- paste0(URGG_PROJECT, "/annotation/results_phase3/s36_cooperative_ko_core_specific")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(KO_LIST_F, EGGNOG_F, CLUSTER_F, GENE_COUNTS_F, SHARING_F)) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

ko_list <- read_excel(KO_LIST_F, skip = 1)
req_cols <- c("KO", "Cooperative_behaviour")
if (!all(req_cols %in% names(ko_list))) {
  stop("Dataset S4 missing expected column(s): ", paste(setdiff(req_cols, names(ko_list)), collapse=", "),
       "\nActual columns: ", paste(names(ko_list), collapse=", "))
}
ko_list <- ko_list %>% mutate(KO = str_trim(KO)) %>% filter(!is.na(KO), KO != "")
cooperative_kos <- unique(ko_list$KO)
ko_to_category <- setNames(ko_list$Cooperative_behaviour, ko_list$KO)

cat(sprintf("Cooperative KO set: %d unique KOs across %d categories (%s)\n",
            length(cooperative_kos), n_distinct(ko_list$Cooperative_behaviour),
            paste(unique(ko_list$Cooperative_behaviour), collapse = ", ")))

if (length(cooperative_kos) != 579) {
  warning(sprintf("Loaded %d unique KOs; Hao et al. report 579. Check for duplicate ",
                   length(cooperative_kos)),
          "KO rows or unexpected NA/blank rows in Dataset S4 before proceeding.")
}

eggnog <- read_tsv(EGGNOG_F, comment = "##", show_col_types = FALSE)
names(eggnog)[1] <- "query"

if (!"KEGG_ko" %in% names(eggnog)) {
  stop("eggNOG output missing 'KEGG_ko' column. Actual columns: ",
       paste(names(eggnog), collapse = ", "))
}

rep_gene_cooperative <- eggnog %>%
  filter(KEGG_ko != "-", !is.na(KEGG_ko)) %>%
  mutate(ko_list_split = str_split(KEGG_ko, ",")) %>%
  mutate(ko_list_split = lapply(ko_list_split, function(x) sub("^ko:", "", x))) %>%
  mutate(
    is_cooperative = vapply(ko_list_split, function(x) any(x %in% cooperative_kos), logical(1)),
    cooperative_categories = vapply(ko_list_split, function(x) {
      hits <- unique(ko_to_category[x[x %in% cooperative_kos]])
      if (length(hits) == 0) NA_character_ else paste(hits, collapse = ";")
    }, character(1))
  ) %>%
  select(rep_gene_id = query, is_cooperative, cooperative_categories)

n_rep_genes_annotated <- nrow(rep_gene_cooperative)
n_rep_genes_cooperative <- sum(rep_gene_cooperative$is_cooperative)
cat(sprintf("NR representative genes with a KO annotation: %s\n", format(n_rep_genes_annotated, big.mark=",")))
cat(sprintf("...of which flagged cooperative: %s (%.2f%%)\n",
            format(n_rep_genes_cooperative, big.mark=","), 100*n_rep_genes_cooperative/n_rep_genes_annotated))

clusters <- fread(CLUSTER_F, header = FALSE, sep = "\t",
                   col.names = c("rep_gene_id", "member_gene_id"))
cat(sprintf("Cluster membership rows (= total raw genes): %s\n", format(nrow(clusters), big.mark=",")))

setDT(rep_gene_cooperative)
clusters <- rep_gene_cooperative[clusters, on = "rep_gene_id"]
clusters[is.na(is_cooperative), is_cooperative := FALSE]

clusters[, genome := sub("~.*$", "", member_gene_id)]

per_genome <- clusters[, .(
  n_cooperative_genes = sum(is_cooperative),
  n_genes_from_cluster_tsv = .N
), by = genome]

gene_counts <- read_tsv(GENE_COUNTS_F, col_names = c("genome", "n_genes_prodigal"), show_col_types = FALSE)

check <- per_genome %>% as_tibble() %>%
  inner_join(gene_counts, by = "genome") %>%
  mutate(mismatch = n_genes_from_cluster_tsv != n_genes_prodigal)

n_mismatch <- sum(check$mismatch)
cat(sprintf("\nGenome gene-count cross-check: %d of %d genomes match cluster.tsv vs. Prodigal count exactly\n",
            sum(!check$mismatch), nrow(check)))
if (n_mismatch > 0) {
  warning(sprintf("%d genomes have a gene-count mismatch between cluster.tsv and per_genome_gene_counts.tsv - ",
                   n_mismatch),
          "see table_gene_count_mismatches.csv for the affected genomes.")
  write_csv(check %>% filter(mismatch), file.path(OUT_DIR, "table_gene_count_mismatches.csv"))
}

per_genome_final <- check %>%
  mutate(proportion_cooperative = n_cooperative_genes / n_genes_prodigal)

sharing <- read_csv(SHARING_F, show_col_types = FALSE) %>%
  select(genome = representative, sharing_category, n_host_species)

analysis <- per_genome_final %>% inner_join(sharing, by = "genome")

n_joined <- nrow(analysis)
n_expected <- nrow(sharing)
cat(sprintf("\nGenomes with a cooperative-gene proportion, joined to sharing tier: %d of %d representative genomes\n",
            n_joined, n_expected))
if (n_joined < n_expected) {
  warning(sprintf("%d representative genomes did not join - check for the same ID-parsing issue seen with dbCAN's overview.tsv before treating this as a benign gap.",
                   n_expected - n_joined))
}

write_csv(analysis, file.path(OUT_DIR, "table_cooperative_proportion_per_genome.csv"))

mod <- glm(cbind(n_cooperative_genes, n_genes_prodigal - n_cooperative_genes) ~ sharing_category,
           data = analysis, family = quasibinomial())
omnibus <- anova(mod, test = "F")
print(omnibus)

pair <- analysis %>% filter(sharing_category %in% c("Core / broad (>3 host species)", "Host-specific (1 host species)")) %>%
  mutate(sharing_category = relevel(factor(sharing_category), ref = "Host-specific (1 host species)"))
mod2 <- glm(cbind(n_cooperative_genes, n_genes_prodigal - n_cooperative_genes) ~ sharing_category,
            data = pair, family = quasibinomial())
est <- summary(mod2)$coefficients
ci  <- suppressMessages(confint(mod2))
coef_row <- grep("sharing_category", rownames(est))

result_summary <- tibble(
  comparison = "Core/broad vs. Host-specific",
  OR = exp(est[coef_row, "Estimate"]),
  OR_ci_low = exp(ci[coef_row, 1]),
  OR_ci_high = exp(ci[coef_row, 2]),
  p_value = est[coef_row, "Pr(>|t|)"],
  median_prop_core_broad = median(analysis$proportion_cooperative[analysis$sharing_category == "Core / broad (>3 host species)"]),
  median_prop_host_specific = median(analysis$proportion_cooperative[analysis$sharing_category == "Host-specific (1 host species)"]),
  median_prop_shared = median(analysis$proportion_cooperative[analysis$sharing_category == "Shared (2-3 host species)"])
)

write_csv(result_summary, file.path(OUT_DIR, "table_cooperative_tier_comparison.csv"))
print(as.data.frame(result_summary))

cat_long <- clusters[is_cooperative == TRUE & !is.na(cooperative_categories)] %>%
  as_tibble() %>%
  mutate(cooperative_categories = str_split(cooperative_categories, ";")) %>%
  tidyr::unnest(cooperative_categories) %>%
  distinct(genome, cooperative_categories) %>%
  inner_join(sharing, by = "genome") %>%
  count(cooperative_categories, sharing_category) %>%
  left_join(sharing %>% count(sharing_category, name = "n_genomes_total"), by = "sharing_category") %>%
  mutate(pct_genomes = round(100 * n / n_genomes_total, 2))

write_csv(cat_long, file.path(OUT_DIR, "table_cooperative_category_by_tier.csv"))

cat("\nAll outputs written to: ", OUT_DIR, "\n", sep = "")
})

if (run_part("coop_categories")) local({
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(readxl); library(stringr)
})

KO_LIST_F      <- paste0(URGG_DB, "/evo_ko_list/ko_dataset_Dataset S4.xlsx")
EGGNOG_F       <- paste0(URGG_PROJECT, "/annotation/eggnog_output/nr_catalog.emapper.annotations")
CLUSTER_F      <- paste0(URGG_PROJECT, "/annotation/gene_catalog/nr_gene_catalog_cluster.tsv")
SHARING_F      <- paste0(URGG_PROJECT, "/annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv")
OUT_DIR        <- paste0(URGG_PROJECT, "/annotation/results_phase3/s36_cooperative_ko_core_specific")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(KO_LIST_F, EGGNOG_F, CLUSTER_F, SHARING_F)) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

peek <- read_excel(KO_LIST_F, col_names = FALSE, n_max = 5)
header_row <- NA_integer_
for (i in seq_len(nrow(peek))) {
  row_vals <- as.character(unlist(peek[i, ]))
  if ("KO" %in% row_vals && "Cooperative_behaviour" %in% row_vals) {
    header_row <- i
    break
  }
}
if (is.na(header_row)) {
  stop("Could not find a row containing both 'KO' and 'Cooperative_behaviour' in the first 5 rows of ",
       KO_LIST_F, ". First 5 rows read as:\n", paste(capture.output(print(peek)), collapse = "\n"),
       "\nInspect the file manually and adjust.")
}
message(sprintf("Detected header on row %d of %s", header_row, basename(KO_LIST_F)))
ko_list <- read_excel(KO_LIST_F, skip = header_row - 1)
req_cols <- c("KO", "Cooperative_behaviour")
if (!all(req_cols %in% names(ko_list))) {
  stop("Dataset S4 missing expected column(s): ", paste(setdiff(req_cols, names(ko_list)), collapse = ", "),
       "\nActual columns: ", paste(names(ko_list), collapse = ", "))
}
ko_list <- ko_list %>% mutate(KO = str_trim(KO)) %>% filter(!is.na(KO), KO != "")
cooperative_kos <- unique(ko_list$KO)
ko_to_category <- setNames(ko_list$Cooperative_behaviour, ko_list$KO)

all_categories <- sort(unique(ko_list$Cooperative_behaviour))
cat(sprintf("Cooperative KO set: %d unique KOs across %d categories (%s)\n",
            length(cooperative_kos), length(all_categories), paste(all_categories, collapse = ", ")))
if (length(cooperative_kos) != 579) {
  warning(sprintf("Loaded %d unique KOs; Hao et al. report 579. Check Dataset S4 for duplicate/NA rows.",
                  length(cooperative_kos)))
}

eggnog <- read_tsv(EGGNOG_F, comment = "##", show_col_types = FALSE)
names(eggnog)[1] <- "query"

if (!"KEGG_ko" %in% names(eggnog)) {
  stop("eggNOG output missing 'KEGG_ko' column. Actual columns: ", paste(names(eggnog), collapse = ", "))
}

rep_gene_categories <- eggnog %>%
  filter(KEGG_ko != "-", !is.na(KEGG_ko)) %>%
  mutate(ko_list_split = str_split(KEGG_ko, ",")) %>%
  mutate(ko_list_split = lapply(ko_list_split, function(x) sub("^ko:", "", x))) %>%
  mutate(
    cooperative_categories = vapply(ko_list_split, function(x) {
      hits <- unique(ko_to_category[x[x %in% cooperative_kos]])
      if (length(hits) == 0) NA_character_ else paste(hits, collapse = ";")
    }, character(1))
  ) %>%
  select(rep_gene_id = query, cooperative_categories)

cat(sprintf("NR representative genes with a KO annotation: %s\n", format(nrow(rep_gene_categories), big.mark = ",")))

clusters <- fread(CLUSTER_F, header = FALSE, sep = "\t", col.names = c("rep_gene_id", "member_gene_id"))
cat(sprintf("Cluster membership rows (= total raw genes): %s\n", format(nrow(clusters), big.mark = ",")))

setDT(rep_gene_categories)
clusters <- rep_gene_categories[clusters, on = "rep_gene_id"]
clusters[, genome := sub("~.*$", "", member_gene_id)]

gene_cat_long <- clusters[!is.na(cooperative_categories), .(genome, cooperative_categories)] %>%
  as_tibble() %>%
  mutate(cooperative_categories = str_split(cooperative_categories, ";")) %>%
  tidyr::unnest(cooperative_categories) %>%
  distinct(genome, cooperative_categories)

presence_wide <- gene_cat_long %>%
  mutate(present = 1L) %>%
  tidyr::pivot_wider(names_from = cooperative_categories, values_from = present, values_fill = 0L)

all_genomes <- unique(clusters$genome)
presence_wide <- tibble(genome = all_genomes) %>%
  left_join(presence_wide, by = "genome") %>%
  mutate(across(all_of(all_categories), ~ replace(., is.na(.), 0L)))

sharing <- read_csv(SHARING_F, show_col_types = FALSE) %>%
  select(genome = representative, sharing_category, n_host_species)

result <- presence_wide %>% inner_join(sharing, by = "genome")

n_joined <- nrow(result)
n_expected <- nrow(sharing)
cat(sprintf("\nGenomes with per-category presence data, joined to sharing tier: %d of %d representative genomes\n",
            n_joined, n_expected))
if (n_joined < n_expected) {
  warning(sprintf("%d representative genomes did not join -- check for the same ID-parsing issue ",
                  n_expected - n_joined),
          "seen with dbCAN's overview.tsv before treating this as a benign gap.")
}

sig_categories <- c("antibiotic_degradation", "siderophores", "biofilm", "secretion_system")
sig_categories <- intersect(sig_categories, all_categories)
for (cat_name in sig_categories) {
  tbl <- result %>%
    group_by(sharing_category) %>%
    summarise(pct_present = round(100 * mean(.data[[cat_name]] == 1), 1), .groups = "drop")
  cat(sprintf("\n%s:\n", cat_name))
  print(as.data.frame(tbl))
}

write_csv(result, file.path(OUT_DIR, "table_cooperative_category_presence_per_genome.csv"))
cat("\nWritten: ", file.path(OUT_DIR, "table_cooperative_category_presence_per_genome.csv"), "\n", sep = "")
cat("Columns: genome, sharing_category, n_host_species, ", paste(all_categories, collapse = ", "), "\n", sep = "")
})

if (run_part("coop_presence")) local({
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(readxl); library(stringr); library(tidyr)
})

KO_LIST_F      <- paste0(URGG_DB, "/evo_ko_list/ko_dataset_Dataset S4.xlsx")
EGGNOG_F       <- paste0(URGG_PROJECT, "/annotation/eggnog_output/nr_catalog.emapper.annotations")
CLUSTER_F      <- paste0(URGG_PROJECT, "/annotation/gene_catalog/nr_gene_catalog_cluster.tsv")
SHARING_F      <- paste0(URGG_PROJECT, "/annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv")
OUT_DIR        <- paste0(URGG_PROJECT, "/annotation/results_phase3/s36_cooperative_ko_core_specific")
COOPERATIVE_JOIN_RDS <- file.path(OUT_DIR, "gene_level_cooperative_join_CHECKPOINT.rds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (file.exists(COOPERATIVE_JOIN_RDS)) {
  clusters <- readRDS(COOPERATIVE_JOIN_RDS)
} else {
  for (f in c(KO_LIST_F, EGGNOG_F, CLUSTER_F, SHARING_F)) {
    if (!file.exists(f)) stop("Required input missing: ", f)
  }

  # ---- cooperative KO set ----
  ko_list <- read_excel(KO_LIST_F, skip = 1)
  cat("Columns found in Dataset S4 after skip=1:", paste(names(ko_list), collapse = " | "), "\n")
  req_cols <- c("KO", "Cooperative_behaviour")
  missing_cols <- setdiff(req_cols, names(ko_list))
  if (length(missing_cols) > 0) {
    stop("Dataset S4 missing expected column(s) after skip=1: ", paste(missing_cols, collapse=", "),
         "\nActual columns: ", paste(names(ko_list), collapse=", "),
         "\nThe header row may not be exactly 1 row below the title, or column names may have changed - check the file directly.")
  }
  ko_list <- ko_list %>% mutate(KO = str_trim(KO)) %>% filter(!is.na(KO), KO != "")
  cooperative_kos <- unique(ko_list$KO)
  ko_to_category <- setNames(ko_list$Cooperative_behaviour, ko_list$KO)
  cat(sprintf("Cooperative KO set: %d unique KOs\n", length(cooperative_kos)))

  eggnog <- read_tsv(EGGNOG_F, comment = "##", show_col_types = FALSE)
  names(eggnog)[1] <- "query"

  rep_gene_cooperative <- eggnog %>%
    filter(KEGG_ko != "-", !is.na(KEGG_ko)) %>%
    mutate(ko_list_split = str_split(KEGG_ko, ",")) %>%
    mutate(ko_list_split = lapply(ko_list_split, function(x) sub("^ko:", "", x))) %>%
    mutate(
      is_cooperative = vapply(ko_list_split, function(x) any(x %in% cooperative_kos), logical(1)),
      cooperative_categories = vapply(ko_list_split, function(x) {
        hits <- unique(ko_to_category[x[x %in% cooperative_kos]])
        if (length(hits) == 0) NA_character_ else paste(hits, collapse = ";")
      }, character(1))
    ) %>%
    select(rep_gene_id = query, is_cooperative, cooperative_categories)

  clusters <- fread(CLUSTER_F, header = FALSE, sep = "\t",
                     col.names = c("rep_gene_id", "member_gene_id"))

  setDT(rep_gene_cooperative)
  clusters <- rep_gene_cooperative[clusters, on = "rep_gene_id"]
  clusters[is.na(is_cooperative), is_cooperative := FALSE]
  clusters[, genome := sub("~.*$", "", member_gene_id)]

  message("Saving checkpoint to avoid re-running this join in future: ", COOPERATIVE_JOIN_RDS)
  saveRDS(clusters, COOPERATIVE_JOIN_RDS)
}

ALL_CATEGORIES <- c("siderophores", "antibiotic_degradation", "biofilm",
                     "quorum_sensing", "secretion_system", "Extracellular")

genome_category_present <- clusters[is_cooperative == TRUE & !is.na(cooperative_categories)] %>%
  as_tibble() %>%
  mutate(cooperative_categories = str_split(cooperative_categories, ";")) %>%
  tidyr::unnest(cooperative_categories) %>%
  distinct(genome, cooperative_categories) %>%
  mutate(present = TRUE)

sharing <- read_csv(SHARING_F, show_col_types = FALSE) %>%
  select(genome = representative, sharing_category, n_host_species)

full_grid <- expand_grid(genome = sharing$genome, cooperative_categories = ALL_CATEGORIES)

presence_long <- full_grid %>%
  left_join(genome_category_present, by = c("genome", "cooperative_categories")) %>%
  mutate(present = coalesce(present, FALSE)) %>%
  left_join(sharing, by = "genome")

n_missing_sharing <- sum(is.na(presence_long$sharing_category)) / length(ALL_CATEGORIES)
if (n_missing_sharing > 0) {
  warning(sprintf("%d genomes did not join to a sharing_category - check for the same ID-parsing issue seen with dbCAN's overview.tsv.", n_missing_sharing))
}

write_csv(presence_long, file.path(OUT_DIR, "table_cooperative_category_presence_per_genome.csv"))

presence_wide <- presence_long %>%
  select(genome, cooperative_categories, present, sharing_category, n_host_species) %>%
  pivot_wider(names_from = cooperative_categories, values_from = present)

write_csv(presence_wide, file.path(OUT_DIR, "table_cooperative_category_presence_per_genome_WIDE.csv"))

cat(sprintf("\nWrote %s rows (long) / %s genomes (wide) to %s\n",
            format(nrow(presence_long), big.mark=","), format(nrow(presence_wide), big.mark=","), OUT_DIR))
})

if (run_part("pathways")) local({
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

BASE <- URGG_PROJECT

gtdbtk_path     <- file.path(BASE, "annotation/gtdbtk_output_v232/gtdbtk.all.summary.tsv")
eggnog_path     <- file.path(BASE, "annotation/eggnog_output/nr_catalog.emapper.annotations")
metadata_path   <- file.path(BASE, "drep_analysis/metadata.tsv")
sharing_path    <- file.path(BASE, "annotation/results_phase3/s33_cluster_sharing/representative_host_range.csv")

out_dir <- file.path(BASE, "phase_3_analysis/results_phase3/s36_methanogenesis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache_path <- file.path(out_dir, "eggnog_ko_long_cache.rds")

n_total_genomes <- 25451L

extract_genome_id <- function(query_id) {
  has_tilde <- str_detect(query_id, "~")
  out <- character(length(query_id))
  out[has_tilde]  <- str_extract(query_id[has_tilde], "^[^~]+")
  out[!has_tilde] <- str_remove(query_id[!has_tilde], "_[0-9]+$")
  out
}
strip_ext <- function(x) str_remove(x, "\\.(fa|fasta|fna)$")

marker_kos <- list(
  methanogenesis = c("K00399", "K00401", "K00402"),
  vfa_scfa = c(
    "K01847", "K01848", "K01849", "K01960", "K01026",
    "K00248", "K00929", "K01034", "K01035", "K00634",
    "K00925", "K00625"
  ),
  vitamin_cofactor = c(
    "K00589", "K02224", "K00768", "K02233",
    "K03147", "K00788", "K00949"
  ),
  nitrogen_metabolism = c(
    "K01428", "K01429", "K01430", "K01915",
    "K00265", "K00266", "K00370", "K00362"
  ),
  amino_acid_biosynthesis = c(
    "K01652", "K00053", "K01714", "K00928", "K01656"
  )
)

if (file.exists(cache_path)) {
  cat(sprintf("Loading cached genome x KO table from %s\n", cache_path))
  eggnog_long <- readRDS(cache_path)
} else {
  eggnog <- fread(eggnog_path, sep = "\t", header = TRUE,
                  na.strings = c("", "-", "NA"),
                  select = c("#query", "KEGG_ko"))
  query_col <- names(eggnog)[1]
  eggnog[, genome_id := extract_genome_id(get(query_col))]
  eggnog <- eggnog[!is.na(KEGG_ko)]

  all_markers <- unique(unlist(marker_kos))
  eggnog_long <- eggnog[, .(ko = str_remove(unlist(str_split(KEGG_ko, ",")), "^ko:")),
                         by = .(genome_id)]
  eggnog_long <- eggnog_long[ko %in% all_markers]
  eggnog_long <- unique(eggnog_long)

  saveRDS(eggnog_long, cache_path)
  cat(sprintf("Cached to %s (%d rows)\n", cache_path, nrow(eggnog_long)))
}

gtdb <- fread(gtdbtk_path, sep = "\t", header = TRUE, na.strings = c("", "NA", "N/A"))
lineage_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species")
gtdb[, (lineage_cols) := tstrsplit(classification, ";", fixed = FALSE)]
for (col in lineage_cols) gtdb[[col]] <- str_remove(gtdb[[col]], "^[a-z]__")

methanogen_genera <- c(
  "Methanobrevibacter", "Methanobrevibacter_A", "Methanobrevibacter_B",
  "Methanosphaera", "Methanomassiliicoccus", "Methanomethylophilus",
  "Methanomassiliicoccaceae", "UBA71", "Ca_Methanomethylophilus",
  "Methanocorpusculum", "Methanoculleus", "Methanomicrobium",
  "Methanimicrococcus", "Methanogenium", "Methanosarcina", "Methanosalsum",
  "Methanobacterium", "Methanobacterium_A", "Methanobacterium_B",
  "Methanobacterium_D", "Methanothermobacter", "Methanothrix",
  "Methanomethylicus", "Methanocatella", "Methanoprimaticola",
  "Methanarcanum", "Methanarmilla"
)
methanogen_families <- c(
  "Methanobacteriaceae", "Methanomassiliicoccaceae", "Methanomicrobiaceae",
  "Methanocorpusculaceae", "Methanosarcinaceae", "Methanomethylophilaceae",
  "Methanomethylicaceae", "Methanotrichaceae"
)
gtdb[, is_archaea := domain == "Archaea"]
gtdb[, is_methanogen := is_archaea &
       (genus %in% methanogen_genera | family %in% methanogen_families)]
gtdb[, join_id := strip_ext(user_genome)]

meta <- fread(metadata_path, sep = "\t", header = TRUE)
meta[, join_id := strip_ext(genome)]
meta[, host_species := str_replace_all(species, "_", " ")]

sharing <- fread(sharing_path, sep = ",", header = TRUE)
print(unique(sharing$sharing_category))
sharing[, join_id := strip_ext(representative)]

base <- gtdb[, .(join_id, user_genome, domain, phylum, class, order, family, genus, is_methanogen)]
base <- merge(base, meta[, .(join_id, host_species, bioproject, sample_site, quality_tier,
                              completeness, contamination)], by = "join_id", all.x = TRUE)
base <- merge(base, sharing[, .(join_id, sharing_category, n_host_species, cluster)],
              by = "join_id", all.x = TRUE)

cat(sprintf("\nBase genome table: %d rows (expect %d)\n", nrow(base), n_total_genomes))
print(sum(is.na(base$sharing_category)))
print(sum(is.na(base$host_species)))

for (cat_name in names(marker_kos)) {
  kos <- marker_kos[[cat_name]]
  n_required <- length(kos)
  hits_per_genome <- eggnog_long[ko %in% kos, .(n_found = uniqueN(ko)), by = genome_id]
  setnames(hits_per_genome, "genome_id", "join_id")

  base <- merge(base, hits_per_genome, by = "join_id", all.x = TRUE)
  base[is.na(n_found), n_found := 0L]
  frac_col <- paste0(cat_name, "_completeness")
  base[[frac_col]] <- base$n_found / n_required
  base[, n_found := NULL]
}

print(names(base)[str_detect(names(base), "_completeness$")])

fwrite(base, file.path(out_dir, "36e_genome_completeness_hosttier.tsv"), sep = "\t")

})

if (run_part("subpathways")) local({
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

BASE <- URGG_PROJECT
out_dir <- file.path(BASE, "phase_3_analysis/results_phase3/s36_methanogenesis")
cache_path <- file.path(out_dir, "eggnog_ko_long_cache.rds")
base_path  <- file.path(out_dir, "36e_genome_completeness_hosttier.tsv")

stopifnot(file.exists(cache_path), file.exists(base_path))

eggnog_long <- readRDS(cache_path)
base <- fread(base_path, sep = "\t", header = TRUE)
cat(sprintf("Loaded base table: %d genomes; cached KO table: %d rows\n",
            nrow(base), nrow(eggnog_long)))

old_cols <- names(base)[str_detect(names(base), "_completeness$")]
base[, (old_cols) := NULL]

subpathways <- list(
  propionate_succinate = c("K01847", "K01848", "K01849", "K01960", "K01026"),

  butyrate_kinase_route      = c("K00248", "K00929", "K00634"),
  butyrate_transferase_route = c("K00248", "K01034", "K01035"),

  # --- Acetate ---
  acetate_production = c("K00925", "K00625"),

  cobalamin_biosynthesis = c("K00589", "K02224", "K00768", "K02233"),
  thiamine_biosynthesis  = c("K03147", "K00788", "K00949"),

  # --- Nitrogen: three independent processes ---
  urea_utilization        = c("K01428", "K01429", "K01430"),
  ammonia_assimilation    = c("K01915", "K00265", "K00266"),
  nitrate_reduction       = c("K00370", "K00362"),

  bcaa_biosynthesis        = c("K01652", "K00053"),
  lysine_biosynthesis      = c("K01714", "K00928"),
  tryptophan_biosynthesis_SINGLE_MARKER_ONLY = c("K01656"),

  methanogenesis = c("K00399", "K00401", "K00402")
)

print(sapply(subpathways, length))

for (sp_name in names(subpathways)) {
  kos <- subpathways[[sp_name]]
  n_required <- length(kos)
  hits <- eggnog_long[ko %in% kos, .(n_found = uniqueN(ko)), by = genome_id]
  setnames(hits, "genome_id", "join_id")
  base <- merge(base, hits, by = "join_id", all.x = TRUE)
  base[is.na(n_found), n_found := 0L]
  col <- paste0(sp_name, "_completeness")
  base[[col]] <- base$n_found / n_required
  base[, n_found := NULL]
}

fwrite(base, file.path(out_dir, "36h_genome_subpathway_completeness.tsv"), sep = "\t")

TIER_HOST   <- "Host-specific (1 host species)"
TIER_SHARED <- "Shared (2-3 host species)"
TIER_CORE   <- "Core / broad (>3 host species)"
tier_levels <- c(TIER_HOST, TIER_SHARED, TIER_CORE)
base[, sharing_category := factor(sharing_category, levels = tier_levels)]

thresholds <- c(0.5, 0.75, 1.0)

run_tests_one_category <- function(d, col, cat_name, thresholds) {
  rbindlist(lapply(thresholds, function(thr) {
    d[, carries := get(col) >= thr]
    tab <- table(d$sharing_category, d$carries)

    if (ncol(tab) < 2 || nrow(tab) < 2 || nrow(d) == 0) {
      return(data.table(
        category = cat_name, threshold = thr, n_genomes_tested = nrow(d),
        omnibus_test = NA_character_, omnibus_p = NA_real_, min_expected = NA_real_,
        core_pct = NA_real_, hostspec_pct = NA_real_,
        OR_core_vs_hostspec = NA_real_, OR_ci_low = NA_real_, OR_ci_high = NA_real_,
        directional_p = NA_real_, note = "degenerate table"
      ))
    }

    csq <- suppressWarnings(chisq.test(tab))
    min_exp <- min(csq$expected)
    if (min_exp < 5) {
      ft <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
      omni_test <- "Fisher (simulated; expected<5)"
      omni_p <- ft$p.value
    } else {
      omni_test <- "Chi-square"
      omni_p <- csq$p.value
    }

    sub <- d[sharing_category %in% c(TIER_CORE, TIER_HOST)]
    sub[, tier2 := factor(as.character(sharing_category), levels = c(TIER_HOST, TIER_CORE))]
    tab2 <- table(sub$tier2, sub$carries)
    if (ncol(tab2) < 2 || nrow(tab2) < 2) {
      or <- ci_l <- ci_h <- dir_p <- NA_real_
    } else {
      ft2 <- fisher.test(tab2)
      or <- unname(ft2$estimate); ci_l <- ft2$conf.int[1]; ci_h <- ft2$conf.int[2]
      dir_p <- ft2$p.value
    }

    core_pct <- 100 * mean(d[sharing_category == TIER_CORE]$carries)
    host_pct <- 100 * mean(d[sharing_category == TIER_HOST]$carries)

    data.table(
      category = cat_name, threshold = thr, n_genomes_tested = nrow(d),
      omnibus_test = omni_test, omnibus_p = omni_p, min_expected = round(min_exp, 1),
      core_pct = round(core_pct, 2), hostspec_pct = round(host_pct, 2),
      OR_core_vs_hostspec = round(or, 3), OR_ci_low = round(ci_l, 3), OR_ci_high = round(ci_h, 3),
      directional_p = dir_p, note = ""
    )
  }))
}

results <- rbindlist(lapply(names(subpathways), function(sp_name) {
  col <- paste0(sp_name, "_completeness")

  if (sp_name == "methanogenesis") {
    d <- base[is_methanogen == TRUE]
    cat(sprintf("\nmethanogenesis: restricting test to %d methanogen genomes (was 25,451)\n", nrow(d)))
  } else {
    d <- base
  }
  run_tests_one_category(d, col, sp_name, thresholds)
}))

results[!is.na(omnibus_p), omnibus_p_adj := p.adjust(omnibus_p, method = "BH")]
results[!is.na(directional_p), directional_p_adj := p.adjust(directional_p, method = "BH")]
results[, direction := fifelse(
  is.na(OR_ci_low) | is.na(OR_ci_high), "untested",
  fifelse(OR_ci_low > 1, "core-enriched",
  fifelse(OR_ci_high < 1, "host-specific-enriched", "no directional claim (CI spans 1)")))]

setcolorder(results, c("category", "threshold", "n_genomes_tested",
                       "omnibus_test", "omnibus_p", "omnibus_p_adj", "min_expected",
                       "core_pct", "hostspec_pct",
                       "OR_core_vs_hostspec", "OR_ci_low", "OR_ci_high",
                       "directional_p", "directional_p_adj", "direction", "note"))

print(results, nrows = 60)
fwrite(results, file.path(out_dir, "36h_subpathway_hosttier_tests.tsv"), sep = "\t")

})
