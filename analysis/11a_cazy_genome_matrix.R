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
#  SCRIPT  : 11a_cazy_genome_matrix.R  -  Genome x CAZy family table from dbCAN
#  AUTHOR  : QRQ
#  INPUT   : dbCAN overview.tsv (argument 1)
#  OUTPUT  : CAZy long table (argument 2 = output folder)
# ==========================================================

suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 11a_cazy_genome_matrix.R <overview.tsv> <output_dir>")
}
overview_path <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Reading dbCAN overview.tsv from: ", overview_path)
overview <- read_tsv(overview_path, show_col_types = FALSE)

required_cols <- c("Gene ID")
missing_cols <- setdiff(required_cols, colnames(overview))
if (length(missing_cols) > 0) {
  stop("overview.tsv missing required column(s): ", paste(missing_cols, collapse = ", "))
}

overview <- overview %>%
  rename(gene_id = `Gene ID`) %>%
  mutate(
    n_tildes_check = lengths(regmatches(gene_id, gregexpr("~", gene_id, fixed = TRUE))),
    genome = ifelse(n_tildes_check > 0,
                    sub("~.*$", "", gene_id),
                    sub("_[0-9]+$", "", gene_id))
  )

n_genomes_recovered <- n_distinct(overview$genome)
message("Distinct genomes recovered after tilde+fallback extraction: ", n_genomes_recovered)
message("(Validated benchmark from original investigation: 22,823 of 25,451 --",
        " compare against this to confirm the same input/logic is being used.)")

extract_family <- function(x) {
  x <- ifelse(is.na(x) | x == "-" | x == "", NA_character_, x)
  x <- sub("\\+.*$", "", x)
  x <- sub("\\(.*$", "", x)
  x <- sub("_e[0-9]+.*$", "", x)
  str_extract(x, "^[A-Z]{2,3}[0-9]+")
}

method_cols <- c("dbCAN_hmm", "dbCAN_sub", "DIAMOND")
missing_methods <- setdiff(method_cols, colnames(overview))
if (length(missing_methods) > 0) {
  stop("overview.tsv missing expected method column(s): ",
       paste(missing_methods, collapse = ", "),
       "\nActual columns found: ", paste(colnames(overview), collapse = ", "))
}

overview <- overview %>%
  mutate(
    fam_hmm = extract_family(dbCAN_hmm),
    fam_sub = extract_family(dbCAN_sub),
    fam_dmd = extract_family(DIAMOND),
    cazy_family = coalesce(fam_hmm, fam_sub, fam_dmd)
  )

n_with_tool_hit <- overview %>%
  filter(!is.na(fam_hmm) | !is.na(fam_sub) | !is.na(fam_dmd) |
         (dbCAN_hmm != "-" & !is.na(dbCAN_hmm)) |
         (dbCAN_sub != "-" & !is.na(dbCAN_sub)) |
         (DIAMOND != "-" & !is.na(DIAMOND))) %>%
  nrow()
n_parseable <- sum(!is.na(overview$cazy_family))
pct_parseable <- if (n_with_tool_hit > 0) n_parseable / n_with_tool_hit else NA

message(sprintf("Genes with a tool hit: %d | Parseable family: %d (%.1f%%)",
                n_with_tool_hit, n_parseable, 100 * pct_parseable))

if (!is.na(pct_parseable) && pct_parseable < 0.5) {
  warning("SANITY CHECK FAILED: fewer than 50% of genes with a tool hit yielded ",
          "a parseable CAZy family. This mirrors the guard in the original ",
          "earlier CAZy tables -- STOP and investigate before ",
          "trusting this matrix. Possible causes: overview.tsv column format ",
          "changed, or a different dbCAN version's output structure.")
}

long_tbl <- overview %>%
  filter(!is.na(cazy_family)) %>%
  distinct(genome, cazy_family)

n_genomes_with_cazy <- n_distinct(long_tbl$genome)
n_families <- n_distinct(long_tbl$cazy_family)
message(sprintf("Genomes with >=1 annotated CAZyme family: %d", n_genomes_with_cazy))
message(sprintf("Distinct CAZy families observed: %d", n_families))
message("(Validated benchmark: 22,823/25,451 genomes with >=1 CAZyme gene;",
        " 424 testable families in the downstream statistical analysis --",
        " 'testable' there likely reflects an additional prevalence filter",
        " applied AFTER this matrix, not necessarily all families seen here.)")

write_csv(long_tbl, file.path(outdir, "genome_cazy_long.csv"))

wide_tbl <- long_tbl %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = cazy_family, values_from = present, values_fill = 0L)

write_csv(wide_tbl, file.path(outdir, "genome_x_cazy_family_presence.csv"))

message("\nOutputs written to: ", normalizePath(outdir))
message("  genome_cazy_long.csv                  (", nrow(long_tbl), " rows)")
message("  genome_x_cazy_family_presence.csv     (", nrow(wide_tbl), " genomes x ",
        ncol(wide_tbl) - 1, " families)")
