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
#  SCRIPT  : 12_gene_catalog_mmseqs2.R  -  Non-redundant gene catalog (MMseqs2, 95% id, 90% cov)
#  AUTHOR  : QRQ
#  INPUT   : GENES/*.faa (step 11)
#  OUTPUT  : CATALOG/nr_gene_catalog_rep_seq.fasta (input for steps 13-15)
#             CATALOG/nr_gene_catalog_cluster.tsv
# ==========================================================
#
# All predicted proteins are pooled, each header prefixed with its genome name so
# every catalog member can be traced to its source genome. Clustered with
# MMseqs2 easy-linclust at 95% identity and 90% coverage.

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

step_gene_catalog <- function(threads = THREADS) {
  ensure_dir(CATALOG)
  all_prot <- file.path(CATALOG, "all_proteins.faa")
  tmp_dir  <- file.path(CATALOG, "tmp")

  faas <- list.files(GENES, pattern = "\\.faa$", full.names = TRUE)
  log_msg("Concatenating proteins from ", length(faas), " genomes")
  con <- file(all_prot, open = "w")
  on.exit(close(con), add = TRUE)
  for (faa in faas) {
    base  <- sub("\\.faa$", "", basename(faa))
    lines <- readLines(faa)
    hdr   <- startsWith(lines, ">")
    lines[hdr] <- paste0(">", base, "~", substring(lines[hdr], 2))
    writeLines(lines, con)
  }
  close(con)
  on.exit()

  log_msg("Running MMseqs2 easy-linclust (95% id, 90% cov)")
  run_tool("mmseqs", c("easy-linclust", all_prot, file.path(CATALOG, "nr_gene_catalog"),
                       tmp_dir,
                       "--min-seq-id", "0.95", "-c", "0.9", "--cov-mode", "1",
                       "--threads", threads))
  unlink(tmp_dir, recursive = TRUE)

  rep_fa <- file.path(CATALOG, "nr_gene_catalog_rep_seq.fasta")
  log_msg("Non-redundant catalog: ", rep_fa, " (",
          sum(startsWith(readLines(rep_fa), ">")), " representatives)")
}

if (sys.nframe() == 0L) step_gene_catalog()
