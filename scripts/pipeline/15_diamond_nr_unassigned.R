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
#  SCRIPT  : 15_diamond_nr_unassigned.R  -  DIAMOND blastp against NCBI nr (unassigned genes)
#  AUTHOR  : QRQ
#  INPUT   : CATALOG/nr_gene_catalog_rep_seq.fasta (step 12), NR_DMND (nr.dmnd)
#  OUTPUT  : ANNOT/diamond_nr/diamond_nr_results.tsv
# ==========================================================
#
# Searches genes without an eggNOG assignment against NCBI nr. nr is about
# 150-200 GB, and the search takes days on the full catalog.

step_diamond_nr <- function(threads = THREADS) {
  out <- ensure_dir(ANNOT, "diamond_nr")
  query <- file.path(CATALOG, "nr_gene_catalog_rep_seq.fasta")

  run_tool("diamond", c("blastp",
                        "-q", query,
                        "-d", NR_DMND,
                        "-o", file.path(out, "diamond_nr_results.tsv"),
                        "--threads", threads,
                        "--max-target-seqs", 5,
                        "--evalue", "1e-5",
                        "--block-size", 12, "--index-chunks", 1,
                        "--outfmt", 6, "qseqid", "sseqid", "pident", "length",
                        "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore", "stitle", "staxids"),
           env = ENV_DREP)
  log_msg("nr annotation done: ", file.path(out, "diamond_nr_results.tsv"))
}
