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
#  SCRIPT  : 13_eggnog_annotation.R  -  Functional annotation (eggNOG-mapper 2.1.15, DB 5.0.2)
#  AUTHOR  : QRQ
#  INPUT   : CATALOG/nr_gene_catalog_rep_seq.fasta (step 12)
#  OUTPUT  : ANNOT/eggnog_output/nr_catalog.emapper.annotations
# ==========================================================
#
# KO, KEGG pathway and module, COG, GO, EC and CAZy columns in one pass.

step_eggnog <- function(threads = THREADS) {
  out <- ensure_dir(ANNOT, "eggnog_output")
  query <- file.path(CATALOG, "nr_gene_catalog_rep_seq.fasta")

  run_tool("emapper.py", c("-i", query,
                           "--itype", "proteins",
                           "-m", "diamond",
                           "--data_dir", EGGNOG_DB,
                           "-o", "nr_catalog",
                           "--output_dir", out,
                           "--cpu", threads,
                           "--override"),
           env = ENV_EGGNOG)
  log_msg("eggNOG-mapper done: ", file.path(out, "nr_catalog.emapper.annotations"))
}
