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
#  SCRIPT  : 14_dbcan_cazyme.R  -  CAZyme annotation (run_dbCAN 5.2.9)
#  AUTHOR  : QRQ
#  INPUT   : CATALOG/nr_gene_catalog_rep_seq.fasta (step 12)
#  OUTPUT  : ANNOT/dbcan_output/overview.txt (consensus CAZyme table)
# ==========================================================
#
# dbCAN consensus of HMMER, DIAMOND and dbCAN-sub calls on the non-redundant catalog.

source(file.path(Sys.getenv("PIPELINE_DIR", "scripts/pipeline"), "config.R"))

step_dbcan <- function(threads = THREADS) {
  out <- ensure_dir(ANNOT, "dbcan_output")
  query <- file.path(CATALOG, "nr_gene_catalog_rep_seq.fasta")

  run_tool("run_dbcan", c(query, "protein",
                          "--out_dir", out,
                          "--db_dir", DBCAN_DB,
                          "--dia_cpu", threads,
                          "--hmm_cpu", threads,
                          "--tf_cpu", threads),
           env = ENV_DBCAN)
  log_msg("dbCAN done: ", file.path(out, "overview.txt"))
}

if (sys.nframe() == 0L) step_dbcan()
