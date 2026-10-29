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
#  SCRIPT  : config.R  -  Shared configuration and helpers for pipeline steps 01-15
#  AUTHOR  : QRQ
#  INPUT   : none (edit the paths below)
#  OUTPUT  : objects and helper functions used by every pipeline step
# ==========================================================

# Run the pipeline from the repository root, e.g.
#   Rscript scripts/pipeline/run_pipeline.R
PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", "scripts/pipeline")

# ---- Root and shared settings (edit) -----------------------------------------
ROOT    <- Sys.getenv("URGG_ROOT", "/data/qrq_metagenome_pipeline")
THREADS <- 48      # threads for single multi-threaded tools (GTDB-Tk, MMseqs2, DIAMOND)
PARALLEL_JOBS <- 24  # concurrent genomes for per-genome tools (Prodigal)

# ---- Input and working folders ----------------------------------------------
RAW_DIR     <- file.path(ROOT, "raw_fastq")          # step 01 input
QC_DIR      <- file.path(ROOT, "qc_fastp")           # step 01 output
REF_DIR     <- file.path(ROOT, "hocort_index")       # step 02 host and plant references
INDEX_DIR   <- file.path(REF_DIR, "Combined_indexes")
HOCORT_DIR  <- file.path(ROOT, "hocort_clean_reads") # step 02 output
ASM_DIR     <- file.path(ROOT, "megahit_assembly")   # step 03 output
BASALT_DIR  <- file.path(ROOT, "basalt_binning")     # step 04 output
QC_OUT      <- file.path(ROOT, "checkm2_quality")    # step 05 output
DREP_DIR    <- file.path(ROOT, "metagenome_analysis", "drep_analysis")
GENOME_DIR  <- file.path(DREP_DIR, "drep_95_species_final", "dereplicated_genomes")
ANNOT       <- file.path(ROOT, "annotation")
GENES       <- file.path(ANNOT, "genes")
CATALOG     <- file.path(ANNOT, "gene_catalog")
DB          <- file.path(ROOT, "databases")

# ---- Reference databases ------------------------------------------------------
GTDBTK_DB  <- file.path(ROOT, "databases", "release226")
NR_DMND    <- file.path(DB, "ncbi_nr", "nr.dmnd")
EGGNOG_DB  <- file.path(DB, "eggnog")
DBCAN_DB   <- file.path(DB, "dbcan")

# ---- Conda environments (NULL = tool is on PATH already) -----------------------
ENV_DREP   <- "drep"
ENV_EGGNOG <- "eggnog"
ENV_DBCAN  <- "dbcan"
ENV_GTDBTK <- NULL

# ---- Helpers --------------------------------------------------------------------

# Timestamped console message
log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
}

# Run an external program. Returns its exit status.
#   tool  : program name, e.g. "fastp"
#   args  : character vector of arguments (one element per argument, no quoting)
#   env   : conda environment to run in (NULL = use PATH)
#   log   : file to capture stdout and stderr (NULL = console)
#   must  : stop on failure when TRUE
run_tool <- function(tool, args = character(), env = NULL, log = NULL, must = TRUE) {
  args <- as.character(args)
  if (!is.null(env)) {
    args <- c("run", "-n", env, tool, args)
    tool <- "conda"
  }
  out <- if (is.null(log)) "" else log
  status <- system2(tool, args, stdout = out, stderr = out)
  if (must && !identical(as.integer(status), 0L)) {
    stop(sprintf("%s failed (exit %s)%s", tool, status,
                 if (is.null(log)) "" else paste0(" - see ", log)))
  }
  invisible(status)
}

# Make a string safe for use in a file name
safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# Ensure a directory exists
ensure_dir <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(d)
}

# Number of files in a folder (0 if the folder does not exist)
count_files <- function(dir, pattern = NULL) {
  if (!dir.exists(dir)) return(0L)
  length(list.files(dir, pattern = pattern))
}
