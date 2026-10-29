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
#  SCRIPT  : config.R  -  Shared paths, settings and helpers for pipeline steps 01-15
#  AUTHOR  : QRQ
#  INPUT   : none
#  OUTPUT  : constants and helper functions used by every pipeline step
# ==========================================================

PIPELINE_DIR <- "scripts/pipeline"

ROOT <- "/data/qrq_metagenome_pipeline"

THREADS        <- 48   # single multi-threaded tools (GTDB-Tk, MMseqs2, DIAMOND)
PARALLEL_JOBS  <- 24   # concurrent genomes for per-genome tools (Prodigal)

RAW_DIR      <- file.path(ROOT, "raw_fastq")
QC_DIR       <- file.path(ROOT, "qc_fastp")
REF_DIR      <- file.path(ROOT, "hocort_index")
INDEX_DIR    <- file.path(REF_DIR, "Combined_indexes")
HOCORT_DIR   <- file.path(ROOT, "hocort_clean_reads")
ASM_DIR      <- file.path(ROOT, "megahit_assembly")
BASALT_DIR   <- file.path(ROOT, "basalt_binning")
QC_OUT       <- file.path(ROOT, "checkm2_quality")
DREP_DIR     <- file.path(ROOT, "metagenome_analysis", "drep_analysis")
GENOME_DIR   <- file.path(DREP_DIR, "drep_95_species_final", "dereplicated_genomes")
ANNOT        <- file.path(ROOT, "annotation")
GENES        <- file.path(ANNOT, "genes")
CATALOG      <- file.path(ANNOT, "gene_catalog")
DB           <- file.path(ROOT, "databases")

GTDBTK_DB <- file.path(DB, "release226")
NR_DMND   <- file.path(DB, "ncbi_nr", "nr.dmnd")
EGGNOG_DB <- file.path(DB, "eggnog")
DBCAN_DB  <- file.path(DB, "dbcan")

# Conda environments the tools are run in (NULL = the tool is on the PATH)
ENV_DREP   <- "drep"
ENV_EGGNOG <- "eggnog"
ENV_DBCAN  <- "dbcan"
ENV_GTDBTK <- NULL

# Timestamped console message
log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
}

# Run an external program and return its exit status.
#   env  : conda environment to run in (NULL = PATH)
#   log  : file for stdout and stderr (NULL = console)
#   must : stop when the program fails
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

# Create a directory (and parents) and return its path
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
