#!/usr/bin/env Rscript
###############################################################################
# 78_para126_plasmid_amr.R
#
# Recomputes the plasmid numbers in the AMR/mobilome paragraph (v14 para 126):
#   (a) plasmid carriage by tier (Platon; per genome; should reproduce the draft)
#   (b) AMR gene on a plasmid contig by tier, defined CONSISTENTLY with the rest of
#       the paper: an AMRFinderPlus AMR-type gene whose contig is a Platon plasmid contig.
#       The draft (2.42% vs 1.14%, OR 2.14) cannot be AMRFinderPlus-based: 1.14% of
#       host-specific genomes is more than the 201 host-specific genomes that carry
#       any AMR-type gene. It most likely used Platon's own "# AMRs" column; that
#       definition is reproduced here for comparison.
#   (c) plasmid-positive vs plasmid-negative genomes: AMR-type gene anywhere (294 definition)
#
#   Rscript 78_para126_plasmid_amr.R
###############################################################################
suppressPackageStartupMessages(library(data.table))
P  <- function(var, default) Sys.getenv(var, unset = default)
A  <- "/disk252/qamar/metagenome_analysis/annotation"
R3 <- file.path(A, "results_phase3")
SHARING_F <- P("SHARING_F", file.path(R3, "s33_cluster_sharing/representative_host_range.csv"))
AMR_F     <- P("AMR_F", file.path(R3, "s65_amr_check/65_amr_presence.tsv"))
MERGED_F  <- P("MERGED_F", file.path(A, "amrfinderplus_merged.tsv"))
PLATON_F  <- P("PLATON_F", file.path(A, "platon_merged.tsv"))
OUT_DIR   <- P("OUT_DIR", file.path(R3, "s78_para126"))
for (f in c(SHARING_F, AMR_F, MERGED_F, PLATON_F)) if (!file.exists(f)) stop("Missing: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "78_run_report.txt"), split = TRUE); options(width = 200)
TIER_HOST <- "Host-specific (1 host species)"; TIER_CORE <- "Core / broad (>3 host species)"

sh <- fread(SHARING_F, select = c("representative", "sharing_category")); setnames(sh, "representative", "genome")
amr_g <- unique(fread(AMR_F)[trait == "any_AMR", genome])
pl <- fread(PLATON_F, sep = "\t", quote = "")
amrcol <- grep("^# ?AMRs$", names(pl), value = TRUE)[1]
mg <- fread(MERGED_F, sep = "\t", quote = "")
tcol <- intersect(c("Type", "Element type"), names(mg))[1]; ccol <- intersect(c("Contig id", "Contig ID", "contig"), names(mg))[1]
am <- unique(mg[get(tcol) == "AMR", .(genome, contig = get(ccol))])
plc <- unique(pl[, .(genome, contig = ID)])
hit <- merge(am, plc, by = c("genome", "contig"))
cat(sprintf("AMR-type gene contigs: %d | on Platon plasmid contigs: %d (in %d genomes)\n", nrow(am), nrow(hit), uniqueN(hit$genome)))
if (!nrow(hit)) cat("WARNING: no AMR contig matched a Platon contig; check contig naming:\n  AMR e.g.", am$contig[1], "\n  Platon e.g.", plc$contig[1], "\n")

sh[, `:=`(plasmid = genome %chin% pl$genome, any_AMR = genome %chin% amr_g,
          amr_on_plasmid = genome %chin% hit$genome,
          platon_amr = if (!is.na(amrcol)) genome %chin% pl[get(amrcol) > 0, genome] else NA)]
tier_row <- function(v, label) {
  a <- sh[sharing_category == TIER_CORE, get(v)]; b <- sh[sharing_category == TIER_HOST, get(v)]
  ft <- fisher.test(matrix(c(sum(a), sum(!a), sum(b), sum(!b)), 2, byrow = TRUE))
  data.table(outcome = label, n_all = sum(sh[[v]]), pct_all = round(100 * mean(sh[[v]]), 2),
             pct_core = round(100 * mean(a), 2), pct_host = round(100 * mean(b), 2),
             OR = signif(unname(ft$estimate), 3), OR_lo = signif(ft$conf.int[1], 3), OR_hi = signif(ft$conf.int[2], 3), p = signif(ft$p.value, 2))
}
cat("\n(a, b) Core/broad vs host-specific\n")
tt <- rbind(tier_row("plasmid", "Plasmid contig (Platon)"),
            tier_row("amr_on_plasmid", "AMR-type gene (AMRFinderPlus) on a Platon plasmid contig"),
            if (!is.na(amrcol)) tier_row("platon_amr", "Platon '# AMRs' > 0 (draft definition, for comparison)"))
print(tt)
cat("Draft: plasmid 47.6 vs 42.0%, OR 1.25 (1.03-1.52), padj 0.058; AMR on plasmid 2.42 vs 1.14%, OR 2.14 (1.05-3.95), padj 0.058\n")

cat("\n(c) AMR-type gene anywhere: plasmid-positive vs plasmid-negative genomes\n")
ft <- fisher.test(table(factor(sh$plasmid, c(TRUE, FALSE)), factor(sh$any_AMR, c(TRUE, FALSE))))
pc <- data.table(pct_amr_plasmid_pos = round(100 * mean(sh[plasmid == TRUE, any_AMR]), 2), pct_amr_plasmid_neg = round(100 * mean(sh[plasmid == FALSE, any_AMR]), 2),
                 n_plasmid_pos = sum(sh$plasmid), OR = signif(unname(ft$estimate), 3), OR_lo = signif(ft$conf.int[1], 3), OR_hi = signif(ft$conf.int[2], 3), p = signif(ft$p.value, 2))
print(pc); cat("Draft: 2.0 vs 0.6%, OR 3.38 (2.61-4.41), p = 2.9e-23\n")
fwrite(tt, file.path(OUT_DIR, "table_s78_plasmid_amr_by_tier.csv")); fwrite(pc, file.path(OUT_DIR, "table_s78_plasmid_vs_amr.csv"))
sink()
