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
#  SCRIPT  : 21a_cooperative_table.R  -  Cooperative-category table with primary-model odds ratios
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s49_trait_presence_taxonomy_glmm/table_s49_coop_categories_adjcompare_presence_models.csv
#  OUTPUT  : annotation/results_phase3/s73_coop_table8/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages(library(data.table))
P  <- function(var, default) Sys.getenv(var, unset = default)
R3 <- paste0(URGG_PROJECT, "/annotation/results_phase3")
IN_F    <- P("IN_F", file.path(R3, "s49_trait_presence_taxonomy_glmm/table_s49_coop_categories_adjcompare_presence_models.csv"))
OUT_DIR <- P("OUT_DIR", file.path(R3, "s73_coop_table8"))
if (!file.exists(IN_F)) stop("Missing: ", IN_F, "\nRun 16a_run_trait_models.sh (coop_categories_adjcompare step) first.")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_DIR, "73_run_report.txt"), split = TRUE); options(width = 220)

ADJ <- c(unadjusted = "unadjusted", size = "adjusted for log genome size", compl = "adjusted for MAG completeness",
         full = "adjusted for log genome size + MAG completeness")
ROLE <- c("Antibiotic degradation" = "Enzymatic neutralization of competitors' antimicrobials",
          "Siderophores" = "Iron scavenging / uptake", "Biofilm formation" = "Attachment / spatial organization",
          "Secretion systems" = "Cell-cell / cell-host interaction machinery", "Quorum sensing" = "Population-density signalling",
          "Extracellular proteins" = "Secreted / surface-exposed proteins")
pretty_cat <- function(x) {
  k <- tolower(gsub("[^a-z]", "", tolower(x))); lab <- names(ROLE); kk <- tolower(gsub("[^a-z]", "", tolower(lab)))
  out <- lab[match(k, kk)]
  for (i in which(is.na(out))) { h <- which(startsWith(kk, substr(k[i], 1, 6))); out[i] <- if (length(h) == 1) lab[h] else x[i] }
  out
}
fmt_p  <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(signif(p, 2), format = "fg", digits = 2)))
fmt_or <- function(o, l, h) ifelse(is.na(o), "NA", sprintf("%.2f (%.2f-%.2f)", o, l, h))

m <- fread(IN_F)
cat("Adjustment levels in file:\n"); print(m[, .N, by = .(adjustment, status_ok = status == "ok")])
m[, category := pretty_cat(trait)]
fwrite(m, file.path(OUT_DIR, "table_s73_coop_all_adjustments.csv"))

base <- m[adjustment == ADJ["unadjusted"]]
if (!nrow(base)) base <- m[adjustment == m$adjustment[1]]
tab <- base[, .(category, role = ROLE[category], core = round(pct_core, 2), shared = round(pct_shared, 2), host = round(pct_host, 2),
                OR_unadj = fmt_or(OR_fisher, OR_fisher_low, OR_fisher_high), p_unadj_BH = fmt_p(p_fisher_BH), dir_unadj = direction_fisher)]
add <- function(tag) {
  x <- m[adjustment == ADJ[tag], .(category, OR = fmt_or(OR_glmm, OR_glmm_low, OR_glmm_high), p = fmt_p(LRT_p_BH),
                                    dir = fifelse(is.na(direction_glmm), "not fitted", direction_glmm), status)]
  setnames(x, c("OR", "p", "dir", "status"), paste0(c("OR_", "LRT_p_BH_", "dir_", "status_"), tag)); x
}
for (tag in c("full", "compl", "unadjusted")) if (any(m$adjustment == ADJ[tag])) tab <- merge(tab, add(tag), by = "category", all.x = TRUE, sort = FALSE)
tab <- tab[order(match(category, names(ROLE)))]
fwrite(tab, file.path(OUT_DIR, "Table8_cooperative_categories.csv"))

cat("\nTable 8 (primary model = taxonomy + log genome size + completeness):\n")
print(tab[, intersect(c("category", "core", "shared", "host", "OR_unadj", "p_unadj_BH", "OR_full", "LRT_p_BH_full", "dir_full",
                        "OR_compl", "dir_compl", "OR_unadjusted", "dir_unadjusted"), names(tab)), with = FALSE])
bad <- m[status != "ok"]
if (nrow(bad)) { cat("\nModels not 'ok' (report as not estimable):\n"); print(bad[, .(category, adjustment, status)]) }
cat("\nDraft unadjusted ORs for comparison: antibiotic degradation 3.16; siderophores 2.76; biofilm 0.47; secretion 0.56; quorum sensing 0.62; extracellular 1.36\n")
cat("Footnote: 'Primary model: logistic GLMM with nested GTDB phylum/class/order/family/genus random intercepts and fixed effects for log",
    "genome size and MAG completeness; LRT p values BH-corrected across the six categories within each model. Unadjusted: Fisher's",
    "exact test, core/broad vs host-specific, BH across the six categories.'\n")
sink()
