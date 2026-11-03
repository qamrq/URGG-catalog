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
#  SCRIPT  : 32a_target_sequence_identity.R  -  Sequence identity of the unclassified sheep genome to its relatives
#  AUTHOR  : QRQ
#  INPUT   : annotation/results_phase3/s32_taxonomy_novelty/target_ml_tree_extrefs/
#  OUTPUT  : annotation/results_phase3/s32_taxonomy_novelty/target_ml_tree_extrefs/seq_identity/
# ==========================================================

URGG_PROJECT <- Sys.getenv("URGG_PROJECT", "/path/to/metagenome_analysis")

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr); library(ggplot2); library(patchwork) })

S <- paste0(URGG_PROJECT, "/annotation/results_phase3")
env <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
WORK     <- env("WORK",     file.path(S, "s32_taxonomy_novelty/target_ml_tree_extrefs"))
ALN      <- env("ALN",      file.path(WORK, "combined_trimmed.fasta"))
ROLES    <- env("ROLES",    file.path(WORK, "tip_roles.tsv"))
LABEL_F  <- env("LABEL_CSV", paste0(URGG_PROJECT, "/phase_3_analysis/scripts/results/lineage_full_classification_table.csv"))
GTDB_TAX <- env("GTDB_TAX", "")
TARGET   <- env("TARGET",   "CNP0001390__Sheep__Feces__T1900044680_TS20190146_bin.8")
OUT_DIR  <- env("OUT_DIR",  file.path(WORK, "seq_identity"))
EXPECT_COLS <- as.integer(env("EXPECT_COLS", "4930"))
WIN <- as.integer(env("WIN", "200")); STEP <- as.integer(env("STEP", "50"))
MIN_OVERLAP <- as.integer(env("MIN_OVERLAP", "500"))
MIN_OCC <- as.numeric(env("MIN_OCC", "0.5"))

stopf <- function(...) stop("CHECKPOINT FAILED: ", sprintf(...), call. = FALSE)
ok <- function(...) cat("  [ok] ", sprintf(...), "\n", sep = "")

cat("== C1 inputs ==\n")
for (f in c(ALN, ROLES, LABEL_F)) if (!file.exists(f)) stopf("missing file: %s", f)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE); ok("inputs found")

cat("== C2 alignment ==\n")
ln <- readLines(ALN); hd <- which(startsWith(ln, ">"))
ids <- sub("^>", "", ln[hd]); ids <- sub("[ \t].*$", "", ids)
seqs <- mapply(function(a, b) paste(ln[(a + 1):b], collapse = ""), hd, c(hd[-1] - 1, length(ln)))
lens <- nchar(seqs)
if (length(unique(lens)) != 1) stopf("sequences differ in length (%s)", paste(unique(lens), collapse = ","))
if (lens[1] != EXPECT_COLS) warning(sprintf("alignment width %d, expected %d (set EXPECT_COLS if intended)", lens[1], EXPECT_COLS))
A <- do.call(rbind, strsplit(toupper(seqs), "")); rownames(A) <- ids
A[A %in% c("-", ".", "X", "*", "?")] <- NA
ok("%d sequences x %d columns", nrow(A), ncol(A))

cat("== C3 groups ==\n")
Lb <- read_csv(LABEL_F, show_col_types = FALSE)
R  <- read_tsv(ROLES, show_col_types = FALSE)
if (!TARGET %in% ids) stopf("target %s not in alignment", TARGET)
grp <- setNames(rep(NA_character_, length(ids)), ids)
grp[ids %in% Lb$genome] <- "RUC lineage"
grp[ids == TARGET] <- "Target"
rr <- R$role[match(ids, R$id)]
grp[is.na(grp) & rr == "outgroup"] <- "Outgroup"
isref <- is.na(grp) & rr == "ref_neighbour"
grp[isref] <- "Reference"
if (nzchar(GTDB_TAX) && file.exists(GTDB_TAX)) {
  tx <- read_tsv(GTDB_TAX, col_names = c("acc", "tax"), col_types = "cc", progress = FALSE)
  o <- str_match(tx$tax[match(ids, tx$acc)], "o__([^;]*)")[, 2]
  grp[isref] <- ifelse(is.na(o[isref]), "Reference", paste0("GTDB ", o[isref]))
} else cat("  (GTDB_TAX not set: references pooled, order panels limited)\n")
if (anyNA(grp)) stopf("%d sequences have no group (first: %s)", sum(is.na(grp)), ids[is.na(grp)][1])
print(table(grp))
if (sum(grp == "RUC lineage") != 76) stopf("expected 76 lineage sequences, found %d", sum(grp == "RUC lineage"))
ok("76 lineage + target present")

ident <- function(a, b) { m <- !is.na(a) & !is.na(b); n <- sum(m); if (n < MIN_OVERLAP) return(c(NA_real_, n)); c(100 * mean(a[m] == b[m]), n) }

cat("== target vs all ==\n")
tv <- t(sapply(ids, function(i) ident(A[TARGET, ], A[i, ]))); colnames(tv) <- c("identity", "n_cols")
TV <- data.frame(id = ids, group = grp, tv)[ids != TARGET, ]
write.csv(TV, file.path(OUT_DIR, "identity_target_vs_all.csv"), row.names = FALSE)
if (any(is.na(TV$identity))) cat("  note: ", sum(is.na(TV$identity)), "sequences share <", MIN_OVERLAP, "columns with the target and are excluded\n")

cat("== yardstick pairs ==\n")
pair_ids <- function(g1, g2, same = (g1 == g2)) {
  a <- ids[grp %in% g1]; b <- ids[grp %in% g2]
  p <- expand.grid(a = a, b = b, stringsAsFactors = FALSE); p <- p[p$a < p$b | (!same & p$a != p$b), ]
  if (!same) p <- p[!duplicated(t(apply(p, 1, sort))), ]
  p
}
ref_grp <- function(pat) unique(grp[grepl(pat, grp)])
sets <- list(
  list("within order",  "RUC vs RUC",             "RUC lineage", "RUC lineage"),
  list("within order",  "RUC vs GTDB UMGS1840",   "RUC lineage", ref_grp("UMGS1840")),
  list("between order", "RUC vs GTDB SK-Y3",      "RUC lineage", ref_grp("SK-Y3")),
  list("between order", "RUC vs GTDB RGIG8785",   "RUC lineage", ref_grp("RGIG8785")),
  list("between order", "GTDB SK-Y3 vs RGIG8785", ref_grp("SK-Y3"), ref_grp("RGIG8785")))
set.seed(42); MAXP <- 3000
PV <- list()
for (s in sets) {
  if (!length(s[[3]]) || !length(s[[4]])) next
  p <- pair_ids(s[[3]], s[[4]], same = identical(s[[3]], s[[4]]))
  if (nrow(p) == 0) next
  if (nrow(p) > MAXP) p <- p[sample(nrow(p), MAXP), ]
  v <- apply(p, 1, function(r) ident(A[r[1], ], A[r[2], ])[1])
  PV[[s[[2]]]] <- data.frame(kind = s[[1]], comparison = s[[2]], identity = v)
}
tg <- function(lab, g) { x <- TV$identity[TV$group %in% g]; if (length(x)) data.frame(kind = "target", comparison = lab, identity = x) }
PV[["t1"]] <- tg("Target vs RUC", "RUC lineage")
PV[["t2"]] <- tg("Target vs GTDB UMGS1840", ref_grp("UMGS1840"))
PV[["t3"]] <- tg("Target vs GTDB SK-Y3", ref_grp("SK-Y3"))
PV[["t4"]] <- tg("Target vs GTDB RGIG8785", ref_grp("RGIG8785"))
PV[["t5"]] <- tg("Target vs outgroup", "Outgroup")
PD <- bind_rows(PV) %>% filter(!is.na(identity))
SM <- PD %>% group_by(kind, comparison) %>%
  summarise(n = n(), median = round(median(identity), 1), q25 = round(quantile(identity, .25), 1),
            q75 = round(quantile(identity, .75), 1), min = round(min(identity), 1), max = round(max(identity), 1), .groups = "drop")
write.csv(SM, file.path(OUT_DIR, "identity_pairs_summary.csv"), row.names = FALSE)
print(as.data.frame(SM), row.names = FALSE)
if (any(SM$n < 2)) cat("  note: some comparisons have < 2 pairs (e.g. RGIG8785 has 2 references); read those as points, not distributions\n")

cat("== sliding window ==\n")
starts <- seq(1, ncol(A) - WIN + 1, by = STEP)
tgt_groups <- c("RUC lineage", ref_grp("UMGS1840"), ref_grp("SK-Y3"), ref_grp("RGIG8785"))
tgt_groups <- unique(tgt_groups)
WI <- do.call(rbind, lapply(tgt_groups, function(g) {
  m <- ids[grp == g]
  do.call(rbind, lapply(starts, function(s) {
    cols <- s:(s + WIN - 1)
    occ <- mean(!is.na(A[TARGET, cols]))
    if (occ < MIN_OCC) return(data.frame(group = g, mid = s + WIN / 2, median = NA_real_, lo = NA_real_, hi = NA_real_, target_occupancy = occ))
    v <- vapply(m, function(i) { a <- A[TARGET, cols]; b <- A[i, cols]; k <- !is.na(a) & !is.na(b); if (sum(k) < WIN / 4) NA_real_ else 100 * mean(a[k] == b[k]) }, 0)
    data.frame(group = g, mid = s + WIN / 2, median = median(v, na.rm = TRUE), lo = quantile(v, .1, na.rm = TRUE), hi = quantile(v, .9, na.rm = TRUE), target_occupancy = occ)
  }))
}))
write.csv(WI, file.path(OUT_DIR, "window_identity.csv"), row.names = FALSE)
nm <- unique(WI$mid[is.na(WI$median)]); cat(sprintf("  %d of %d windows masked (target has < %.0f%% residues there)\n", length(nm), length(starts), 100 * MIN_OCC))

cat("== figure ==\n")
th <- theme_classic(base_size = 9) + theme(legend.position = "bottom", legend.title = element_blank())
PD$comparison <- factor(PD$comparison, levels = unique(PD$comparison))
pa <- ggplot(PD, aes(comparison, identity, fill = kind)) +
  geom_boxplot(outlier.size = .6, linewidth = .3, width = .65) +
  scale_fill_manual(values = c("within order" = "#B39DDB", "between order" = "#FDB462", target = "#E31A1C")) +
  coord_flip() + labs(x = NULL, y = "Amino-acid identity over aligned marker columns (%)", tag = "a") + th
pb <- ggplot(WI, aes(mid, median, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .12, colour = NA) + geom_line(linewidth = .5) +
  labs(x = "Alignment position (trimmed bac120 markers)", y = "Target identity per window (%)", tag = "b") + th
fig <- pa / pb + plot_layout(heights = c(1, 1))
ggsave(file.path(OUT_DIR, "fig_target_sequence_identity.pdf"), fig, width = 7.2, height = 7.2)
ggsave(file.path(OUT_DIR, "fig_target_sequence_identity.png"), fig, width = 7.2, height = 7.2, dpi = 300)
ok("wrote figure and tables to %s", OUT_DIR)

cat("\nHOW TO READ THIS (guidance, not a result):\n",
    "  - the within-order and between-order boxes are the yardsticks; the red target boxes are read against them;\n",
    "  - a target-vs-UMGS1840 box that overlaps the between-order boxes and sits below the within-order ones is\n",
    "    consistent with the target lying outside UMGS1840; it does NOT say which order it belongs to;\n",
    "  - a flat window profile means the divergence is spread across all markers (no sign of a chimeric or\n",
    "    contaminated bin); dips confined to a few windows would need checking;\n",
    "  - values are identities over conserved markers, not genome-wide AAI.\n", sep = "")
