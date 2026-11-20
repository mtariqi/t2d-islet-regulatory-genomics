#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/validation/replicate_GSE153855.R
# Description:  Pre-registered replication in GSE153855: beta-cell pseudobulk and frozen gene-set scores.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================
#
# Run from the repository root:
#   Rscript scripts/validation/replicate_GSE153855.R 2>&1 | tee logs/replicate_GSE153855.log
#
# Follows docs/Independent_Validation_Plan.md (tag prereg-validation-v1) and
# docs/prereg_validation_v2.md with its addendum (tag prereg-validation-v2).
#
# Frozen gene sets (checked against their recorded SHA256):
#   proteasome29     expected LOWER in T2D  (v1)
#   core279_strict   direction-aware score expected HIGHER in T2D  (v2 primary)
#   kchannel5        expected HIGHER in T2D  (v2 secondary)
# Donor sets (fixed before expression was examined):
#   primary >= 10 beta cells; sensitivity >= 20; sensitivity all donors.
# No age/sex/BMI in GEO for this cohort, so models are unadjusted: score ~ disease.
# Scores: TMM log2 CPM of beta-cell pseudobulk, z-scored per gene across the donors
# in each donor set. Genes are matched by symbol; undetected genes are dropped, not imputed.

suppressPackageStartupMessages(library(edgeR))
set.seed(20261119)

d_dir   <- "data/raw/validation/GSE153855"
ann_f   <- file.path(d_dir, "GSE153855_Cell_annotation.txt.gz")
mat_f   <- file.path(d_dir, "GSE153855_Expression_counts_HQ_allsamples.txt.gz")
out_dir <- "results/validation/GSE153855"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- Frozen gene sets and checksums ----------------
check_hash <- function(file, hash_file) {
  rec <- read.table(hash_file, col.names = c("hash", "path"), stringsAsFactors = FALSE)
  want <- rec$hash[basename(rec$path) == basename(file)]
  got <- strsplit(system2("sha256sum", file, stdout = TRUE), " +")[[1]][1]
  if (length(want) != 1 || got != want) stop("Checksum mismatch for ", file, call. = FALSE)
  cat("Checksum OK:", file, "\n")
}
check_hash("config/gene_sets/proteasome29.tsv", "config/gene_sets/proteasome29.sha256")
check_hash("config/gene_sets/core279_strict.tsv", "config/gene_sets/prereg_v2.sha256")
check_hash("config/gene_sets/kchannel5.tsv", "config/gene_sets/prereg_v2.sha256")
prot <- read.delim("config/gene_sets/proteasome29.tsv")
core <- read.delim("config/gene_sets/core279_strict.tsv")
kch  <- read.delim("config/gene_sets/kchannel5.tsv")

# ---------------- Annotation (positional link to matrix columns) ----------------
# The source file ends with one line using CR separators; strip CR before parsing
# so that line is read as one row (otherwise R splits it into three).
strip_cr <- function(f) pipe(paste("zcat", shQuote(f), "| tr -d '\\r'"))
ann <- read.delim(strip_cr(ann_f), stringsAsFactors = FALSE, strip.white = TRUE)
ann[] <- lapply(ann, trimws)
stopifnot(identical(names(ann), c("Donor", "Disease", "CellType")))
first <- readLines(strip_cr(mat_f), n = 1)
n_fields <- length(strsplit(first, "\t", fixed = TRUE)[[1]])
cat("Annotated cells:", nrow(ann), "| matrix cell columns:", n_fields - 1, "\n")
stopifnot("matrix columns must equal annotation rows" = n_fields - 1 == nrow(ann))
ann$disease <- c(normal = "ND", `type II diabetes` = "T2D")[ann$Disease]
stopifnot(!anyNA(ann$disease))

beta <- which(ann$CellType == "Beta")
cat("Beta cells:", length(beta), "\n")

# Read only gene names and beta-cell columns
cls <- rep("NULL", n_fields); cls[1] <- "character"; cls[beta + 1] <- "integer"
mat <- read.delim(strip_cr(mat_f), header = FALSE, colClasses = cls, quote = "")
genes <- mat[[1]]; mat <- as.matrix(mat[, -1]); storage.mode(mat) <- "double"  # avoid integer overflow in sums
stopifnot("counts must be non-negative integers" = all(mat >= 0))

# Pseudobulk: sum beta-cell counts per donor; duplicated symbols summed
donor_of_cell <- ann$Donor[beta]
pb <- t(rowsum(t(mat), group = donor_of_cell))
pb <- rowsum(pb, group = genes)
stopifnot("count conservation" = sum(pb) == sum(mat))
donors <- data.frame(donor = colnames(pb),
                     disease = ann$disease[match(colnames(pb), ann$Donor)],
                     n_beta = as.integer(table(donor_of_cell)[colnames(pb)]),
                     total_counts = colSums(pb))
write.table(donors, file.path(out_dir, "donors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDonors:\n"); print(donors, row.names = FALSE)

# ---------------- Analysis per donor set ----------------
donor_sets <- list(primary_min10 = donors$donor[donors$n_beta >= 10],
                   sens_min20 = donors$donor[donors$n_beta >= 20],
                   sens_all = donors$donor)

score_test <- function(score, dis, expected) {
  f <- lm(score ~ dis); cf <- summary(f)$coefficients; ci <- confint(f)
  data.frame(estimate_T2D_minus_ND = cf[2, 1], ci_low = ci[2, 1], ci_high = ci[2, 2],
             p_two_sided = cf[2, 4], expected = expected,
             direction_as_expected = sign(cf[2, 1]) == ifelse(expected == "higher", 1, -1))
}

res <- list(); per_gene <- list()
for (set in names(donor_sets)) {
  ds <- donor_sets[[set]]
  dis <- factor(donors$disease[match(ds, donors$donor)], levels = c("ND", "T2D"))
  design <- model.matrix(~ dis)
  y <- DGEList(pb[, ds])
  keep <- filterByExpr(y, design = design)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  lcpm <- cpm(y, log = TRUE, prior.count = 2)
  z <- t(scale(t(lcpm))); z <- z[is.finite(rowSums(z)), , drop = FALSE]
  cat("\n=== ", set, ": ", sum(dis == "ND"), " ND vs ", sum(dis == "T2D"), " T2D; genes tested ",
      nrow(y), " ===\n", sep = "")

  up <- intersect(core$gene_symbol[core$expected_direction == "up"], rownames(z))
  dn <- intersect(core$gene_symbol[core$expected_direction == "down"], rownames(z))
  pg <- intersect(prot$gene_symbol, rownames(z))
  kg <- intersect(kch$gene_symbol, rownames(z))
  cat("Genes used: core up", length(up), "/", sum(core$expected_direction == "up"),
      "| core down", length(dn), "/", sum(core$expected_direction == "down"),
      "| proteasome", length(pg), "/ 29 | K-channel", length(kg), "/ 5\n")

  scores <- list(
    core279_direction_aware = colMeans(z[up, , drop = FALSE]) - colMeans(z[dn, , drop = FALSE]),
    kchannel5 = colMeans(z[kg, , drop = FALSE]),
    proteasome29 = colMeans(z[pg, , drop = FALSE]))
  expect <- c(core279_direction_aware = "higher", kchannel5 = "higher", proteasome29 = "lower")
  n_used <- c(core279_direction_aware = length(up) + length(dn), kchannel5 = length(kg), proteasome29 = length(pg))
  for (s in names(scores)) {
    r <- score_test(scores[[s]], dis, expect[[s]])
    res[[paste(set, s)]] <- cbind(donor_set = set, n_ND = sum(dis == "ND"), n_T2D = sum(dis == "T2D"),
                                  score = s, genes_used = n_used[[s]], r)
  }

  # fry gene-set tests (v1: proteasome module; descriptive for the others)
  y <- estimateDisp(y, design, robust = TRUE)
  fr <- fry(y, index = list(proteasome29 = pg, kchannel5 = kg, core279_up = up, core279_down = dn),
            design = design, contrast = 2)
  fr$gene_set <- rownames(fr); fr$donor_set <- set
  res[[paste(set, "fry")]] <- fr

  # Per-gene direction agreement (descriptive)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  cg <- core[core$gene_symbol %in% rownames(tt), ]
  agree <- sign(tt[cg$gene_symbol, "logFC"]) == ifelse(cg$expected_direction == "up", 1, -1)
  per_gene[[set]] <- data.frame(donor_set = set, gene_symbol = cg$gene_symbol,
                                expected = cg$expected_direction, discovery_log2FC = cg$discovery_log2FC,
                                replication_logFC = tt[cg$gene_symbol, "logFC"],
                                replication_P = tt[cg$gene_symbol, "PValue"], agrees = agree)
  cat("Core genes with expected direction:", sum(agree), "/", length(agree),
      sprintf("(%.1f%%)", 100 * mean(agree)),
      "| Spearman discovery vs replication logFC:",
      round(cor(cg$discovery_log2FC, tt[cg$gene_symbol, "logFC"], method = "spearman"), 3), "\n")
}

scores_out <- do.call(rbind, res[!grepl("fry$", names(res))])
fry_out <- do.call(rbind, res[grepl("fry$", names(res))])
write.table(scores_out, file.path(out_dir, "replication_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fry_out, file.path(out_dir, "replication_fry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, per_gene), file.path(out_dir, "replication_per_gene_core279.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== SCORE TESTS (T2D minus ND, z units) ===\n")
print(scores_out[, c("donor_set", "score", "genes_used", "estimate_T2D_minus_ND", "ci_low", "ci_high",
                     "p_two_sided", "direction_as_expected")], row.names = FALSE, digits = 3)
cat("\n=== fry ===\n")
print(fry_out[, c("donor_set", "gene_set", "NGenes", "Direction", "PValue")], row.names = FALSE, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nSingle-cohort results only. Pooled meta-analysis across cohorts is the pre-registered endpoint.\n")
