#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/influence/assess_core285_lodo.R
# Description:  Classify leave-one-donor-out support for the 285 core genes (strict >=90% and all-refit tiers).
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
# Run from the repository root, after scripts/influence/leave_one_donor_out.R:
#   Rscript scripts/influence/assess_core285_lodo.R 2>&1 | tee logs/assess_core285_lodo.log
#
# Rules (both methods, DESeq2 and edgeR):
#   complete_refits        all 34 refits gave a finite log2FC and an FDR
#   stable_direction       same sign as the full model in every refit
#   high_FDR_retention     FDR < 0.05 in at least 90% of refits (31 of 34)
#   strict_LODO_supported  complete_refits & stable_direction & high_FDR_retention
#   all_refits_significant complete_refits & stable_direction & FDR < 0.05 in every refit
#
# Writes core285_lodo_assessment.tsv (all 285 genes) and core285_lodo_review.tsv
# (genes not strictly supported). If earlier versions exist, they are compared first
# and kept as *.previous.tsv.

core_path    <- "metadata/GSE221156_expression/core_supported_genes.tsv"
summary_path <- "results/influence/GSE221156/gene_influence_summary.tsv"
out_dir      <- "results/influence/GSE221156"
retention_threshold <- 0.90

for (p in c(core_path, summary_path)) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)
core <- read.delim(core_path, check.names = FALSE)
lodo <- read.delim(summary_path, check.names = FALSE)
if (!"gene_id" %in% names(core)) stop(core_path, " lacks gene_id", call. = FALSE)

n_refits <- max(lodo$valid_FDR_DESeq2, na.rm = TRUE)
message("Refits per gene: ", n_refits)

# Keep core columns; add LODO columns that the core table does not already have
lodo_cols <- c("gene_id", setdiff(names(lodo), names(core)))
a <- merge(core, lodo[, lodo_cols], by = "gene_id", all.x = TRUE, sort = FALSE)
if (any(is.na(a$direction_consistency_DESeq2)))
  stop("Some core genes are missing from the LODO summary", call. = FALSE)

a$complete_refits <- with(a, valid_logFC_DESeq2 == n_refits & valid_logFC_edgeR == n_refits &
                             valid_FDR_DESeq2 == n_refits & valid_FDR_edgeR == n_refits)
a$stable_direction <- with(a, direction_consistency_DESeq2 == 1 & direction_consistency_edgeR == 1)
a$high_FDR_retention <- with(a, fraction_FDR05_DESeq2 >= retention_threshold &
                                fraction_FDR05_edgeR >= retention_threshold)
a$strict_LODO_supported <- with(a, complete_refits & stable_direction & high_FDR_retention)
a$all_refits_significant <- with(a, complete_refits & stable_direction &
                                    fraction_FDR05_DESeq2 == 1 & fraction_FDR05_edgeR == 1)

cat("Core genes:", nrow(a), "\n")
cat("Same direction in every refit:", sum(a$complete_refits & a$stable_direction), "\n")
cat("Strictly LODO-supported (FDR < 0.05 in >=", retention_threshold * 100, "% of refits):",
    sum(a$strict_LODO_supported), "\n")
cat("Significant in all", n_refits, "refits:", sum(a$all_refits_significant), "\n")

sym <- if ("gene_symbol" %in% names(a)) a$gene_symbol else a$gene_id
cat("Not strictly supported:", paste(sym[!a$strict_LODO_supported], collapse = ", "), "\n")
cat("Strict but not significant in every refit:",
    paste(sym[a$strict_LODO_supported & !a$all_refits_significant], collapse = ", "), "\n")

# Compare with any earlier (interactively created) assessment before overwriting
old_path <- file.path(out_dir, "core285_lodo_assessment.tsv")
if (file.exists(old_path)) {
  old <- read.delim(old_path, check.names = FALSE)
  if ("strict_LODO_supported" %in% names(old)) {
    same <- setequal(old$gene_id[as.logical(old$strict_LODO_supported)],
                     a$gene_id[a$strict_LODO_supported])
    cat("Matches earlier strict_LODO_supported gene set:", same, "\n")
    if (!same) warning("Strict gene set differs from the earlier file; review before using.")
  }
  file.copy(old_path, sub("\\.tsv$", ".previous.tsv", old_path), overwrite = TRUE)
}
review_path <- file.path(out_dir, "core285_lodo_review.tsv")
if (file.exists(review_path))
  file.copy(review_path, sub("\\.tsv$", ".previous.tsv", review_path), overwrite = TRUE)

write.table(a, old_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(a[!a$strict_LODO_supported, ], review_path, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_assess_core285_lodo.txt"))
message("Wrote ", old_path, " and ", review_path)
