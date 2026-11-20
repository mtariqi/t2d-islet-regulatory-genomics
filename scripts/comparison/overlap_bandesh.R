#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/comparison/overlap_bandesh.R
# Description:  Overlap and direction agreement between the 285 core genes and Bandesh et al. T2D beta-cell DEGs.
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
# Usage (from the repository root):
#   Rscript scripts/comparison/overlap_bandesh.R <their_table.xlsx|.csv|.tsv> [sheet]
#
# The source table must contain the T2D vs ND beta-cell differential expression results
# from Bandesh et al. (EMBO J 2026; doi:10.1038/s44318-026-00744-w). Column names are
# detected automatically (gene symbol, log fold-change, FDR); the detected columns are
# printed so they can be checked. Their DEG set is defined as FDR < 0.05, as in the paper.
#
# Background for the enrichment test: genes tested in our primary analysis that also
# appear in their table (if their table lists all tested genes), otherwise our tested genes.
# The test is descriptive: the two analyses share donors, so overlap is expected and does
# not count as independent replication.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript scripts/comparison/overlap_bandesh.R <table> [sheet]", call. = FALSE)
f <- args[1]; sheet <- if (length(args) > 1) args[2] else 1
if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)

if (grepl("\\.xlsx?$", f, ignore.case = TRUE)) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Install readxl first: install.packages('readxl')", call. = FALSE)
  cat("Sheets in file:", paste(readxl::excel_sheets(f), collapse = " | "), "\n")
  their <- as.data.frame(readxl::read_excel(f, sheet = if (grepl("^[0-9]+$", sheet)) as.integer(sheet) else sheet))
} else {
  their <- read.delim(f, sep = if (grepl("\\.csv$", f, ignore.case = TRUE)) "," else "\t", check.names = FALSE)
}
cat("Their table:", nrow(their), "rows; columns:", paste(names(their), collapse = ", "), "\n")

pick <- function(patterns, what) {
  for (p in patterns) {
    hit <- grep(p, names(their), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  stop("Could not find a ", what, " column. Columns are: ", paste(names(their), collapse = ", "), call. = FALSE)
}
g_col  <- pick(c("^gene$", "symbol", "gene_name", "^gene", "^name$", "^\\.\\.\\.1$", "^X$"), "gene symbol")
fc_col <- pick(c("^logFC$", "log2FC", "log2.?fold", "logfc", "fold"), "log fold-change")
q_col  <- pick(c("^FDR$", "padj", "adj", "q.?val", "fdr"), "FDR")
cat("Using columns: gene =", g_col, "| logFC =", fc_col, "| FDR =", q_col, "\n")

their <- data.frame(gene = as.character(their[[g_col]]), logFC = as.numeric(their[[fc_col]]),
                    FDR = as.numeric(their[[q_col]]))
their <- their[!is.na(their$gene) & nzchar(their$gene) & is.finite(their$logFC), ]
their <- their[!duplicated(their$gene), ]
deg <- their[!is.na(their$FDR) & their$FDR < 0.05, ]
cat("Their DEGs at FDR < 0.05:", nrow(deg), "(paper reports 746)\n")
if (nrow(deg) < 100) warning("Fewer DEGs than expected: is this the full T2D vs ND beta-cell table?")

core <- read.delim("results/influence/GSE221156/core285_lodo_assessment.tsv")
core <- data.frame(gene = core$gene_symbol, our_log2FC = core$primary_DESeq2_log2FC,
                   strict = as.logical(core$strict_LODO_supported))
tested <- unique(read.delim("results/expression/GSE221156/primary_DESeq2_edgeR.tsv")$gene_symbol)

full_table <- nrow(their) > 2 * nrow(deg)
universe <- if (full_table) intersect(tested, their$gene) else tested
cat("Background:", length(universe), "genes",
    if (full_table) "(tested in both analyses)" else "(our tested genes; their table appears to list DEGs only)", "\n")

core_u <- core[core$gene %in% universe, ]
deg_u <- deg[deg$gene %in% universe, ]
ov <- merge(core_u, deg_u, by = "gene")
ov$same_direction <- sign(ov$our_log2FC) == sign(ov$logFC)

a <- nrow(ov); b <- nrow(core_u) - a; c <- nrow(deg_u) - a; d <- length(universe) - a - b - c
ft <- fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater")
expected <- nrow(core_u) * nrow(deg_u) / length(universe)

cat("\n=== OVERLAP ===\n")
cat("Core genes in background:", nrow(core_u), "of", nrow(core), "\n")
cat("Core genes that are Bandesh DEGs:", a, sprintf("(%.1f%%; %.1f expected by chance)", 100 * a / nrow(core_u), expected), "\n")
cat("Fisher's exact test (one-sided): odds ratio", signif(ft$estimate, 3), "| p =", signif(ft$p.value, 3), "\n")
cat("Same direction among overlapping genes:", sum(ov$same_direction), "of", a, "\n")
cat("Strictly LODO-supported core genes among the overlap:", sum(ov$strict), "of", a, "\n")
cat("Core genes not among their DEGs:", nrow(core_u) - a, "\n")
if (any(!ov$same_direction)) {
  cat("Opposite direction:", paste(ov$gene[!ov$same_direction], collapse = ", "), "\n")
}
cat("Potassium-channel genes in their DEG list:",
    paste(intersect(c("KCNA1", "KCNAB1", "KCNAB2", "KCNH8", "KCNJ2"), deg$gene), collapse = ", "), "\n")

dir.create("results/comparison", recursive = TRUE, showWarnings = FALSE)
core_all <- merge(core, their, by = "gene", all.x = TRUE)
core_all$bandesh_DEG <- core_all$gene %in% deg$gene
write.table(core_all, "results/comparison/core285_vs_bandesh.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(sprintf("source_table\t%s", f), sprintf("their_DEGs\t%d", nrow(deg)),
             sprintf("background\t%d", length(universe)), sprintf("core_in_background\t%d", nrow(core_u)),
             sprintf("overlap\t%d", a), sprintf("expected\t%.1f", expected),
             sprintf("odds_ratio\t%.3g", ft$estimate), sprintf("fisher_p\t%.3g", ft$p.value),
             sprintf("same_direction\t%d", sum(ov$same_direction))),
           "results/comparison/core285_vs_bandesh_summary.tsv")
cat("\nWrote results/comparison/core285_vs_bandesh.tsv and _summary.tsv\n")
