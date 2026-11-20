#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/robustness/disease_chemistry_interaction.R
# Description:  Test whether the T2D effect differs between V2 and V3 chemistry (core genes and proteasome module).
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
#   Rscript scripts/robustness/disease_chemistry_interaction.R 2>&1 | tee logs/disease_chemistry_interaction.log
#
# Model (17 ND vs 17 T2D donors, same as the primary analysis plus one term):
#   ~ chemistry + sex + age_z + disease_status + chemistry:disease_status
# Reference levels: chemistry V3, disease ND. Positive log2FC = higher in T2D.
# Genes: filterByExpr on the PRIMARY design, so the gene set matches the primary analysis.
#
# Outputs (results/robustness/GSE221156/chemistry_interaction/):
#   gene_interaction.tsv        all tested genes, both methods
#   core285_interaction.tsv     the 285 core genes, interaction FDR within the core set
#   module_interaction.tsv      proteasome module: fry gene-set test and donor-level module score
#
# Power note: V2 has only 11 donors (4 ND, 7 T2D). A non-significant interaction is
# not evidence of no interaction; read the effect estimates, not only the p-values.

suppressPackageStartupMessages({
  library(edgeR)
  library(DESeq2)
})

counts_path  <- "data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz"
samples_path <- "metadata/GSE221156_beta_pseudobulk_samples.tsv"
core_path    <- "metadata/GSE221156_expression/core_supported_genes.tsv"
module_path  <- "config/gene_sets/proteasome29.tsv"
primary_path <- "results/expression/GSE221156/primary_DESeq2_edgeR.tsv"
out_dir      <- "results/robustness/GSE221156/chemistry_interaction"
set.seed(20261119)

for (p in c(counts_path, samples_path, core_path, module_path, primary_path))
  if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- Data ----------------
counts <- read.delim(counts_path, check.names = FALSE)
rownames(counts) <- counts[[1]]
counts <- as.matrix(counts[, -1, drop = FALSE])
m <- read.delim(samples_path, check.names = FALSE)
m <- m[m$disease_status %in% c("ND", "T2D"), ]
stopifnot("expected 34 ND/T2D donors" = nrow(m) == 34L)
stopifnot("donors missing from count matrix" = all(m$donor_id %in% colnames(counts)))
counts <- counts[, m$donor_id]
m$disease_status <- factor(m$disease_status, levels = c("ND", "T2D"))
m$chemistry <- factor(m$chemistry, levels = c("V3", "V2"))
m$sex <- factor(m$sex)
m$age_z <- as.numeric(scale(as.numeric(m$age)))
rownames(m) <- m$donor_id
cat("Donors by chemistry and disease:\n"); print(table(m$chemistry, m$disease_status))

core_ids <- read.delim(core_path)$gene_id
module_ids <- read.delim(module_path)$gene_id
primary <- read.delim(primary_path)
sym <- setNames(primary$gene_symbol, primary$gene_id)

# ---------------- Designs and filtering ----------------
primary_design <- model.matrix(~ chemistry + sex + age_z + disease_status, data = m)
int_design <- model.matrix(~ chemistry + sex + age_z + disease_status + chemistry:disease_status, data = m)
int_coef <- "chemistryV2:disease_statusT2D"
dis_coef <- grep("^disease_status", colnames(int_design), value = TRUE)
dis_coef <- dis_coef[!grepl(":", dis_coef)]
stopifnot(int_coef %in% colnames(int_design), length(dis_coef) == 1)
stopifnot("interaction design not full rank" = qr(int_design)$rank == ncol(int_design))
cat("Residual df (interaction model):", nrow(int_design) - ncol(int_design), "\n")

y <- DGEList(counts)
keep <- filterByExpr(y, design = primary_design)
y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
cat("Genes tested:", nrow(y), "(primary analysis reported 17,352)\n")
missing_core <- setdiff(core_ids, rownames(y))
if (length(missing_core)) cat("Core genes not retained by filter:", length(missing_core), "\n")

# ---------------- edgeR robust QL ----------------
y <- estimateDisp(y, int_design, robust = TRUE)
fit <- glmQLFit(y, int_design, robust = TRUE)
ci <- which(colnames(int_design) == int_coef)
di <- which(colnames(int_design) == dis_coef)
e_int <- topTags(glmQLFTest(fit, coef = ci), n = Inf, sort.by = "none")$table
e_v3 <- topTags(glmQLFTest(fit, coef = di), n = Inf, sort.by = "none")$table
con_v2 <- numeric(ncol(int_design)); con_v2[c(di, ci)] <- 1
e_v2 <- topTags(glmQLFTest(fit, contrast = con_v2), n = Inf, sort.by = "none")$table

# ---------------- DESeq2 Wald ----------------
dds <- DESeqDataSetFromMatrix(round(y$counts), colData = m,
                              design = ~ chemistry + sex + age_z + disease_status + chemistry:disease_status)
dds <- DESeq(dds, quiet = TRUE)
rn <- resultsNames(dds)
d_int_name <- grep("chemistryV2.disease_statusT2D", rn, value = TRUE)
d_dis_name <- grep("^disease_status_T2D_vs_ND$", rn, value = TRUE)
stopifnot(length(d_int_name) == 1, length(d_dis_name) == 1)
d_int <- results(dds, name = d_int_name, independentFiltering = FALSE)
d_v3 <- results(dds, name = d_dis_name, independentFiltering = FALSE)
d_v2 <- results(dds, contrast = list(c(d_dis_name, d_int_name)), independentFiltering = FALSE)

g <- rownames(y)
genes <- data.frame(
  gene_id = g, gene_symbol = unname(sym[g]),
  is_core = g %in% core_ids, is_module = g %in% module_ids,
  edgeR_T2D_effect_V3 = e_v3$logFC, edgeR_T2D_effect_V2 = e_v2$logFC,
  edgeR_interaction_logFC = e_int$logFC, edgeR_interaction_P = e_int$PValue, edgeR_interaction_FDR = e_int$FDR,
  DESeq2_T2D_effect_V3 = d_v3[g, "log2FoldChange"], DESeq2_T2D_effect_V2 = d_v2[g, "log2FoldChange"],
  DESeq2_interaction_log2FC = d_int[g, "log2FoldChange"], DESeq2_interaction_P = d_int[g, "pvalue"],
  DESeq2_interaction_FDR = d_int[g, "padj"]
)
write.table(genes, file.path(out_dir, "gene_interaction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------- Core genes ----------------
core <- genes[genes$is_core, ]
core$edgeR_interaction_FDR_core <- p.adjust(core$edgeR_interaction_P, "BH")
core$DESeq2_interaction_FDR_core <- p.adjust(core$DESeq2_interaction_P, "BH")
core$same_sign_V2_V3 <- sign(core$edgeR_T2D_effect_V2) == sign(core$edgeR_T2D_effect_V3) &
                        sign(core$DESeq2_T2D_effect_V2) == sign(core$DESeq2_T2D_effect_V3)
core$interaction_both_methods <- core$edgeR_interaction_FDR_core < 0.05 & core$DESeq2_interaction_FDR_core < 0.05
write.table(core, file.path(out_dir, "core285_interaction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== CORE GENES (n =", nrow(core), ") ===\n")
cat("Same direction of T2D effect in V2 and V3 (both methods):", sum(core$same_sign_V2_V3, na.rm = TRUE), "\n")
cat("Interaction FDR < 0.05 within core set, edgeR:", sum(core$edgeR_interaction_FDR_core < 0.05, na.rm = TRUE),
    "| DESeq2:", sum(core$DESeq2_interaction_FDR_core < 0.05, na.rm = TRUE),
    "| both:", sum(core$interaction_both_methods, na.rm = TRUE), "\n")
cat("Genome-wide interaction FDR < 0.05, edgeR:", sum(genes$edgeR_interaction_FDR < 0.05, na.rm = TRUE),
    "| DESeq2:", sum(genes$DESeq2_interaction_FDR < 0.05, na.rm = TRUE), "\n")
flag <- core[core$interaction_both_methods %in% TRUE | !core$same_sign_V2_V3 %in% TRUE, ]
if (nrow(flag)) {
  cat("Core genes with an interaction or a sign difference:\n")
  print(flag[, c("gene_symbol", "edgeR_T2D_effect_V3", "edgeR_T2D_effect_V2", "edgeR_interaction_FDR_core")],
        row.names = FALSE, digits = 3)
}

# ---------------- Proteasome module ----------------
mod_in <- intersect(module_ids, rownames(y))
cat("\n=== PROTEASOME MODULE (", length(mod_in), "/29 genes tested) ===\n", sep = "")
fry_res <- fry(y, index = list(proteasome29 = mod_in), design = int_design, contrast = ci)
cat("fry interaction test: direction", as.character(fry_res$Direction),
    "| PValue", signif(fry_res$PValue, 3), "| PValue.Mixed", signif(fry_res$PValue.Mixed, 3), "\n")

lcpm <- cpm(y, log = TRUE, prior.count = 2)
z <- t(scale(t(lcpm[mod_in, , drop = FALSE])))
m$module_score <- colMeans(z)
lmfit <- lm(module_score ~ chemistry + sex + age_z + disease_status + chemistry:disease_status, data = m)
cf <- summary(lmfit)$coefficients
ci_lm <- confint(lmfit)
int_row <- grep(":", rownames(cf), value = TRUE)
v3_row <- grep("^disease_status", rownames(cf), value = TRUE); v3_row <- v3_row[!grepl(":", v3_row)]
L <- setNames(numeric(nrow(cf)), rownames(cf)); L[c(v3_row, int_row)] <- 1
v2_est <- sum(L * coef(lmfit)); v2_se <- sqrt(drop(t(L) %*% vcov(lmfit) %*% L))
v2_p <- 2 * pt(-abs(v2_est / v2_se), df = lmfit$df.residual)
tcrit <- qt(0.975, lmfit$df.residual)

module_out <- data.frame(
  quantity = c("T2D effect within V3", "T2D effect within V2", "Interaction (V2 minus V3)", "fry interaction gene-set test"),
  estimate = c(cf[v3_row, 1], v2_est, cf[int_row, 1], NA),
  ci_low = c(ci_lm[v3_row, 1], v2_est - tcrit * v2_se, ci_lm[int_row, 1], NA),
  ci_high = c(ci_lm[v3_row, 2], v2_est + tcrit * v2_se, ci_lm[int_row, 2], NA),
  p_value = c(cf[v3_row, 4], v2_p, cf[int_row, 4], fry_res$PValue),
  note = c("module score units (mean z)", "module score units (mean z)", "module score units (mean z)",
           paste("direction", fry_res$Direction, "; descriptive, genes selected from these donors"))
)
write.table(module_out, file.path(out_dir, "module_interaction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(module_out[, 1:5], row.names = FALSE, digits = 3)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nWrote results to", out_dir, "\n")
