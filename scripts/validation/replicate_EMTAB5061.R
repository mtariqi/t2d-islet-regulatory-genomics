#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/validation/replicate_EMTAB5061.R
# Description:  Pre-registered replication in E-MTAB-5061: beta-cell pseudobulk and frozen gene-set scores.
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
# Run from the repository root, after the count file has fully downloaded:
#   Rscript scripts/validation/replicate_EMTAB5061.R 2>&1 | tee logs/replicate_EMTAB5061.log
#
# Rules: docs/prereg_validation_v2.md, addendum "E-MTAB-5061" (commit a61a07f).
#   cells:   single cell quality == "OK" and inferred cell type == "beta cell"
#   link:    SDRF "Source Name" to count-file column names (by name, not position)
#   counts:  only the read-count block of pancreas_refseq_rpkms_counts_3514sc.txt
#            (row = symbol, refseq, n RPKM columns, n count columns); RPKMs never used
#   donors:  >= 10 beta cells (primary; all 10 donors). >= 20 set: descriptive only.
#   model:   primary score ~ disease + sex + age_z; sensitivity unadjusted score ~ disease
# Output columns match results/validation/GSE153855/replication_scores.tsv.

suppressPackageStartupMessages(library(edgeR))
set.seed(20261119)

d_dir   <- "data/raw/validation/E-MTAB-5061"
sdrf_f  <- file.path(d_dir, "E-MTAB-5061.sdrf.txt")
mat_f   <- file.path(d_dir, "pancreas_refseq_rpkms_counts_3514sc.txt")
out_dir <- "results/validation/E-MTAB-5061"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c(sdrf_f, mat_f)) if (!file.exists(f)) stop("Missing file: ", f, call. = FALSE)
expected_bytes <- 722957491
if (file.size(mat_f) != expected_bytes)
  stop("Count file is ", file.size(mat_f), " bytes, expected ", expected_bytes,
       ": download incomplete", call. = FALSE)

# ---------------- Frozen gene sets ----------------
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

# ---------------- Annotation ----------------
s <- read.delim(sdrf_f, check.names = FALSE, stringsAsFactors = FALSE)
s <- s[!duplicated(s[["Source Name"]]), ]
stopifnot("Source Name must be unique" = !anyDuplicated(s[["Source Name"]]))
cell <- data.frame(id = s[["Source Name"]],
                   donor = s[["Characteristics [individual]"]],
                   disease = s[["Characteristics [disease]"]],
                   sex = s[["Characteristics [sex]"]],
                   age = as.numeric(s[["Characteristics [age]"]]),
                   quality = s[["Characteristics [single cell quality]"]],
                   type = s[["Characteristics [inferred cell type]"]])
cell$disease <- ifelse(cell$disease == "normal", "ND",
                ifelse(cell$disease == "type II diabetes mellitus", "T2D", NA))
stopifnot(!anyNA(cell$disease))
beta_cells <- cell[cell$quality == "OK" & cell$type == "beta cell", ]
cat("Beta cells (quality OK):", nrow(beta_cells), "\n")
stopifnot("expected 270 beta cells, as in the addendum" = nrow(beta_cells) == 270)

# ---------------- Count block of the combined file ----------------
header <- strsplit(readLines(mat_f, n = 1), "\t", fixed = TRUE)[[1]]
names_hdr <- header[-1]
n <- length(names_hdr)
row1 <- strsplit(readLines(mat_f, n = 2)[2], "\t", fixed = TRUE)[[1]]
cat("Cells in header:", n, "| fields per data row:", length(row1), "\n")
stopifnot("data rows must be symbol + refseq + n RPKM + n counts" = length(row1) == 2 + 2 * n)
pos <- match(beta_cells$id, names_hdr)
stopifnot("every beta cell must be found in the count file header" = !anyNA(pos))
cls <- rep("NULL", 2 + 2 * n); cls[1] <- "character"; cls[2 + n + pos] <- "numeric"
mat <- read.delim(mat_f, header = FALSE, skip = 1, colClasses = cls, quote = "", comment.char = "")
genes <- mat[[1]]; mat <- as.matrix(mat[, -1])
colnames(mat) <- names_hdr[sort(pos)]
mat <- mat[, beta_cells$id]
stopifnot("counts must be non-negative" = all(mat >= 0, na.rm = TRUE), !anyNA(mat))
frac <- mean(abs(mat - round(mat)) > 1e-6)
cat("Fraction of non-integer count values:", signif(frac, 3), "\n")
if (frac > 0.001) stop("Count block does not look like read counts; check the column layout", call. = FALSE)
mat <- round(mat)
keep_rows <- !grepl("^ERCC-", genes)
cat("ERCC spike-in rows removed:", sum(!keep_rows), "\n")
mat <- mat[keep_rows, ]; genes <- genes[keep_rows]

# Pseudobulk by donor; duplicated symbols summed
pb <- t(rowsum(t(mat), group = beta_cells$donor))
pb <- rowsum(pb, group = genes)
stopifnot("count conservation" = abs(sum(pb) - sum(mat)) < 1e-6)
dinfo <- unique(beta_cells[, c("donor", "disease", "sex", "age")])
donors <- data.frame(donor = colnames(pb))
donors <- merge(donors, dinfo, by = "donor", sort = FALSE)
donors <- donors[match(colnames(pb), donors$donor), ]
donors$n_beta <- as.integer(table(beta_cells$donor)[donors$donor])
donors$total_counts <- colSums(pb)
write.table(donors, file.path(out_dir, "donors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDonors:\n"); print(donors, row.names = FALSE)

# ---------------- Analysis ----------------
ds <- donors$donor[donors$n_beta >= 10]
m <- donors[match(ds, donors$donor), ]
m$dis <- factor(m$disease, levels = c("ND", "T2D"))
m$sex <- factor(m$sex)
m$age_z <- as.numeric(scale(m$age))
cat("\nPrimary donor set:", sum(m$dis == "ND"), "ND vs", sum(m$dis == "T2D"), "T2D\n")

y <- DGEList(pb[, ds])
keep <- filterByExpr(y, design = model.matrix(~ dis, data = m))
y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
lcpm <- cpm(y, log = TRUE, prior.count = 2)
z <- t(scale(t(lcpm))); z <- z[is.finite(rowSums(z)), , drop = FALSE]
cat("Genes tested:", nrow(y), "\n")

up <- intersect(core$gene_symbol[core$expected_direction == "up"], rownames(z))
dn <- intersect(core$gene_symbol[core$expected_direction == "down"], rownames(z))
pg <- intersect(prot$gene_symbol, rownames(z))
kg <- intersect(kch$gene_symbol, rownames(z))
cat("Genes used: core up", length(up), "/", sum(core$expected_direction == "up"),
    "| core down", length(dn), "/", sum(core$expected_direction == "down"),
    "| proteasome", length(pg), "/ 29 | K-channel", length(kg), "/ 5\n")
m$core279_direction_aware <- colMeans(z[up, , drop = FALSE]) - colMeans(z[dn, , drop = FALSE])
m$kchannel5 <- colMeans(z[kg, , drop = FALSE])
m$proteasome29 <- colMeans(z[pg, , drop = FALSE])
write.table(m, file.path(out_dir, "donor_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

expect <- c(core279_direction_aware = "higher", kchannel5 = "higher", proteasome29 = "lower")
n_used <- c(core279_direction_aware = length(up) + length(dn), kchannel5 = length(kg), proteasome29 = length(pg))
fit_score <- function(s, adjusted) {
  f <- if (adjusted) lm(m[[s]] ~ dis + sex + age_z, data = m) else lm(m[[s]] ~ dis, data = m)
  cf <- summary(f)$coefficients; ci <- confint(f); r <- "disT2D"
  data.frame(donor_set = if (adjusted) "primary_min10" else "sens_unadjusted",
             n_ND = sum(m$dis == "ND"), n_T2D = sum(m$dis == "T2D"), score = s, genes_used = n_used[[s]],
             estimate_T2D_minus_ND = cf[r, 1], ci_low = ci[r, 1], ci_high = ci[r, 2],
             p_two_sided = cf[r, 4], expected = expect[[s]],
             direction_as_expected = sign(cf[r, 1]) == ifelse(expect[[s]] == "higher", 1, -1))
}
scores_out <- do.call(rbind, c(lapply(names(expect), fit_score, adjusted = TRUE),
                               lapply(names(expect), fit_score, adjusted = FALSE)))
write.table(scores_out, file.path(out_dir, "replication_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# fry and per-gene agreement, adjusted design
design <- model.matrix(~ dis + sex + age_z, data = m)
y <- estimateDisp(y, design, robust = TRUE)
fr <- fry(y, index = list(proteasome29 = pg, kchannel5 = kg, core279_up = up, core279_down = dn),
          design = design, contrast = 2)
fr$gene_set <- rownames(fr); fr$donor_set <- "primary_min10"
write.table(fr, file.path(out_dir, "replication_fry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
fit <- glmQLFit(y, design, robust = TRUE)
tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
cg <- core[core$gene_symbol %in% rownames(tt), ]
agree <- sign(tt[cg$gene_symbol, "logFC"]) == ifelse(cg$expected_direction == "up", 1, -1)
write.table(data.frame(gene_symbol = cg$gene_symbol, expected = cg$expected_direction,
                       discovery_log2FC = cg$discovery_log2FC, replication_logFC = tt[cg$gene_symbol, "logFC"],
                       replication_P = tt[cg$gene_symbol, "PValue"], agrees = agree),
            file.path(out_dir, "replication_per_gene_core279.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# >= 20 beta cells: descriptive only (4 ND vs 1 T2D), as pre-registered
d20 <- m[m$n_beta >= 20, ]
desc20 <- aggregate(d20[, names(expect)], by = list(disease = d20$dis), FUN = mean)
write.table(desc20, file.path(out_dir, "descriptive_min20.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nCore genes with expected direction:", sum(agree), "/", length(agree),
    sprintf("(%.1f%%)", 100 * mean(agree)), "| Spearman discovery vs replication logFC:",
    round(cor(cg$discovery_log2FC, tt[cg$gene_symbol, "logFC"], method = "spearman"), 3), "\n")
cat("\n=== SCORE TESTS (T2D minus ND, z units) ===\n")
print(scores_out[, c("donor_set", "score", "genes_used", "estimate_T2D_minus_ND", "ci_low", "ci_high",
                     "p_two_sided", "direction_as_expected")], row.names = FALSE, digits = 3)
cat("\n=== fry (adjusted) ===\n"); print(fr[, c("gene_set", "NGenes", "Direction", "PValue")], row.names = FALSE, digits = 3)
cat("\n=== >= 20 beta cells, descriptive group means (", sum(d20$dis == "ND"), "ND,",
    sum(d20$dis == "T2D"), "T2D) ===\n"); print(desc20, row.names = FALSE, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nSingle-cohort results. Pooling depends on the independence check in the addendum.\n")
