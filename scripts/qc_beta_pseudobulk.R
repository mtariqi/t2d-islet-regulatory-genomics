suppressPackageStartupMessages(library(DESeq2))

out <- "results/qc/GSE221156"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(
  "metadata/GSE221156_beta_pseudobulk_samples.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts <- as.matrix(read.delim(
  gzfile("data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz"),
  row.names = 1,
  check.names = FALSE
))

stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(samples$donor_id),
  setequal(colnames(counts), samples$donor_id),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == floor(counts)),
  max(counts) <= .Machine$integer.max
)

samples <- samples[match(colnames(counts), samples$donor_id), ]
rownames(samples) <- samples$donor_id
stopifnot(
  identical(colnames(counts), rownames(samples)),
  all(colSums(counts) == samples$total_counts),
  sum(counts) == 2081765364
)

storage.mode(counts) <- "integer"

samples$disease_status <- factor(
  samples$disease_status, levels = c("ND", "PD", "T2D")
)
samples$chemistry <- factor(samples$chemistry)
samples$sex <- factor(samples$sex)
samples$age <- as.numeric(scale(as.numeric(samples$age)))

# Exploratory QC filter, applied across all 48 donors.
# Differential-testing filters will be documented separately.
keep <- rowSums(counts >= 10) >= 10
if (sum(keep) < 2) stop("Too few expressed genes for QC")

dds <- DESeqDataSetFromMatrix(
  countData = counts[keep, ],
  colData = samples,
  design = ~ chemistry + sex + age + disease_status
)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
transformed <- assay(vsd)

variances <- apply(transformed, 1, var)
top <- head(order(variances, decreasing = TRUE), 500)
pca <- prcomp(t(transformed[top, ]), center = TRUE, scale. = FALSE)
variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)

scores <- data.frame(
  donor_id = rownames(pca$x),
  disease_status = samples$disease_status,
  chemistry = samples$chemistry,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
write.table(
  scores, file.path(out, "PCA_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

distances <- as.matrix(dist(t(transformed)))
median_distances <- apply(
  distances, 1,
  function(x) median(x[x > 0])
)
distance_mad <- mad(median_distances)
review_flag <- if (distance_mad > 0) {
  median_distances >
    median(median_distances) + 3 * distance_mad
} else {
  rep(FALSE, nrow(samples))
}

qc <- data.frame(
  donor_id = samples$donor_id,
  disease_status = samples$disease_status,
  chemistry = samples$chemistry,
  n_beta_cells = samples$n_beta_cells,
  total_counts = colSums(counts),
  detected_genes = colSums(counts > 0),
  counts_per_beta_cell = colSums(counts) / samples$n_beta_cells,
  size_factor = sizeFactors(dds),
  median_VST_distance = median_distances,
  review_flag = review_flag
)
write.table(
  qc, file.path(out, "donor_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

colors <- c(ND = "#2878B5", PD = "#E6A000", T2D = "#C83E4D")
point_colors <- colors[as.character(samples$disease_status)]
point_shapes <- ifelse(samples$chemistry == "V2", 17, 16)

pdf(file.path(out, "beta_pseudobulk_QC.pdf"), width = 11, height = 8.5)

par(mar = c(5, 5, 3, 2))
plot(
  pca$x[, 1], pca$x[, 2],
  col = point_colors, pch = point_shapes,
  xlab = sprintf("PC1 (%.1f%%)", variance_percent[1]),
  ylab = sprintf("PC2 (%.1f%%)", variance_percent[2]),
  main = "Beta-cell donor PCA - top 500 variable genes"
)
text(
  pca$x[, 1], pca$x[, 2],
  labels = samples$donor_id, pos = 3, cex = 0.65
)
legend(
  "topright", legend = names(colors),
  col = colors, pch = 16, bty = "n"
)
legend(
  "bottomright", legend = c("V2", "V3"),
  pch = c(17, 16), bty = "n"
)

par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
boxplot(
  total_counts ~ disease_status, data = qc,
  log = "y", ylab = "Total counts", xlab = "Disease"
)
boxplot(
  n_beta_cells ~ disease_status, data = qc,
  ylab = "Beta cells", xlab = "Disease"
)
plot(
  qc$n_beta_cells, qc$total_counts,
  log = "y", col = point_colors, pch = point_shapes,
  xlab = "Beta cells", ylab = "Total counts"
)
plot(
  qc$n_beta_cells, qc$counts_per_beta_cell,
  col = point_colors, pch = point_shapes,
  xlab = "Beta cells", ylab = "Counts per beta cell"
)

par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
plot(
  hclust(as.dist(distances)),
  main = "Donor clustering using VST expression",
  xlab = "", sub = "", cex = 0.7
)
dev.off()

saveRDS(vsd, file.path(out, "beta_QC_vst.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(out, "sessionInfo.txt")
)

cat("QC donors:", ncol(counts), "\n")
cat("Input genes:", nrow(counts), "\n")
cat("Genes retained for QC:", sum(keep), "\n")
cat("PC1 variance:", round(variance_percent[1], 2), "%\n")
cat("PC2 variance:", round(variance_percent[2], 2), "%\n")
cat("Donors flagged for review:",
    paste(qc$donor_id[qc$review_flag], collapse = ", "), "\n")
cat("No donors were excluded.\n")
cat("Saved:", file.path(out, "beta_pseudobulk_QC.pdf"), "\n")
