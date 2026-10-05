input <- "results/expression/GSE221156"
output <- file.path(input, "robustness")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

analyses <- c("primary", "BMI_sensitivity", "V3_only")

tables <- setNames(lapply(analyses, function(name) {
  path <- file.path(input, paste0(name, "_DESeq2_edgeR.tsv"))
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

  required <- c(
    "gene_id", "gene_symbol",
    "DESeq2_log2FC", "DESeq2_FDR_no_independent_filter",
    "edgeR_log2FC", "edgeR_FDR"
  )
  stopifnot(all(required %in% names(x)))
  stopifnot(!anyNA(x$gene_id), !anyDuplicated(x$gene_id))

  x$supported <- with(
    x,
    is.finite(DESeq2_FDR_no_independent_filter) &
    is.finite(edgeR_FDR) &
    is.finite(DESeq2_log2FC) &
    is.finite(edgeR_log2FC) &
    DESeq2_FDR_no_independent_filter < 0.05 &
    edgeR_FDR < 0.05 &
    sign(DESeq2_log2FC) == sign(edgeR_log2FC)
  )

  rownames(x) <- x$gene_id
  x
}), analyses)

# Keep all genes tested in any analysis; untested genes remain NA.
genes <- sort(unique(unlist(lapply(tables, function(x) x$gene_id))))
result <- data.frame(gene_id = genes, stringsAsFactors = FALSE)

for (name in analyses) {
  x <- tables[[name]]
  idx <- match(genes, x$gene_id)

  if (name == "primary") {
    result$gene_symbol <- x$gene_symbol[idx]
  }

  for (column in c(
    "DESeq2_log2FC", "DESeq2_FDR_no_independent_filter",
    "edgeR_log2FC", "edgeR_FDR", "supported"
  )) {
    result[[paste(name, column, sep = "_")]] <- x[[column]][idx]
  }
}

support_columns <- paste0(analyses, "_supported")
support <- as.matrix(result[, support_columns])
result$n_analyses_supported <- rowSums(support, na.rm = TRUE)

fc_columns <- unlist(lapply(analyses, function(name) {
  paste(name, c("DESeq2_log2FC", "edgeR_log2FC"), sep = "_")
}))
fc <- as.matrix(result[, fc_columns])

result$consistent_direction_all_six <- apply(fc, 1, function(values) {
  all(is.finite(values)) &&
    (all(values > 0) || all(values < 0))
})

result$core_supported <- (
  result$n_analyses_supported == 3 &
  result$consistent_direction_all_six
)

fdr_columns <- unlist(lapply(analyses, function(name) {
  paste(name, c("DESeq2_FDR_no_independent_filter", "edgeR_FDR"), sep = "_")
}))
fdr <- as.matrix(result[, fdr_columns])

result$maximum_FDR_all_six <- apply(fdr, 1, function(values) {
  if (all(is.finite(values))) max(values) else NA_real_
})

core <- result[result$core_supported, ]
core <- core[order(core$maximum_FDR_all_six, core$gene_id), ]

write.table(
  result, file.path(output, "all_genes_robustness.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  core, file.path(output, "core_supported_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary <- data.frame(
  analysis = analyses,
  genes_tested = vapply(tables, nrow, integer(1)),
  supported_both_methods = vapply(
    tables, function(x) sum(x$supported), integer(1)
  )
)
write.table(
  summary, file.path(output, "robustness_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output, "sessionInfo.txt")
)

print(summary, row.names = FALSE)
cat("\nCore genes supported across all six tests:", nrow(core), "\n")
cat("Saved outputs to:", output, "\n")
