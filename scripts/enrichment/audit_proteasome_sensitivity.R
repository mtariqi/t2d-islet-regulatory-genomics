base <- "results/expression/GSE221156"
enrichment <- "results/enrichment/GSE221156"

files <- c(
  primary = file.path(base, "primary_DESeq2_edgeR.tsv"),
  BMI = file.path(base, "BMI_sensitivity_DESeq2_edgeR.tsv"),
  V3 = file.path(base, "V3_only_DESeq2_edgeR.tsv")
)

shared <- read.delim(
  file.path(
    enrichment,
    "Reactome_negative18_shared29_gene_characterization.tsv"
  ),
  check.names = FALSE
)

genes <- unique(shared$gene_id)

stopifnot(length(genes) == 29L)

required <- c(
  "gene_id",
  "DESeq2_log2FC",
  "DESeq2_FDR_native",
  "edgeR_log2FC",
  "edgeR_FDR"
)

results <- list()

for (model in names(files)) {

  d <- read.delim(
    files[[model]],
    check.names = FALSE
  )

  stopifnot(all(required %in% names(d)))
  stopifnot(!anyDuplicated(d$gene_id))

  z <- d[match(genes, d$gene_id), ]

  z$model <- model

  results[[model]] <- z
}

combined <- do.call(rbind, results)

combined$negative_both <- with(
  combined,
  !is.na(DESeq2_log2FC) &
  !is.na(edgeR_log2FC) &
  DESeq2_log2FC < 0 &
  edgeR_log2FC < 0
)

combined$significant_both <- with(
  combined,
  !is.na(DESeq2_FDR_native) &
  !is.na(edgeR_FDR) &
  DESeq2_FDR_native < 0.05 &
  edgeR_FDR < 0.05
)

cat("\n=== PROTEASOME SENSITIVITY AUDIT ===\n")

for (model in names(files)) {

  z <- combined[combined$model == model, ]

  cat("\nModel:", model, "\n")

  cat(
    "Genes matched:",
    sum(!is.na(z$DESeq2_log2FC)),
    "/29\n"
  )

  cat(
    "Negative in both methods:",
    sum(z$negative_both),
    "/29\n"
  )

  cat(
    "Significant in both methods:",
    sum(z$significant_both),
    "/29\n"
  )

  cat(
    "Median DESeq2 log2FC:",
    median(z$DESeq2_log2FC, na.rm = TRUE),
    "\n"
  )

  cat(
    "Median edgeR log2FC:",
    median(z$edgeR_log2FC, na.rm = TRUE),
    "\n"
  )
}

# Gene-level stability across all three models
stability <- aggregate(
  negative_both ~ gene_id,
  data = combined,
  FUN = function(v) length(v) == 3L && all(v)
)

names(stability)[2] <- "negative_all_models"

cat("\n=== CROSS-MODEL STABILITY ===\n")

cat(
  "Genes negative in both methods across all models:",
  sum(stability$negative_all_models),
  "/29\n"
)

# Preserve original gene symbols
symbols <- shared[, c("gene_id", "gene_symbol")]

stability <- merge(
  stability,
  symbols,
  by = "gene_id",
  all.x = TRUE
)

cat("\nGenes not consistently negative:\n")

print(
  stability[!stability$negative_all_models, ],
  row.names = FALSE
)

write.table(
  combined,
  file.path(
    enrichment,
    "proteasome29_sensitivity_gene_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  stability,
  file.path(
    enrichment,
    "proteasome29_cross_model_stability.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSensitivity audit completed.\n")
