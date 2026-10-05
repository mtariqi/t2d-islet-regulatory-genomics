suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))

out <- "results/expression/GSE221156"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(
  "metadata/GSE221156_beta_pseudobulk_samples.tsv",
  stringsAsFactors = FALSE, check.names = FALSE
)
counts <- as.matrix(read.delim(
  gzfile("data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz"),
  row.names = 1, check.names = FALSE
))
annotations <- read.delim(
  "metadata/GSE221156_beta_gene_annotations.tsv",
  stringsAsFactors = FALSE
)

stopifnot(
  !anyDuplicated(samples$donor_id),
  !anyDuplicated(rownames(counts)),
  setequal(colnames(counts), samples$donor_id),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == floor(counts)),
  max(counts) <= .Machine$integer.max
)
samples <- samples[match(colnames(counts), samples$donor_id), ]
rownames(samples) <- samples$donor_id
stopifnot(all(colSums(counts) == samples$total_counts))
storage.mode(counts) <- "integer"

configs <- list(
  primary = list(
    selected = samples$disease_status %in% c("ND", "T2D"),
    formula = ~ chemistry + sex + age_z + disease_status
  ),
  BMI_sensitivity = list(
    selected = samples$disease_status %in% c("ND", "T2D"),
    formula = ~ chemistry + sex + age_z + BMI_z + disease_status
  ),
  V3_only = list(
    selected = samples$disease_status %in% c("ND", "T2D") &
      samples$chemistry == "V3",
    formula = ~ sex + age_z + disease_status
  )
)

summaries <- list()
pdf(file.path(out, "DESeq2_edgeR_effect_comparison.pdf"),
    width = 8, height = 7)

for (name in names(configs)) {
  config <- configs[[name]]
  m <- samples[config$selected, , drop = FALSE]
  m$disease_status <- factor(m$disease_status, levels = c("ND", "T2D"))
  m$chemistry <- factor(m$chemistry)
  m$sex <- factor(m$sex)
  m$age_z <- as.numeric(scale(as.numeric(m$age)))
  m$BMI_z <- as.numeric(scale(as.numeric(m$BMI)))

  design <- model.matrix(config$formula, data = m)
  stopifnot(
    !anyNA(design),
    nrow(design) == nrow(m),
    qr(design)$rank == ncol(design)
  )

  cts <- counts[, rownames(m), drop = FALSE]

  # Same expression filter for both methods within this analysis.
  keep <- edgeR::filterByExpr(
    edgeR::DGEList(counts = cts), design = design
  )
  cts <- cts[keep, , drop = FALSE]
  if (!nrow(cts)) stop("No genes retained")

  write.table(
    data.frame(gene_id = rownames(cts)),
    file.path(out, paste0(name, "_tested_genes.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    m, file.path(out, paste0(name, "_samples.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # DESeq2: primary method; disable automatic count replacement.
  dds <- DESeqDataSetFromMatrix(
    countData = cts, colData = m, design = config$formula
  )
  dds <- DESeq(dds, minReplicatesForReplace = Inf)

  native <- results(
    dds, contrast = c("disease_status", "T2D", "ND"), alpha = 0.05
  )
  comparison <- results(
    dds,
    contrast = c("disease_status", "T2D", "ND"),
    alpha = 0.05,
    independentFiltering = FALSE
  )

  # edgeR: TMM normalization and robust quasi-likelihood model.
  y <- edgeR::DGEList(counts = cts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  coefficient <- match("disease_statusT2D", colnames(design))
  stopifnot(!is.na(coefficient))
  test <- edgeR::glmQLFTest(fit, coef = coefficient)
  er <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  er <- er[rownames(cts), , drop = FALSE]

  tab <- data.frame(
    gene_id = rownames(cts),
    gene_symbol = annotations$gene_symbol[
      match(rownames(cts), annotations$gene_id)
    ],
    DESeq2_baseMean = native$baseMean,
    DESeq2_log2FC = native$log2FoldChange,
    DESeq2_pvalue = native$pvalue,
    DESeq2_FDR_native = native$padj,
    DESeq2_FDR_no_independent_filter = comparison$padj,
    edgeR_log2FC = er$logFC,
    edgeR_logCPM = er$logCPM,
    edgeR_pvalue = er$PValue,
    edgeR_FDR = er$FDR
  )

  valid <- is.finite(tab$DESeq2_pvalue) &
    is.finite(tab$edgeR_pvalue) &
    is.finite(tab$DESeq2_log2FC) &
    is.finite(tab$edgeR_log2FC)

  ds_sig <- !is.na(tab$DESeq2_FDR_no_independent_filter) &
    tab$DESeq2_FDR_no_independent_filter < 0.05
  er_sig <- !is.na(tab$edgeR_FDR) & tab$edgeR_FDR < 0.05
  both <- ds_sig & er_sig
  either <- ds_sig | er_sig

  tab$significant_both <- both
  tab$same_effect_direction <-
    sign(tab$DESeq2_log2FC) == sign(tab$edgeR_log2FC)

  write.table(
    tab, file.path(out, paste0(name, "_DESeq2_edgeR.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  summaries[[name]] <- data.frame(
    analysis = name,
    n_donors = nrow(m),
    genes_tested = nrow(tab),
    genes_with_finite_results_both = sum(valid),
    DESeq2_native_FDR05 = sum(
      !is.na(tab$DESeq2_FDR_native) & tab$DESeq2_FDR_native < 0.05
    ),
    DESeq2_comparison_FDR05 = sum(ds_sig),
    edgeR_FDR05 = sum(er_sig),
    overlap_FDR05 = sum(both),
    overlap_same_direction = sum(both & tab$same_effect_direction),
    significance_Jaccard = if (sum(either)) sum(both) / sum(either) else NA,
    log2FC_Spearman = cor(
      tab$DESeq2_log2FC[valid], tab$edgeR_log2FC[valid],
      method = "spearman"
    ),
    effect_direction_agreement = mean(tab$same_effect_direction[valid]),
    pvalue_rank_Spearman = cor(
      tab$DESeq2_pvalue[valid], tab$edgeR_pvalue[valid],
      method = "spearman"
    )
  )

  plot(
    tab$DESeq2_log2FC[valid], tab$edgeR_log2FC[valid],
    pch = 16, cex = 0.35,
    col = ifelse(both[valid], "#C83E4D", "#77777755"),
    xlab = "DESeq2 log2 fold change",
    ylab = "edgeR log2 fold change",
    main = paste(name, "- T2D versus ND")
  )
  abline(0, 1, col = "#2878B5", lty = 2)
}
dev.off()

summary <- do.call(rbind, summaries)
write.table(
  summary, file.path(out, "method_comparison_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()), file.path(out, "sessionInfo.txt")
)
print(summary, row.names = FALSE)
