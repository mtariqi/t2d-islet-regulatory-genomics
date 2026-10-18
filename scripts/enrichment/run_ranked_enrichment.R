suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

set.seed(20261018)

input <- "results/expression/GSE221156/primary_DESeq2_edgeR.tsv"
outdir <- "results/enrichment/GSE221156"
figdir <- "figures/phase8"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

de <- read.delim(input, check.names = FALSE)

required <- c(
  "gene_id", "gene_symbol",
  "DESeq2_log2FC", "DESeq2_pvalue",
  "edgeR_log2FC", "edgeR_pvalue"
)

stopifnot(all(required %in% names(de)))
stopifnot(!anyDuplicated(de$gene_id))

# Signed two-sided p-value transformation.
# Scores are NOT original DESeq2 Wald / edgeR QL statistics.
make_rank <- function(logfc, pvalue, gene_id) {
  valid <- is.finite(logfc) &
    is.finite(pvalue) &
    pvalue >= 0 & pvalue <= 1 &
    !is.na(gene_id) & nzchar(gene_id)

  logfc <- logfc[valid]
  pvalue <- pvalue[valid]
  gene_id <- gene_id[valid]

  positive <- pvalue[pvalue > 0]
  floor_p <- if (length(positive)) {
    min(positive) / 10
  } else {
    .Machine$double.xmin
  }

  pvalue <- pmax(pvalue, floor_p)

  score <- sign(logfc) * qnorm(
    pvalue / 2,
    lower.tail = FALSE
  )

  names(score) <- gene_id
  score <- sort(score, decreasing = TRUE)

  stopifnot(!anyDuplicated(names(score)))
  stopifnot(all(is.finite(score)))

  score
}

ranks <- list(
  DESeq2 = make_rank(
    de$DESeq2_log2FC, de$DESeq2_pvalue, de$gene_id
  ),
  edgeR = make_rank(
    de$edgeR_log2FC, de$edgeR_pvalue, de$gene_id
  )
)

# Use MSigDB human gene sets, with Ensembl mapping.
get_sets <- function(collection, subcollection = NULL) {
  args <- list(
    species = "Homo sapiens",
    collection = collection
  )

  if (!is.null(subcollection)) {
    args$subcollection <- subcollection
  }

  ms <- do.call(msigdbr::msigdbr, args)

  stopifnot(
    all(c("gs_name", "ensembl_gene") %in% names(ms))
  )

  ms <- ms[
    !is.na(ms$ensembl_gene) &
      nzchar(ms$ensembl_gene),
    c("gs_name", "ensembl_gene")
  ]

  ms <- unique(ms)

  split(ms$ensembl_gene, ms$gs_name)
}

collections <- list(
  GO_BP = get_sets("C5", "GO:BP"),
  Reactome = get_sets("C2", "CP:REACTOME")
)

all_results <- list()

for (collection_name in names(collections)) {
  pathways <- collections[[collection_name]]

  for (method in names(ranks)) {
    stats <- ranks[[method]]

    cat(
      "\nRunning:", collection_name, method,
      "Ranked genes:", length(stats),
      "\n"
    )

    fg <- fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      minSize = 15,
      maxSize = 500,
      eps = 0
    )

    fg <- as.data.frame(fg)
    if (!nrow(fg)) next

    fg$method <- method
    fg$collection <- collection_name

    # Convert leading-edge vectors into TSV-safe strings.
    fg$leadingEdge <- vapply(
      fg$leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )

    fg <- fg[
      order(fg$padj, -abs(fg$NES)),
    ]

    filename <- file.path(
      outdir,
      paste0(collection_name, "_", method, "_fgsea.tsv")
    )

    write.table(
      fg, filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    all_results[[paste(collection_name, method)]] <- fg

    cat(
      "Tested pathways:", nrow(fg),
      "FDR < 0.05:",
      sum(fg$padj < 0.05, na.rm = TRUE),
      "\n"
    )
  }
}

combined <- do.call(rbind, all_results)

if (is.null(combined) || !nrow(combined)) {
  stop("No enrichment results returned.")
}

write.table(
  combined,
  file.path(outdir, "all_ranked_enrichment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Identify pathways significant in both methods,
# with the same direction of enrichment.
consensus <- list()

for (collection_name in names(collections)) {
  a <- subset(
    combined,
    collection == collection_name & method == "DESeq2"
  )
  b <- subset(
    combined,
    collection == collection_name & method == "edgeR"
  )

  merged <- merge(
    a[, c("pathway", "NES", "padj", "size")],
    b[, c("pathway", "NES", "padj", "size")],
    by = "pathway",
    suffixes = c("_DESeq2", "_edgeR")
  )

  merged$direction_consistent <- (
    sign(merged$NES_DESeq2) ==
      sign(merged$NES_edgeR)
  )

  merged$consensus_supported <- (
    merged$padj_DESeq2 < 0.05 &
    merged$padj_edgeR < 0.05 &
    merged$direction_consistent
  )

  merged$collection <- collection_name

  consensus[[collection_name]] <- merged

  cat(
    "Consensus", collection_name, ":",
    sum(merged$consensus_supported, na.rm = TRUE),
    "\n"
  )
}

consensus <- do.call(rbind, consensus)

write.table(
  consensus,
  file.path(outdir, "pathway_method_consensus.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Visualize up to 20 significant consensus pathways per collection.
for (collection_name in unique(consensus$collection)) {
  x <- subset(
    consensus,
    collection == collection_name & consensus_supported
  )

  if (!nrow(x)) next

  x <- x[
    order(pmax(x$padj_DESeq2, x$padj_edgeR)),
  ]

  x <- head(x, 20)

  plotdata <- data.frame(
    pathway = rep(x$pathway, 2),
    method = rep(c("DESeq2", "edgeR"), each = nrow(x)),
    NES = c(x$NES_DESeq2, x$NES_edgeR)
  )

  plotdata$pathway <- factor(
    plotdata$pathway,
    levels = rev(x$pathway)
  )

  p <- ggplot(
    plotdata,
    aes(x = NES, y = pathway, color = method)
  ) +
    geom_point(size = 2.8) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(
      title = paste(collection_name, "cross-method enrichment"),
      x = "Normalized enrichment score",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(
      figdir,
      paste0(collection_name, "_consensus_NES.pdf")
    ),
    p, width = 11, height = 8
  )

  ggsave(
    file.path(
      figdir,
      paste0(collection_name, "_consensus_NES.png")
    ),
    p, width = 11, height = 8, dpi = 300
  )
}

sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\nPhase 8 primary enrichment completed.\n")
