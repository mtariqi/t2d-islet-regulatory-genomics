base <- "results/enrichment/GSE221156"

x <- read.delim(
  file.path(base, "pathway_redundancy_clusters_v2.tsv")
)

de <- read.delim(
  "results/expression/GSE221156/primary_DESeq2_edgeR.tsv"
)

stopifnot(all(c(
  "gene_id",
  "gene_symbol",
  "DESeq2_log2FC",
  "DESeq2_FDR_native",
  "edgeR_log2FC",
  "edgeR_FDR"
) %in% names(de)))

z <- x[x$cluster_v2 == "Reactome.negative_18", ]

stopifnot(nrow(z) == 38L)

parse_genes <- function(s) {
  unique(trimws(
    strsplit(s, ";", fixed = TRUE)[[1]]
  ))
}

sets <- lapply(z$leadingEdge, parse_genes)

freq <- sort(
  table(unlist(sets)),
  decreasing = TRUE
)

shared <- names(freq)[freq == nrow(z)]

stopifnot(length(shared) == 29L)

out <- data.frame(
  gene_id = shared,
  pathway_frequency = as.integer(freq[shared])
)

out <- merge(
  out,
  de,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

out <- out[match(shared, out$gene_id), ]

out$significant_both_recomputed <- with(
  out,
  !is.na(DESeq2_FDR_native) &
  !is.na(edgeR_FDR) &
  DESeq2_FDR_native < 0.05 &
  edgeR_FDR < 0.05
)

out$negative_both <- with(
  out,
  DESeq2_log2FC < 0 &
  edgeR_log2FC < 0
)

outfile <- file.path(
  base,
  "Reactome_negative18_shared29_gene_characterization.tsv"
)

write.table(
  out,
  outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("=== SHARED REACTOME GENE CHARACTERIZATION ===\n")

cat("Shared genes:", nrow(out), "\n")

cat(
  "Matched to primary DE table:",
  sum(!is.na(out$DESeq2_log2FC)),
  "\n"
)

cat(
  "Significant in both methods:",
  sum(out$significant_both_recomputed, na.rm = TRUE),
  "\n"
)

cat(
  "Negative log2FC in both methods:",
  sum(out$negative_both, na.rm = TRUE),
  "\n"
)

cat("\nGene-level results:\n")

print(
  out[, c(
    "gene_id",
    "gene_symbol",
    "DESeq2_log2FC",
    "DESeq2_FDR_native",
    "edgeR_log2FC",
    "edgeR_FDR"
  )],
  row.names = FALSE
)

cat("\nSaved:", outfile, "\n")
