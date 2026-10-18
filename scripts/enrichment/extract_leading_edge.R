de <- read.delim(
  "results/expression/GSE221156/primary_DESeq2_edgeR.tsv"
)

core <- read.delim(
  "results/influence/GSE221156/core285_lodo_assessment.tsv"
)

fg <- read.delim(
  "results/enrichment/GSE221156/all_ranked_enrichment.tsv"
)

selected <- c(
  "REACTOME_TRANSLATION",
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_POTASSIUM_CHANNELS",
  "REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS",
  "GOBP_SYNAPSE_ASSEMBLY",
  "GOBP_ACTION_POTENTIAL"
)

fg <- fg[
  fg$pathway %in% selected &
  fg$method == "DESeq2",
]

stopifnot(nrow(fg) == length(selected))

out <- list()

for (i in seq_len(nrow(fg))) {
  ids <- strsplit(
    fg$leadingEdge[i], ";",
    fixed = TRUE
  )[[1]]

  z <- de[match(ids, de$gene_id), ]

  z$pathway <- fg$pathway[i]
  z$NES <- fg$NES[i]
  z$pathway_FDR <- fg$padj[i]

  z$core285 <- z$gene_id %in% core$gene_id
  z$strict_LODO_supported <- z$gene_id %in%
    core$gene_id[core$strict_LODO_supported]

  out[[i]] <- z
}

out <- do.call(rbind, out)

write.table(
  out,
  "results/enrichment/GSE221156/selected_pathway_leading_edges.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n=== LEADING-EDGE SUMMARY ===\n")

summary <- aggregate(
  cbind(core285, strict_LODO_supported) ~ pathway,
  data = out,
  FUN = sum
)

summary$leading_edge_genes <- as.integer(
  table(factor(
    out$pathway,
    levels = summary$pathway
  ))
)

print(summary, row.names = FALSE)

cat("\n=== CORE GENES IN SELECTED PATHWAYS ===\n")

z <- out[out$core285, c(
  "pathway", "gene_symbol",
  "DESeq2_log2FC",
  "DESeq2_FDR_no_independent_filter",
  "strict_LODO_supported"
)]

print(z, row.names = FALSE)
