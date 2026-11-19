# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/enrichment/audit_phase8.R
# Description:  Audit Phase 8 enrichment outputs for completeness and internal consistency.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================

x <- read.delim(
  "results/enrichment/GSE221156/pathway_method_consensus.tsv"
)

de <- read.delim(
  "results/expression/GSE221156/primary_DESeq2_edgeR.tsv"
)

stopifnot(
  all(c(
    "collection", "pathway",
    "NES_DESeq2", "NES_edgeR",
    "padj_DESeq2", "padj_edgeR",
    "consensus_supported"
  ) %in% names(x))
)

cat("\n=== PHASE 8 VALIDATION AUDIT ===\n")

cat("Primary tested genes:", nrow(de), "\n")
cat("Unique Ensembl IDs:", length(unique(de$gene_id)), "\n")
cat("Consensus pathways:", sum(x$consensus_supported), "\n")

cat("\nConsensus pathways by collection:\n")
print(table(x$collection[x$consensus_supported]))

cat("\nEnrichment direction:\n")
supported <- x[x$consensus_supported, ]

print(table(
  collection = supported$collection,
  direction = ifelse(
    supported$NES_DESeq2 > 0,
    "Positive",
    "Negative"
  )
))

cat("\nAgreement in NES across all tested pathways:\n")
for (collection in unique(x$collection)) {
  z <- x[x$collection == collection, ]

  cat(
    collection,
    "Pearson r =",
    round(cor(
      z$NES_DESeq2,
      z$NES_edgeR,
      use = "complete.obs"
    ), 4),
    "\n"
  )
}

cat("\nPrimary contrast direction check:\n")
print(head(
  de[, c(
    "gene_symbol",
    "DESeq2_log2FC",
    "edgeR_log2FC"
  )],
  10
))

cat("\nTop 10 positive consensus pathways:\n")
pos <- supported[supported$NES_DESeq2 > 0, ]

pos <- pos[order(
  pmax(pos$padj_DESeq2, pos$padj_edgeR)
), ]

print(head(
  pos[, c(
    "collection", "pathway",
    "NES_DESeq2", "NES_edgeR",
    "padj_DESeq2", "padj_edgeR"
  )],
  10
), row.names = FALSE)
