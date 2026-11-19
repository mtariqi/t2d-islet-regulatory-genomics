# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/enrichment/audit_leading_edge_membership.R
# Description:  Audit leading-edge gene membership across the pathways of each redundancy cluster.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================

suppressPackageStartupMessages({
  library(msigdbr)
})

base <- "results/enrichment/GSE221156"

clusters <- read.delim(
  file.path(base, "pathway_redundancy_clusters_v2.tsv"),
  check.names = FALSE
)

fg <- read.delim(
  file.path(base, "all_ranked_enrichment.tsv"),
  check.names = FALSE
)

z <- clusters[
  clusters$cluster_v2 == "Reactome.negative_18",
]

stopifnot(nrow(z) == 38L)

parse_genes <- function(s) {
  if (is.na(s) || !nzchar(s)) return(character())

  unique(trimws(
    strsplit(s, ";", fixed = TRUE)[[1]]
  ))
}

# Obtain Reactome definitions using the same
# human MSigDB collection used in Phase 8.
msig <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)

stopifnot(all(c(
  "gs_name", "ensembl_gene"
) %in% names(msig)))

pathway_sets <- split(
  msig$ensembl_gene,
  msig$gs_name
)

original <- fg[
  fg$collection == "Reactome" &
  fg$method == "DESeq2",
]

original <- original[
  match(z$pathway, original$pathway),
]

stopifnot(
  nrow(original) == nrow(z),
  identical(
    as.character(original$pathway),
    as.character(z$pathway)
  )
)

audit <- lapply(seq_len(nrow(z)), function(i) {
  pathway <- z$pathway[i]

  leading <- parse_genes(z$leadingEdge[i])
  original_leading <- parse_genes(original$leadingEdge[i])

  members <- unique(pathway_sets[[pathway]])

  data.frame(
    pathway = pathway,
    leading_edge_n = length(leading),
    original_matches = setequal(
      leading,
      original_leading
    ),
    pathway_found = !is.null(members),
    nonmember_n = if (is.null(members)) NA_integer_
                  else sum(!leading %in% members),
    stringsAsFactors = FALSE
  )
})

audit <- do.call(rbind, audit)

cat("=== LEADING-EDGE MEMBERSHIP AUDIT ===\n")

cat("Pathways:", nrow(audit), "\n")

cat(
  "Original leading-edge matches:",
  sum(audit$original_matches),
  "/", nrow(audit), "\n"
)

cat(
  "Pathways found in MSigDB:",
  sum(audit$pathway_found),
  "/", nrow(audit), "\n"
)

cat(
  "Pathways with nonmember genes:",
  sum(audit$nonmember_n > 0, na.rm = TRUE),
  "\n"
)

cat("\nProblematic pathways:\n")

problems <- audit[
  !audit$original_matches |
  !audit$pathway_found |
  (!is.na(audit$nonmember_n) &
   audit$nonmember_n > 0),
]

print(problems, row.names = FALSE)

write.table(
  audit,
  file.path(
    base,
    "Reactome_negative18_membership_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
