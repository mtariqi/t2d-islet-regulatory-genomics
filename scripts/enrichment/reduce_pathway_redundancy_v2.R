base <- "results/enrichment/GSE221156"

x <- read.delim(
  file.path(base, "pathway_redundancy_clusters.tsv"),
  check.names = FALSE
)

stopifnot(all(c(
  "collection", "pathway", "leadingEdge",
  "NES_DESeq2", "NES_edgeR", "worst_FDR"
) %in% names(x)))

parse_genes <- function(s) {
  if (is.na(s) || !nzchar(s)) return(character())
  unique(trimws(strsplit(s, ";", fixed = TRUE)[[1]]))
}

cluster_group <- function(z, threshold = 0.5) {
  n <- nrow(z)

  if (n == 1L) {
    z$new_cluster <- 1L
    return(z)
  }

  sets <- lapply(z$leadingEdge, parse_genes)

  sim <- diag(1, n)

  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      a <- sets[[i]]
      b <- sets[[j]]

      u <- length(union(a, b))

      value <- if (u == 0L) 0 else {
        length(intersect(a, b)) / u
      }

      sim[i, j] <- value
      sim[j, i] <- value
    }
  }

  d <- as.dist(1 - sim)

  hc <- hclust(d, method = "complete")

  z$new_cluster <- cutree(
    hc,
    h = 1 - threshold
  )

  z
}

x$direction <- ifelse(
  x$NES_DESeq2 > 0, "positive", "negative"
)

stopifnot(all(
  sign(x$NES_DESeq2) == sign(x$NES_edgeR)
))

groups <- split(
  x,
  interaction(x$collection, x$direction, drop = TRUE)
)

out <- lapply(names(groups), function(group_name) {
  z <- cluster_group(groups[[group_name]])

  z$cluster_v2 <- paste0(
    group_name, "_", z$new_cluster
  )

  z
})

out <- do.call(rbind, out)

out <- out[order(out$worst_FDR), ]

out$representative_v2 <- !duplicated(out$cluster_v2)

write.table(
  out,
  file.path(base, "pathway_redundancy_clusters_v2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

representatives <- out[out$representative_v2, ]

write.table(
  representatives,
  file.path(base, "representative_pathways_v2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("=== PHASE 8 COMPLETE-LINKAGE REDUNDANCY ===\n")

cat("Input pathways:", nrow(out), "\n")
cat("Clusters:", length(unique(out$cluster_v2)), "\n")

sizes <- sort(
  table(out$cluster_v2),
  decreasing = TRUE
)

cat("Largest cluster:", max(sizes), "\n")
cat("Singleton clusters:", sum(sizes == 1), "\n")

cat("\nClusters by collection and direction:\n")

print(
  aggregate(
    cluster_v2 ~ collection + direction,
    data = out,
    FUN = function(v) length(unique(v))
  ),
  row.names = FALSE
)

cat("\nLargest 15 clusters:\n")
print(head(sizes, 15))

cat("\nTop positive representatives:\n")

p <- representatives[
  representatives$direction == "positive",
]

p <- p[order(p$worst_FDR), ]

print(
  head(
    p[, c(
      "collection", "pathway",
      "NES_DESeq2", "NES_edgeR",
      "worst_FDR"
    )],
    15
  ),
  row.names = FALSE
)
