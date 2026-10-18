base <- "results/enrichment/GSE221156"

x <- read.delim(
  file.path(base, "pathway_method_consensus.tsv")
)

fg <- read.delim(
  file.path(base, "all_ranked_enrichment.tsv")
)

stopifnot(all(c(
  "collection", "pathway", "consensus_supported",
  "padj_DESeq2", "padj_edgeR"
) %in% names(x)))

stopifnot(all(c(
  "collection", "pathway", "method", "leadingEdge"
) %in% names(fg)))

x <- x[!is.na(x$consensus_supported) &
         x$consensus_supported, ]

fg <- fg[fg$method == "DESeq2", ]

x <- merge(
  x, fg[, c("collection", "pathway", "leadingEdge")],
  by = c("collection", "pathway")
)

x$worst_FDR <- pmax(x$padj_DESeq2, x$padj_edgeR)

all_clusters <- list()

for (collection in unique(x$collection)) {

  z <- x[x$collection == collection, ]
  z <- z[order(z$worst_FDR), ]

  sets <- strsplit(
    as.character(z$leadingEdge),
    ";", fixed = TRUE
  )

  cluster <- rep(NA_integer_, nrow(z))
  cluster_id <- 0L

  for (i in seq_len(nrow(z))) {

    if (!is.na(cluster[i])) next

    cluster_id <- cluster_id + 1L
    cluster[i] <- cluster_id

    queue <- i

    while (length(queue)) {

      j <- queue[1]
      queue <- queue[-1]

      remaining <- which(is.na(cluster))

      for (k in remaining) {

        intersection <- length(intersect(
          sets[[j]], sets[[k]]
        ))

        union <- length(union(
          sets[[j]], sets[[k]]
        ))

        similarity <- if (union > 0) {
          intersection / union
        } else {
          0
        }

        if (similarity >= 0.5) {
          cluster[k] <- cluster_id
          queue <- c(queue, k)
        }
      }
    }
  }

  z$cluster_id <- paste0(collection, "_", cluster)
  z$representative <- !duplicated(z$cluster_id)

  all_clusters[[collection]] <- z

  cat(
    collection,
    ": ", nrow(z),
    "consensus pathways -> ",
    length(unique(z$cluster_id)),
    "overlap clusters\n",
    sep = ""
  )
}

out <- do.call(rbind, all_clusters)

write.table(
  out,
  file.path(base, "pathway_redundancy_clusters.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

representatives <- out[out$representative, ]

write.table(
  representatives,
  file.path(base, "representative_pathways.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nTop 20 representative pathways:\n")

representatives <- representatives[
  order(representatives$worst_FDR),
]

print(
  head(representatives[, c(
    "collection", "pathway",
    "cluster_id", "worst_FDR"
  )], 20),
  row.names = FALSE
)
