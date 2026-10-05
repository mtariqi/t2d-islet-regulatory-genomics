suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(uwot))

out <- "results/qc/GSE221156"
vsd <- readRDS(file.path(out, "beta_QC_vst.rds"))
expression <- assay(vsd)
metadata <- as.data.frame(colData(vsd))

colors <- c(ND = "#2878B5", PD = "#E6A000", T2D = "#C83E4D")
coordinates <- list()

pdf(file.path(out, "donor_UMAP_sensitivity.pdf"), width = 9, height = 7)

for (subset_name in c("all_donors", "V3_only")) {
  selected <- if (subset_name == "all_donors") {
    rep(TRUE, ncol(expression))
  } else {
    metadata$chemistry == "V3"
  }

  x <- expression[, selected, drop = FALSE]
  m <- metadata[selected, , drop = FALSE]
  variable_genes <- head(
    order(apply(x, 1, var), decreasing = TRUE), 500
  )
  pcs <- prcomp(
    t(x[variable_genes, , drop = FALSE]),
    center = TRUE, scale. = FALSE
  )$x[, 1:10, drop = FALSE]

  for (neighbors in c(5, 10, 15)) {
    for (seed in c(11, 42, 2026)) {
      set.seed(seed)
      embedding <- uwot::umap(
        pcs,
        n_neighbors = neighbors,
        min_dist = 0.3,
        metric = "euclidean",
        n_threads = 1,
        n_sgd_threads = 1,
        verbose = FALSE
      )

      key <- paste(subset_name, neighbors, seed, sep = "_")
      coordinates[[key]] <- data.frame(
        subset = subset_name,
        n_neighbors = neighbors,
        seed = seed,
        donor_id = rownames(m),
        disease_status = m$disease_status,
        chemistry = m$chemistry,
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2]
      )

      plot(
        embedding,
        col = colors[as.character(m$disease_status)],
        pch = ifelse(m$chemistry == "V2", 17, 16),
        xlab = "UMAP1", ylab = "UMAP2",
        main = sprintf(
          "%s | neighbors=%d | seed=%d",
          subset_name, neighbors, seed
        )
      )
      legend(
        "topright", legend = names(colors),
        col = colors, pch = 16, bty = "n"
      )
      if (subset_name == "all_donors") {
        legend(
          "bottomright", legend = c("V2", "V3"),
          pch = c(17, 16), bty = "n"
        )
      }
    }
  }
}
dev.off()

write.table(
  do.call(rbind, coordinates),
  file.path(out, "donor_UMAP_coordinates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out, "UMAP_sessionInfo.txt")
)
cat("Saved donor_UMAP_sensitivity.pdf and coordinates.\n")
