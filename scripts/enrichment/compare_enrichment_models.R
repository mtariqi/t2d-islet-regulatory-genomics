base <- "results/enrichment/GSE221156"

models <- c("primary", "BMI_sensitivity", "V3_only")

read_model <- function(model) {
  directory <- if (model == "primary") {
    base
  } else {
    file.path(base, model)
  }

  f <- file.path(directory, "pathway_method_consensus.tsv")
  stopifnot(file.exists(f))

  d <- read.delim(f, check.names = FALSE)

  required <- c(
    "collection", "pathway",
    "NES_DESeq2", "NES_edgeR",
    "padj_DESeq2", "padj_edgeR",
    "consensus_supported"
  )

  stopifnot(all(required %in% names(d)))

  d$key <- paste(d$collection, d$pathway, sep = "::")
  stopifnot(!anyDuplicated(d$key))

  d
}

datasets <- lapply(models, read_model)
names(datasets) <- models

# Preserve the union of tested pathways.
keys <- unique(unlist(lapply(datasets, function(d) d$key)))

comparison <- data.frame(
  key = keys,
  stringsAsFactors = FALSE
)

for (model in models) {
  d <- datasets[[model]]
  idx <- match(comparison$key, d$key)

  for (field in c(
    "NES_DESeq2",
    "NES_edgeR",
    "padj_DESeq2",
    "padj_edgeR",
    "consensus_supported"
  )) {
    comparison[[paste(model, field, sep = "_")]] <- d[[field]][idx]
  }
}

primary <- datasets$primary

comparison$collection <- primary$collection[
  match(comparison$key, primary$key)
]

comparison$pathway <- primary$pathway[
  match(comparison$key, primary$key)
]

# Fill identifiers for pathways absent from primary.
missing <- is.na(comparison$collection)

comparison$collection[missing] <- sub(
  "::.*$", "", comparison$key[missing]
)

comparison$pathway[missing] <- sub(
  "^[^:]+::", "", comparison$key[missing]
)

supported <- function(model) {
  x <- comparison[[paste0(model, "_consensus_supported")]]
  !is.na(x) & x
}

comparison$primary_sig <- supported("primary")
comparison$BMI_sig <- supported("BMI_sensitivity")
comparison$V3_sig <- supported("V3_only")

comparison$significant_all_three <-
  comparison$primary_sig &
  comparison$BMI_sig &
  comparison$V3_sig

# Check direction across both methods and all models.
nes_columns <- unlist(lapply(
  models,
  function(m) paste0(m, "_NES_", c("DESeq2", "edgeR"))
))

nes <- as.matrix(comparison[, nes_columns])

comparison$direction_stable_all <- apply(nes, 1, function(v) {
  all(is.finite(v)) &&
    (all(v > 0) || all(v < 0))
})

cat("=== CROSS-MODEL ENRICHMENT ROBUSTNESS ===\n")

for (collection in c("GO_BP", "Reactome")) {
  x <- comparison[comparison$collection == collection, ]

  cat("\n", collection, "\n", sep = "")

  cat("Primary consensus:", sum(x$primary_sig), "\n")

  cat(
    "Primary retained in BMI:",
    sum(x$primary_sig & x$BMI_sig),
    "\n"
  )

  cat(
    "Primary retained in V3:",
    sum(x$primary_sig & x$V3_sig),
    "\n"
  )

  cat(
    "Significant in all three:",
    sum(x$significant_all_three),
    "\n"
  )

  cat(
    "Primary pathways with stable NES direction:",
    sum(x$primary_sig & x$direction_stable_all),
    "\n"
  )
}

# Analyze the 38-pathway Reactome cluster.
cluster_file <- file.path(
  base,
  "pathway_redundancy_clusters_v2.tsv"
)

clusters <- read.delim(cluster_file, check.names = FALSE)

cluster38 <- clusters[
  clusters$cluster_v2 == "Reactome.negative_18",
]

stopifnot(nrow(cluster38) == 38L)

target_keys <- paste(
  cluster38$collection,
  cluster38$pathway,
  sep = "::"
)

target <- comparison[
  match(target_keys, comparison$key),
]

stopifnot(nrow(target) == 38L)
stopifnot(!anyNA(target$key))

cat("\n=== PROTEASOME-ASSOCIATED CLUSTER ===\n")

for (model in models) {
  sig <- target[[paste0(model, "_consensus_supported")]]

  de_nes <- target[[paste0(model, "_NES_DESeq2")]]
  er_nes <- target[[paste0(model, "_NES_edgeR")]]

  cat("\nModel:", model, "\n")

  cat(
    "Significant in both methods:",
    sum(sig, na.rm = TRUE),
    "/38\n"
  )

  cat(
    "Negative NES in both methods:",
    sum(
      is.finite(de_nes) &
      is.finite(er_nes) &
      de_nes < 0 &
      er_nes < 0
    ),
    "/38\n"
  )

  cat(
    "Median DESeq2 NES:",
    median(de_nes, na.rm = TRUE),
    "\n"
  )

  cat(
    "Median edgeR NES:",
    median(er_nes, na.rm = TRUE),
    "\n"
  )
}

outfile <- file.path(
  base,
  "cross_model_pathway_robustness.tsv"
)

write.table(
  comparison,
  outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  target,
  file.path(base, "proteasome38_cross_model_robustness.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\nSaved:", outfile, "\n")
cat("Cross-model comparison completed.\n")
