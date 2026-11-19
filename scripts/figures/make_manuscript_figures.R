#!/usr/bin/env Rscript
# Data-driven manuscript figures for t2d-islet-regulatory-genomics.
# Run from the repository root:
#   Rscript scripts/figures/make_manuscript_figures.R 2>&1 | tee logs/manuscript_figures.log
#
# Produces (PDF, PNG and TIFF at 300 dpi) in figures/manuscript/:
#   Figure2_volcano_and_method_concordance
#   Figure4B_excitability_donor_expression
#   Figure5B_proteasome_per_gene_heatmap
#   Figure5C_proteasome_module_score_by_donor
#   FigureS1_donor_PCA
#
# The script stops with a clear message if a file or column is missing,
# rather than guessing. Edit the PATHS block if your files live elsewhere.

suppressPackageStartupMessages({
  library(ggplot2)
  library(edgeR)
})

# ---------------- PATHS ----------------
counts_path  <- "data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz"
samples_path <- "metadata/GSE221156_beta_pseudobulk_samples.tsv"
res_dir      <- "results/expression/GSE221156"
res_files    <- c(Primary = "primary_DESeq2_edgeR.tsv",
                  `BMI-adjusted` = "BMI_sensitivity_DESeq2_edgeR.tsv",
                  `V3 only` = "V3_only_DESeq2_edgeR.tsv")
core_path    <- "metadata/GSE221156_expression/core_supported_genes.tsv"
module_path  <- "config/gene_sets/proteasome29.tsv"
out_dir      <- "figures/manuscript"
excit_genes  <- c("KCNA1", "KCNAB1", "KCNAB2", "KCNH8", "KCNJ2", "NTRK2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
need <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)
for (p in c(counts_path, samples_path, core_path, module_path,
            file.path(res_dir, res_files))) need(p)

# Okabe-Ito palette, matching the summary figures
pal_dis  <- c(ND = "#0072B2", PD = "#009E73", T2D = "#D55E00")
shape_ch <- c(V2 = 17, V3 = 16)
theme_pub <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        plot.tag = element_text(face = "bold", size = 11),
        legend.key.size = unit(3, "mm"))

save_fig <- function(p, name, w_mm, h_mm) {
  f <- file.path(out_dir, name)
  ggsave(paste0(f, ".pdf"), p, width = w_mm, height = h_mm, units = "mm")
  ggsave(paste0(f, ".png"), p, width = w_mm, height = h_mm, units = "mm", dpi = 300)
  ggsave(paste0(f, ".tiff"), p, width = w_mm, height = h_mm, units = "mm", dpi = 300,
         device = "tiff", compression = "lzw")
  message("Wrote ", f, ".{pdf,png,tiff}")
}

# ---------------- Load results ----------------
read_res <- function(f) {
  d <- read.delim(file.path(res_dir, f), check.names = FALSE)
  req <- c("gene_id", "gene_symbol", "DESeq2_log2FC", "DESeq2_FDR_no_independent_filter",
           "edgeR_log2FC", "edgeR_FDR")
  miss <- setdiff(req, names(d))
  if (length(miss)) stop(f, " lacks columns: ", paste(miss, collapse = ", "), call. = FALSE)
  d
}
res <- lapply(res_files, read_res)

core <- read.delim(core_path, check.names = FALSE)
if (!"gene_id" %in% names(core)) stop(core_path, " lacks a gene_id column", call. = FALSE)
module <- read.delim(module_path, check.names = FALSE)
stopifnot("module file must hold 29 genes" = nrow(module) == 29L)

# ---------------- Figure 2: volcano + method concordance ----------------
pr <- res$Primary
pr$joint <- pr$DESeq2_FDR_no_independent_filter < 0.05 & pr$edgeR_FDR < 0.05 &
            sign(pr$DESeq2_log2FC) == sign(pr$edgeR_log2FC)
pr$class <- ifelse(pr$gene_id %in% core$gene_id, "Core (285)",
            ifelse(pr$joint, "Primary, both methods", "Not significant"))
pr$class <- factor(pr$class, levels = c("Not significant", "Primary, both methods", "Core (285)"))
pr <- pr[order(pr$class), ]
cls_col <- c("Not significant" = "grey80", "Primary, both methods" = "#56B4E9", "Core (285)" = "#D55E00")
message("Primary jointly significant (recomputed): ", sum(pr$joint, na.rm = TRUE),
        "; core genes found in table: ", sum(pr$class == "Core (285)"))

# Label: named candidates plus the 3 most significant core genes
top_core <- head(pr[pr$class == "Core (285)", ][order(pr$DESeq2_FDR_no_independent_filter[pr$class == "Core (285)"]), ], 3)
lab <- unique(rbind(pr[pr$gene_symbol %in% c(excit_genes, "MPP1", "PHLDB2", "DACT2"), ], top_core))
message("Volcano labels: ", paste(lab$gene_symbol, collapse = ", "))
# x-range from significant genes, with outliers outside it drawn at the edge
xlim_v <- max(abs(pr$DESeq2_log2FC[pr$class != "Not significant"]), na.rm = TRUE) * 1.1
pr$x_plot <- pmax(pmin(pr$DESeq2_log2FC, xlim_v), -xlim_v)
lab$x_plot <- lab$DESeq2_log2FC
n_clip <- sum(abs(pr$DESeq2_log2FC) > xlim_v, na.rm = TRUE)
if (n_clip) message(n_clip, " non-significant gene(s) drawn at the x-axis edge")

has_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!has_repel) message("ggrepel not installed: labels may overlap. install.packages('ggrepel')")
lab_layer <- if (has_repel) {
  ggrepel::geom_text_repel(data = lab, aes(label = gene_symbol), size = 2.3, fontface = "italic",
                           colour = "black", segment.size = 0.25, segment.colour = "grey40",
                           min.segment.length = 0, box.padding = 0.4, max.overlaps = Inf, seed = 1,
                           bg.color = "white", bg.r = 0.15)
} else {
  geom_text(data = lab, aes(label = gene_symbol), size = 2.3, fontface = "italic",
            colour = "black", vjust = -0.7)
}
pA <- ggplot(pr, aes(x_plot, -log10(DESeq2_FDR_no_independent_filter), colour = class)) +
  geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.25, linetype = 2, colour = "grey40") +
  geom_point(data = lab, size = 1.2, shape = 21, fill = NA, colour = "black", stroke = 0.3) +
  lab_layer +
  scale_colour_manual(values = cls_col, name = NULL) +
  scale_x_continuous(limits = c(-xlim_v, xlim_v)) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  labs(x = expression(log[2]~"fold-change, T2D vs ND (DESeq2)"),
       y = expression(-log[10]~"FDR (DESeq2)"), tag = "A") +
  theme_pub + theme(legend.position = "top")

rng <- c(-xlim_v, xlim_v)
pr$y_plot <- pmax(pmin(pr$edgeR_log2FC, xlim_v), -xlim_v)
pB <- ggplot(pr, aes(x_plot, y_plot, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = 2, colour = "grey40") +
  geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = cls_col, guide = "none") +
  coord_equal(xlim = rng, ylim = rng) +
  annotate("text", x = rng[1], y = rng[2], hjust = 0, vjust = 1, size = 2.5,
           label = sprintf("Spearman rho = %.3f",
                           cor(pr$DESeq2_log2FC, pr$edgeR_log2FC, method = "spearman",
                               use = "complete.obs"))) +
  labs(x = expression(log[2]~"fold-change (DESeq2)"),
       y = expression(log[2]~"fold-change (edgeR)"), tag = "B") +
  theme_pub

if (requireNamespace("patchwork", quietly = TRUE)) {
  save_fig(patchwork::wrap_plots(pA, pB, widths = c(1.2, 1)),
           "Figure2_volcano_and_method_concordance", 180, 85)
} else {
  message("patchwork not installed: saving Figure 2 panels separately")
  save_fig(pA, "Figure2A_volcano", 95, 85)
  save_fig(pB, "Figure2B_method_concordance", 85, 85)
}

# ---------------- Donor-level data ----------------
counts <- read.delim(counts_path, check.names = FALSE)
gene_col <- names(counts)[1]
rownames(counts) <- counts[[gene_col]]
counts <- as.matrix(counts[, -1, drop = FALSE])
samples <- read.delim(samples_path, check.names = FALSE)

# Detect donor, disease and chemistry columns from their contents
hit <- function(pred) names(samples)[vapply(samples, pred, logical(1))]
donor_col <- hit(function(x) all(colnames(counts) %in% as.character(x)))
disease_col <- hit(function(x) all(c("ND", "T2D") %in% as.character(x)))
chem_col <- hit(function(x) any(grepl("V2|v2", as.character(x))) && any(grepl("V3|v3", as.character(x))))
for (nm in c("donor_col", "disease_col", "chem_col")) {
  v <- get(nm)
  if (length(v) != 1) stop("Could not identify a unique ", nm, " in ", samples_path,
                           " (found: ", paste(v, collapse = ", "),
                           "). Set it by hand below this line.", call. = FALSE)
}
message("Using columns: donor=", donor_col, ", disease=", disease_col, ", chemistry=", chem_col)
samples <- samples[match(colnames(counts), samples[[donor_col]]), ]
stopifnot("sample order mismatch" = identical(as.character(samples[[donor_col]]), colnames(counts)))
meta <- data.frame(donor = colnames(counts),
                   disease = factor(samples[[disease_col]], levels = c("ND", "PD", "T2D")),
                   chemistry = ifelse(grepl("2", samples[[chem_col]]), "V2", "V3"))

# TMM log-CPM on genes passing the exploratory QC filter (>=10 counts in >=10 donors)
keep <- rowSums(counts >= 10) >= 10
y <- normLibSizes(DGEList(counts[keep, ]))
lcpm <- cpm(y, log = TRUE, prior.count = 2)
sym <- setNames(res$Primary$gene_symbol, res$Primary$gene_id)

long_expr <- function(ids) {
  do.call(rbind, lapply(ids, function(g) data.frame(meta, gene = sym[[g]], logCPM = lcpm[g, ])))
}

# ---------------- Figure 4B: excitability genes per donor ----------------
ids <- names(sym)[sym %in% excit_genes]
ids <- ids[ids %in% rownames(lcpm)]
if (length(ids) < length(excit_genes))
  message("Not all excitability genes passed the QC filter; plotting ", length(ids))
d4 <- long_expr(ids)
d4$gene <- factor(d4$gene, levels = excit_genes)
p4 <- ggplot(d4, aes(disease, logCPM, colour = disease)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, colour = "grey35") +
  geom_point(aes(shape = chemistry), position = position_jitter(width = 0.15, seed = 1),
             size = 1.1, alpha = 0.9) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal_dis, guide = "none") +
  scale_shape_manual(values = shape_ch, name = "Chemistry") +
  labs(x = NULL, y = "TMM log2 CPM", tag = "B",
       caption = "Each point is one donor. PD donors are shown for context and were not part of the tested contrast.") +
  theme_pub + theme(legend.position = "right", plot.caption = element_text(size = 6, hjust = 0))
save_fig(p4, "Figure4B_excitability_donor_expression", 180, 62)

# ---------------- Figure 5B: per-gene heatmap across models and methods ----------------
mod_ids <- module$gene_id
hm <- do.call(rbind, lapply(names(res), function(m) {
  d <- res[[m]][match(mod_ids, res[[m]]$gene_id), ]
  rbind(
    data.frame(gene = d$gene_symbol, model = m, method = "DESeq2", lfc = d$DESeq2_log2FC,
               sig = d$DESeq2_FDR_no_independent_filter < 0.05),
    data.frame(gene = d$gene_symbol, model = m, method = "edgeR", lfc = d$edgeR_log2FC,
               sig = d$edgeR_FDR < 0.05))
}))
ord <- res$Primary$gene_symbol[match(mod_ids, res$Primary$gene_id)][
  order(res$Primary$DESeq2_log2FC[match(mod_ids, res$Primary$gene_id)])]
hm$gene <- factor(hm$gene, levels = rev(ord))
hm$col <- factor(paste(hm$model, hm$method, sep = "\n"),
                 levels = as.vector(outer(names(res), c("DESeq2", "edgeR"), paste, sep = "\n"))[c(1, 4, 2, 5, 3, 6)])
lim <- max(abs(hm$lfc), na.rm = TRUE)
p5b <- ggplot(hm, aes(col, gene, fill = lfc)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = hm[hm$sig %in% TRUE, ], shape = 8, size = 0.8, colour = "black") +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       limits = c(-lim, lim), name = "log2 FC") +
  labs(x = NULL, y = NULL, tag = "B", caption = "* FDR < 0.05 in that method and model") +
  theme_pub + theme(axis.text.y = element_text(face = "italic", size = 6),
                    axis.line = element_blank(), axis.ticks = element_blank(),
                    plot.caption = element_text(size = 6))
save_fig(p5b, "Figure5B_proteasome_per_gene_heatmap", 95, 120)

# ---------------- Figure 5C: per-donor module score ----------------
ok <- mod_ids[mod_ids %in% rownames(lcpm)]
message("Module genes passing QC filter: ", length(ok), "/29")
z <- t(scale(t(lcpm[ok, , drop = FALSE])))
meta$module_score <- colMeans(z)
write.table(meta, file.path(out_dir, "Figure5C_module_scores_by_donor.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
p5c <- ggplot(meta, aes(disease, module_score, colour = disease)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, colour = "grey35") +
  geom_point(aes(shape = chemistry), position = position_jitter(width = 0.15, seed = 2),
             size = 1.4, alpha = 0.9) +
  scale_colour_manual(values = pal_dis, guide = "none") +
  scale_shape_manual(values = shape_ch, name = "Chemistry") +
  labs(x = NULL, y = "Module score\n(mean z-scored log2 CPM, 29 genes)", tag = "C",
       caption = "Descriptive only: genes were selected from these donors.") +
  theme_pub + theme(plot.caption = element_text(size = 6, hjust = 0))
save_fig(p5c, "Figure5C_proteasome_module_score_by_donor", 85, 75)

# ---------------- Figure S1: donor PCA ----------------
v <- apply(lcpm, 1, var)
pc <- prcomp(t(lcpm[order(v, decreasing = TRUE)[1:500], ]), center = TRUE, scale. = FALSE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
dpc <- data.frame(meta, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
pS1 <- ggplot(dpc, aes(PC1, PC2, colour = disease, shape = chemistry)) +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_colour_manual(values = pal_dis, name = "Disease") +
  scale_shape_manual(values = shape_ch, name = "Chemistry") +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
       caption = "TMM log2 CPM, 500 most variable genes, 48 donors.") +
  theme_pub + theme(plot.caption = element_text(size = 6, hjust = 0))
save_fig(pS1, "FigureS1_donor_PCA", 110, 85)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_manuscript_figures.txt"))
message("All figures written to ", out_dir)
