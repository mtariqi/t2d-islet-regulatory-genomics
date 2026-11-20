#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/figures/make_figure6_replication.R
# Description:  Figure 6: forest plot of pre-registered replication scores by cohort, with pooled estimate.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================
#
# Run from the repository root, after the replication and meta-analysis scripts:
#   Rscript scripts/figures/make_figure6_replication.R
# Rows: each cohort's primary and sensitivity analyses, then the pre-registered pooled
# random-effects estimate (results/validation/meta/pooled_scores.tsv) if it exists.

suppressPackageStartupMessages(library(ggplot2))

files <- Sys.glob("results/validation/*/replication_scores.tsv")
if (!length(files)) stop("No results/validation/*/replication_scores.tsv found", call. = FALSE)
d <- do.call(rbind, lapply(files, function(f) cbind(cohort = basename(dirname(f)), read.delim(f))))

set_lab <- c(primary_min10 = "At least 10 beta cells (primary)",
             sens_min20 = "At least 20 beta cells",
             sens_all = "All donors",
             sens_unadjusted = "All donors, unadjusted")
d$row_lab <- unname(set_lab[d$donor_set])
adj <- d$cohort == "E-MTAB-5061" & d$donor_set == "primary_min10"
d$row_lab[adj] <- "At least 10 beta cells (primary),\nadjusted for sex and age"
d$row_lab <- paste0(d$row_lab, "\n", sprintf("%d ND / %d T2D", d$n_ND, d$n_T2D))
d$kind <- ifelse(d$donor_set == "primary_min10", "primary", "sensitivity")
d$order <- match(d$donor_set, names(set_lab))
keep <- c("cohort", "score", "estimate_T2D_minus_ND", "ci_low", "ci_high", "row_lab", "kind", "order")
e <- d[, keep]
e$I2 <- NA_real_

pooled_f <- "results/validation/meta/pooled_scores.tsv"
if (file.exists(pooled_f)) {
  p <- read.delim(pooled_f)
  e <- rbind(e, data.frame(cohort = "Pooled", score = p$score, estimate_T2D_minus_ND = p$RE_estimate,
                           ci_low = p$RE_ci_low, ci_high = p$RE_ci_high,
                           row_lab = sprintf("Random effects, %d cohorts", p$k),
                           kind = "pooled", order = 1, I2 = p$I2))
}

score_lab <- c(core279_direction_aware = "Core 279\n(direction-aware; expected higher)",
               kchannel5 = "Potassium channels, 5 genes\n(expected higher)",
               proteasome29 = "Proteasome module, 29 genes\n(expected lower)")
e$score_f <- factor(score_lab[e$score], levels = score_lab)
e$cohort <- factor(e$cohort, levels = intersect(c("GSE153855", "E-MTAB-5061", "Pooled"), unique(e$cohort)))
e$key <- paste0(e$cohort, "|", e$row_lab)
e$key <- factor(e$key, levels = unique(e$key[order(e$order, decreasing = TRUE)]))

g <- ggplot(e, aes(estimate_T2D_minus_ND, key)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.2,
                linewidth = 0.45, colour = "grey30") +
  geom_point(aes(shape = kind, fill = kind), size = 2.6, stroke = 0.4) +
  geom_label(data = e[e$kind == "pooled", ], fill = "white", label.size = 0, label.padding = unit(0.08, "lines"),
            aes(x = estimate_T2D_minus_ND, label = sprintf("I^2 == '%.0f%%'", I2)),
            parse = TRUE, vjust = -1.1, size = 2.3, colour = "grey30") +
  scale_shape_manual(values = c(primary = 21, sensitivity = 21, pooled = 23), guide = "none") +
  scale_fill_manual(values = c(primary = "#D55E00", sensitivity = "white", pooled = "black"), guide = "none") +
  scale_y_discrete(labels = function(x) sub("^[^|]*\\|", "", x)) +
  facet_grid(cohort ~ score_f, scales = "free_y", space = "free_y") +
  labs(x = "Score difference, T2D minus ND (z units, 95% CI)", y = NULL,
       caption = "Filled circles: primary analysis. Open circles: sensitivity analyses. Diamonds: pre-registered pooled estimate.") +
  theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.spacing.x = unit(4, "mm"), panel.spacing.y = unit(3, "mm"),
        axis.line = element_line(linewidth = 0.3), plot.caption = element_text(size = 6, hjust = 0))

out <- "figures/manuscript/Figure6_replication_forest"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
h <- 30 + 14 * length(unique(e$key))
ggsave(paste0(out, ".pdf"), g, width = 180, height = h, units = "mm")
ggsave(paste0(out, ".png"), g, width = 180, height = h, units = "mm", dpi = 300)
ggsave(paste0(out, ".tiff"), g, width = 180, height = h, units = "mm", dpi = 300,
       device = "tiff", compression = "lzw")
message("Wrote ", out, ".{pdf,png,tiff}: ", length(files), " cohort(s)",
        if (file.exists(pooled_f)) " + pooled estimate" else "")
