#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/figures/make_figure6_replication.R
# Description:  Figure 6: forest plot of pre-registered replication scores by cohort and donor set.
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
# Run from the repository root, after the replication script(s):
#   Rscript scripts/figures/make_figure6_replication.R
# Reads every results/validation/<cohort>/replication_scores.tsv it finds, so a second
# cohort is added to the figure automatically once its results exist.

suppressPackageStartupMessages(library(ggplot2))

files <- Sys.glob("results/validation/*/replication_scores.tsv")
if (!length(files)) stop("No results/validation/*/replication_scores.tsv found", call. = FALSE)
d <- do.call(rbind, lapply(files, function(f) cbind(cohort = basename(dirname(f)), read.delim(f))))

score_lab <- c(core279_direction_aware = "Core 279\n(direction-aware; expected higher)",
               kchannel5 = "Potassium channels, 5 genes\n(expected higher)",
               proteasome29 = "Proteasome module, 29 genes\n(expected lower)")
set_lab <- c(primary_min10 = "At least 10 beta cells (primary)",
             sens_min20 = "At least 20 beta cells",
             sens_all = "All donors",
             sens_unadjusted = "All donors, unadjusted")
d$score_f <- factor(score_lab[d$score], levels = score_lab)
d$label <- sprintf("%d ND / %d T2D", d$n_ND, d$n_T2D)
# Donor numbers go in the row label, not inside the panel, so they never overlap the intervals
d$y_lab <- paste0(set_lab[d$donor_set], "\n", d$label)
lv <- unique(d$y_lab[order(match(d$donor_set, names(set_lab)), decreasing = TRUE)])
d$set_f <- factor(d$y_lab, levels = lv)
d$primary <- d$donor_set == "primary_min10"

p <- ggplot(d, aes(estimate_T2D_minus_ND, set_f)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18, linewidth = 0.45, colour = "grey30") +
  geom_point(aes(fill = primary), shape = 21, size = 2.4, stroke = 0.4) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white"), guide = "none") +
  facet_grid(cohort ~ score_f, scales = "free_y") +
  labs(x = "Score difference, T2D minus ND (z units, 95% CI)", y = NULL,
       caption = "Pre-registered gene sets (tags prereg-validation-v1 and v2). Filled points: primary analysis. GSE153855: unadjusted (no covariates in GEO). E-MTAB-5061: primary adjusted for sex and age.") +
  theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.spacing.x = unit(4, "mm"), axis.line = element_line(linewidth = 0.3),
        plot.caption = element_text(size = 6, hjust = 0))

out <- "figures/manuscript/Figure6_replication_forest"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
h <- 30 + 32 * length(unique(d$cohort))
ggsave(paste0(out, ".pdf"), p, width = 180, height = h, units = "mm")
ggsave(paste0(out, ".png"), p, width = 180, height = h, units = "mm", dpi = 300)
ggsave(paste0(out, ".tiff"), p, width = 180, height = h, units = "mm", dpi = 300,
       device = "tiff", compression = "lzw")
message("Wrote ", out, ".{pdf,png,tiff} from ", length(files), " cohort(s)")
