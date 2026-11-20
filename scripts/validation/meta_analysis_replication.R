#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/validation/meta_analysis_replication.R
# Description:  Pre-registered random-effects meta-analysis of replication scores across cohorts.
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
# Run from the repository root, after both replication scripts:
#   Rscript scripts/validation/meta_analysis_replication.R 2>&1 | tee logs/meta_analysis_replication.log
#
# Written and committed before the E-MTAB-5061 result and before any pooled estimate existed.
#
# Inputs: results/validation/<cohort>/replication_scores.tsv, primary rows (donor_set primary_min10).
# Effect: T2D minus ND difference in score (z units, scored within each cohort).
# SE: recovered from each cohort's 95% CI with t quantiles on that model's residual df:
#     df = donors - n_params, n_params = 2 (unadjusted: intercept, disease) or
#                                        4 (adjusted: intercept, disease, sex, age).
# Primary: random-effects meta-analysis (metafor::rma, REML) per score.
# Also reported: fixed-effect estimate, I^2, and Hartung-Knapp CI. With two cohorts,
# between-cohort heterogeneity cannot be estimated reliably; the HK interval shows that uncertainty.
#
# Independence gate (prereg addendum, commit a61a07f): pooled estimates are reported only if
# INDEPENDENCE_CONFIRMED is TRUE in config/validation_independence.txt. Otherwise the script
# writes the cohorts side by side and states that pooling was not done.

if (!requireNamespace("metafor", quietly = TRUE))
  stop("Install metafor first: install.packages('metafor')", call. = FALSE)

cohort_params <- c(GSE153855 = 2, `E-MTAB-5061` = 4)  # model parameters in each cohort's primary model
gate_file <- "config/validation_independence.txt"
out_dir <- "results/validation/meta"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

independent <- FALSE
if (file.exists(gate_file)) {
  g <- trimws(readLines(gate_file, warn = FALSE))
  independent <- any(g == "INDEPENDENCE_CONFIRMED=TRUE")
}
cat("Independence confirmed (", gate_file, "):", independent, "\n")

files <- Sys.glob("results/validation/*/replication_scores.tsv")
files <- files[basename(dirname(files)) %in% names(cohort_params)]
if (length(files) < 2) stop("Need results from both cohorts; found: ", paste(files, collapse = ", "), call. = FALSE)
d <- do.call(rbind, lapply(files, function(f) cbind(cohort = basename(dirname(f)), read.delim(f))))
d <- d[d$donor_set == "primary_min10", ]
d$n <- d$n_ND + d$n_T2D
d$df <- d$n - cohort_params[d$cohort]
d$se <- (d$ci_high - d$ci_low) / (2 * qt(0.975, d$df))
write.table(d, file.path(out_dir, "cohort_inputs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nCohort inputs (primary analyses):\n")
print(d[, c("cohort", "score", "n_ND", "n_T2D", "df", "estimate_T2D_minus_ND", "se", "p_two_sided")],
      row.names = FALSE, digits = 3)

if (!independent) {
  cat("\nPooling NOT done: independence between cohorts is not confirmed.\n",
      "Report the cohorts side by side (cohort_inputs.tsv).\n", sep = "")
  quit(save = "no", status = 0)
}

expect <- c(core279_direction_aware = "higher", kchannel5 = "higher", proteasome29 = "lower")
res <- lapply(names(expect), function(s) {
  x <- d[d$score == s, ]
  re <- metafor::rma(yi = x$estimate_T2D_minus_ND, sei = x$se, method = "REML")
  fe <- metafor::rma(yi = x$estimate_T2D_minus_ND, sei = x$se, method = "FE")
  hk <- metafor::rma(yi = x$estimate_T2D_minus_ND, sei = x$se, method = "REML", test = "knha")
  sgn <- ifelse(expect[[s]] == "higher", 1, -1)
  excludes_zero <- re$ci.lb > 0 || re$ci.ub < 0
  verdict <- if (sign(re$b[1]) == sgn && excludes_zero) "Replicated" else
             if (sign(re$b[1]) == sgn) "Directionally consistent, not replicated" else "Not replicated"
  data.frame(score = s, expected = expect[[s]], k = nrow(x),
             RE_estimate = re$b[1], RE_ci_low = re$ci.lb, RE_ci_high = re$ci.ub, RE_p = re$pval,
             FE_estimate = fe$b[1], FE_ci_low = fe$ci.lb, FE_ci_high = fe$ci.ub,
             HK_ci_low = hk$ci.lb, HK_ci_high = hk$ci.ub, I2 = re$I2, tau2 = re$tau2,
             verdict_prereg = verdict)
})
res <- do.call(rbind, res)
write.table(res, file.path(out_dir, "pooled_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== POOLED (primary: random effects, REML) ===\n")
print(res, row.names = FALSE, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
