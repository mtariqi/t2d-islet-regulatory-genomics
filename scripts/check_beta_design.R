samples <- read.delim(
  "metadata/GSE221156_beta_pseudobulk_samples.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(samples) == 48,
  !anyDuplicated(samples$donor_id),
  all(samples$n_beta_cells >= 50),
  all(samples$total_counts > 0)
)

samples$disease_status <- factor(
  samples$disease_status, levels = c("ND", "PD", "T2D")
)
samples$chemistry <- factor(samples$chemistry)
samples$sex <- factor(samples$sex)
samples$age <- as.numeric(samples$age)
samples$BMI <- as.numeric(samples$BMI)

cat("Chemistry by disease\n")
print(with(samples, table(disease_status, chemistry)))

cat("\nSex by disease\n")
print(with(samples, table(disease_status, sex)))

cat("\nBeta-cell counts per donor\n")
print(summary(samples$n_beta_cells))

primary <- droplevels(
  samples[samples$disease_status %in% c("ND", "T2D"), ]
)
v3 <- droplevels(primary[primary$chemistry == "V3", ])

check_design <- function(name, data, formula) {
  variables <- all.vars(formula)
  if (anyNA(data[, variables, drop = FALSE])) {
    stop(paste(name, "has missing model covariates"))
  }

  design <- model.matrix(formula, data = data)
  rank <- qr(design)$rank

  data.frame(
    analysis = name,
    formula = paste(deparse(formula), collapse = " "),
    n_donors = nrow(design),
    n_parameters = ncol(design),
    rank = rank,
    residual_df = nrow(design) - rank,
    full_rank = rank == ncol(design)
  )
}

checks <- rbind(
  check_design(
    "primary_T2D_vs_ND", primary,
    ~ chemistry + sex + age + disease_status
  ),
  check_design(
    "BMI_sensitivity", primary,
    ~ chemistry + sex + age + BMI + disease_status
  ),
  check_design(
    "V3_only_sensitivity", v3,
    ~ sex + age + disease_status
  )
)

cat("\nDesign checks\n")
print(checks, row.names = FALSE)

write.table(
  checks, "metadata/GSE221156_beta_design_checks.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (!all(checks$full_rank) || any(checks$residual_df <= 0)) {
  stop("At least one design is not estimable")
}

cat("\nPASS: candidate designs are estimable.\n")
