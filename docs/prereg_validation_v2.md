# Pre-registration v2: core-gene and potassium-channel replication

Written after the discovery analyses and the disease x chemistry interaction test,
and before any replication expression data were analysed.

## Frozen gene sets
- config/gene_sets/core279_strict.tsv: 279 strictly LODO-supported core genes with expected direction
- config/gene_sets/kchannel5.tsv: KCNA1, KCNAB1, KCNAB2, KCNH8, KCNJ2 (expected up in T2D)
- SHA256 checksums: config/gene_sets/prereg_v2.sha256

## Hypotheses
- H-R1 (primary): a direction-aware core score (mean z of up genes minus mean z of down genes,
  donor-level) is higher in T2D than ND donors.
- H-R2 (secondary): the potassium-channel score (mean z of the five genes) is higher in T2D.

## Cohorts and analysis
Same cohorts, donor-overlap audit, donor-level aggregation and covariate rules as
docs/Independent_Validation_Plan.md. Genes not detected in a cohort are dropped, not imputed;
the number used is reported. Primary test per cohort: linear model of the score on disease,
adjusted for sex and age where available. Pooled across cohorts by random-effects meta-analysis.
Secondary: per-gene direction agreement with discovery, reported descriptively.

## Decision rules
- Replicated: pooled T2D effect positive with 95% CI excluding 0.
- Directionally consistent, not replicated: pooled effect positive, CI includes 0.
- Not replicated: pooled effect negative or near 0; reported prominently.
Platform differences (Smart-seq2 vs 10x) are noted; a null result does not prove absence of effect.

## Addendum (before any GSE153855 expression was examined)
Cell annotation inspected only (donor, disease, cell type counts). Beta cells per donor: 9-100.
The discovery 50-cell minimum is not used for Smart-seq2 data.
- Primary: donors with >= 10 annotated beta cells (excludes Sample107; 6 ND vs 4 T2D).
- Sensitivity 1: >= 20 beta cells (5 ND vs 3 T2D).
- Sensitivity 2: all 11 donors.
Cells are linked to annotation by column position (no cell IDs in either file);
the script checks that matrix columns equal annotation rows before use.

## Addendum: E-MTAB-5061 (before any expression was examined)
Written after inspecting the SDRF annotation only; the count file had not finished downloading.
- Cells: "Characteristics [single cell quality]" == "OK" exactly (authors' final QC; excludes
  "OK, filtered" and "not OK") and "Characteristics [inferred cell type]" == "beta cell". 270 beta cells.
- Cells are linked to counts by name (SDRF "Source Name", 3,514 unique). Only the read-count
  columns of pancreas_refseq_rpkms_counts_3514sc.txt are used; RPKM columns are never summed.
- Donor threshold: >= 10 beta cells (primary) keeps all 10 donors (6 ND, 4 T2D).
  The >= 20 set leaves 4 ND vs 1 T2D and is reported descriptively only (not tested).
- Primary model: score ~ disease + sex + age (age z-scored), as the plan requires where
  covariates exist. Sensitivity: unadjusted score ~ disease. Age is higher in T2D donors
  (37-57 vs 22-48), so the adjusted and unadjusted estimates are both reported.
- Independence from GSE153855: both studies are from Karolinska Institutet. If the GSE153855
  publications (PubMed 35482056, 39919147) state that their donors are distinct from the
  E-MTAB-5061 donors (PubMed 27667667), both cohorts enter the pooled random-effects
  meta-analysis. If they share donors, or this cannot be confirmed, E-MTAB-5061 is reported
  as a separate cohort and not pooled.
