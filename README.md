# T2D Islet Regulatory Genomics

[![Status](https://img.shields.io/badge/Status-expression_analysis_complete-2878B5)](#current-results)
[![Genetic integration](https://img.shields.io/badge/Genetic_integration-planned-orange)](#research-plan)
[![Python](https://img.shields.io/badge/Python-h5py_%7C_NumPy-3776AB?logo=python&logoColor=white)](#tools-and-methods)
[![R](https://img.shields.io/badge/R-Bioconductor-276DC3?logo=r&logoColor=white)](#tools-and-methods)
[![DESeq2](https://img.shields.io/badge/DESeq2-primary_analysis-2878B5)](#statistical-design)
[![edgeR](https://img.shields.io/badge/edgeR-robust_QL_analysis-C83E4D)](#statistical-design)
[![UMAP](https://img.shields.io/badge/UMAP-parameter_sensitivity-7B4AB5)](#quality-control-and-umap)
[![GEO](https://img.shields.io/badge/GEO-GSE221156-008080)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156)
[![Donors](https://img.shields.io/badge/Donors-48-008080)](#dataset-and-provenance)
[![Beta cells](https://img.shields.io/badge/Beta_cells-99%2C029-008080)](#dataset-and-provenance)
[![Core genes](https://img.shields.io/badge/Sensitivity_supported_genes-285-008080)](#current-results)
[![Study DOI](https://img.shields.io/badge/Source_DOI-10.1038%2Fs44318--026--00744--w-blue)](https://doi.org/10.1038/s44318-026-00744-w)
<img width="1280" height="353" alt="image" src="https://github.com/user-attachments/assets/754195db-ecd0-425c-8d02-ccbccd579ab7" />



**Integrating T2D genetics and human islet functional genomics to prioritize beta-cell regulatory mechanisms.**

This repository implements an audited donor-level analysis of public human beta-cell transcriptomes. The broader research plan connects T2D susceptibility variants to regulatory elements, genes, and beta-cell responses through GWAS, islet eQTL colocalization, chromatin annotation, and independent validation.

**Current scope:** data audit, count validation, donor pseudobulk, exploratory QC, DESeq2–edgeR comparison, and model-sensitivity analysis are complete. Genetic integration and independent replication remain planned. Findings are observational associations; they do not establish genetic causality.

**Maintainer:** Md. Tariqul Islam · [GitHub](https://github.com/mtariqi) · [ORCID](https://orcid.org/0009-0009-6545-8040)  
**Documentation updated:** October 5, 2026. Analysis milestone: October 4, 2026, America/Toronto.

The badges above are descriptive project labels, not automated checks. No CI, coverage, license, DOI for this repository, container, or workflow-engine certification is claimed.


## Contents

- [Research objective and rationale](#research-objective-and-rationale)
- [Research gap](#research-gap)
- [Hypotheses](#hypotheses)
- [Research questions](#research-questions)
- [Research plan](#research-plan)
- [Dataset and provenance](#dataset-and-provenance)
- [Statistical design](#statistical-design)
- [Quality control and UMAP](#quality-control-and-umap)
- [Current results](#current-results)
- [Tools and methods](#tools-and-methods)
- [Repository structure](#repository-structure)
- [Running the current pipeline](#running-the-current-pipeline)
- [Interpretation and limitations](#interpretation-and-limitations)
- [Reproducibility and data policy](#reproducibility-and-data-policy)
- [References](#references)

## Research objective and rationale

The overall objective is to identify plausible beta-cell regulatory mechanisms associated with T2D susceptibility by integrating genetic association and functional genomic evidence.

The immediate objective is to identify donor-level beta-cell expression changes associated with T2D that remain consistent across two count-based statistical methods, BMI adjustment, and sequencing-chemistry restriction.

The unit of biological replication is the **donor**. Cells and repeated libraries from the same donor do not constitute independent biological replicates. Assigned beta-cell counts are aggregated across libraries by Islet.

This project focuses on T2D. T1D HLA and autoimmune analyses, gene–environment interactions, clinical risk prediction, and cross-ancestry comparisons have not been implemented.

## Research gap

GWAS associations frequently identify loci containing several possible effector genes and regulatory elements. Disease-associated expression can reflect consequences of disease or donor differences rather than an inherited causal mechanism. A simple overlap between differential-expression and GWAS/eQTL gene lists therefore leaves important mechanistic uncertainty.

The source GSE221156 study already investigates beta-cell expression defects and genetic/regulatory interpretation. Repeating its differential-expression analysis alone would not establish a novel contribution. The proposed contribution must be demonstrated through a source-study comparison, updated signal-aware genetic integration, explicit technical and model sensitivities, and independent replication.

The current analysis supplies a reproducible expression-evidence layer for that broader investigation.

## Hypotheses

| ID | Working hypothesis | Evaluation | Status |
|---|---|---|---|
| H1 | A subset of donor-level beta-cell transcripts differs between T2D and ND after adjusting for measured covariates. | Adjusted DESeq2 and edgeR tests; BH FDR below 0.05. | Evaluated; associated genes identified. |
| H2 | A subset retains effect direction and statistical support across methods, BMI adjustment, and V3 restriction. | Intersect per-analysis consensus sets; compare effects on common genes. | Evaluated; 285 sensitivity-supported genes; 279 strictly supported in leave-one-donor-out refits (270 significant in every refit). |
| H3 | Technical factors associated with chemistry contribute to donor expression structure. | PCA/UMAP annotated by chemistry; adjusted modeling and V3 sensitivity. | Chemistry-associated structure observed. Disease × chemistry interaction test done: 283 of 285 core genes consistent across chemistries; proteasome module decrease concentrated in V2 (interaction p = 0.049). |
| H4 | Some expression candidates will also have compatible T2D GWAS, islet eQTL, and chromatin evidence. | Harmonized, signal-aware colocalization and regulatory annotation. | Planned; not part of the current manuscript. |
| H5 | Prioritized candidates will reproduce in an independent cohort and respond to targeted perturbation in a suitable model. | Independent expression replication and controlled functional experiments. | Pre-registered replication in GSE153855 and E-MTAB-5061 (Smart-seq2, tags `prereg-validation-v1` and `-v2`): pooled random-effects analysis replicated the potassium-channel set and, at the margin, the core signature (I² = 79%); proteasome module not replicated. Functional experiments planned. |
| H6 | Robust T2D-associated genes converge on coherent biological pathways. | Ranked fgsea (Reactome, GO BP) across three models; redundancy clustering; leading-edge analysis. | Evaluated; excitability and potassium-channel genes higher in T2D; a 29-gene proteasome–ubiquitin module lower in T2D (small, concentrated in V2 chemistry; did not replicate). |

For H1, the per-gene statistical null is no disease-associated expression effect conditional on the fitted covariates. H2 is a descriptive sensitivity assessment, not a separate calibrated joint hypothesis test. H6 is exploratory; the proteasome–ubiquitin module was selected from discovery data and is tested prospectively under H5. H4 remains untested.

## Research questions

1. Which beta-cell genes differ between T2D and ND at the donor level?
2. How closely do DESeq2 and edgeR agree in effect estimates, significance, and direction?
3. Which associations remain supported following BMI adjustment and restriction to V3 libraries?
4. How do chemistry, recovered beta-cell numbers, and count depth relate to exploratory donor structure?
5. Are candidate changes broadly distributed across donors or strongly influenced by individual donors?
6. Which candidates share a genetic association signal between T2D GWAS and islet expression QTLs?
7. Do islet regulatory elements and cell-type evidence provide a compatible mechanism?
8. Which findings reproduce independently and merit experimental testing?

Prediabetes contrasts and pathway interpretation are planned secondary analyses.

## Research plan

| Phase | Activities and deliverables | Status |
|---|---|---|
| 1. Source audit | Inventory GEO files; resolve pooled libraries; record donor groups, asset URLs, and source hashes. | Complete |
| 2. Matrix validation | Inspect H5AD; validate identifiers, mappings, finite values, integer counts, and provenance. | Complete |
| 3. Donor pseudobulk | Sum assigned beta cells across donor libraries; check count conservation and donor coverage. | Complete |
| 4. Design and QC | Check full-rank models; inspect PCA, count metrics, clustering, and repeated UMAP. | Complete |
| 5. Differential expression | Fit DESeq2 and robust edgeR QL models with matched within-analysis gene filters. | Complete |
| 6. Robustness | Compare methods and sensitivity models; preserve complete and core gene tables. | Complete |
| 7. Candidate influence | Produce donor plots and leave-one-donor-out model refits. | Descriptive summaries reviewed; plots/refits pending |
| 8. Pathways and PD | Ranked pathway analysis with versioned gene sets; explicitly defined PD contrasts. | Planned |
| 9. Genetics | Audit GWAS/eQTL sources; harmonize builds/alleles; assess LD and multiple signals. | Planned |
| 10. Colocalization/chromatin | Signal-aware colocalization, prior sensitivity, and compatible regulatory-element evidence. | Planned |
| 11. Replication/validation | Independent donor cohort, candidate prioritization, and functional follow-up. | Planned |

Future functional work should specify the beta-cell model, independent biological replicates, perturbation and negative controls, molecular readouts, and functional endpoints before experimentation. No wet-laboratory experiment was performed in this computational phase.

## Dataset and provenance

| Property | Audited value |
|---|---|
| GEO accession | [GSE221156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156) |
| Source study | [Bandesh et al., 2026](https://doi.org/10.1038/s44318-026-00744-w) |
| Donors | 48: ND 17, PD 14, T2D 17 |
| GEO libraries | 54: 42 single-donor and 12 pooled |
| Listed GEO supplementary files | 162 |
| Beta-cell object | 99,029 cells × 26,936 genes |
| CELLxGENE schema | 7.1.0 |
| Count input | Author-curated unnormalized counts from `raw/X` |
| Output matrix | 26,936 genes × 48 donors |
| Conserved total counts | 2,081,765,364 |
| Minimum recovered beta cells | 582 per donor |
| Proposed 50-cell threshold | No donor below threshold |

**Selected H5AD:** [CELLxGENE beta-cell asset](https://datasets.cellxgene.cziscience.com/faacd98e-af54-459d-8a7a-74427766bc81.h5ad)  
**Collection:** [Human pancreatic islet transcriptomes](https://cellxgene.cziscience.com/collections/58e85c2f-d52e-4c19-8393-b854b84d516e)

Input SHA256:

```text
c0d7e00dc1bf91fb7682413730f7b6cec43dfb5a5d977c7bc3b517197adda39a
```

`X` contains fractional transformed expression. `raw/X` passes finite, nonnegative, integer-value checks despite float32 storage. The upstream correction history has not been independently reconstructed; the name `raw` does not establish untouched original counts.

## Statistical design

| Analysis | Donors | Formula | Residual df |
|---|---:|---|---:|
| Primary T2D vs ND | 34 | `~ chemistry + sex + age_z + disease_status` | 29 |
| BMI sensitivity | 34 | `~ chemistry + sex + age_z + BMI_z + disease_status` | 28 |
| V3-only sensitivity | 23 | `~ sex + age_z + disease_status` | 19 |

All designs were full rank. Age and BMI are standardized within each DE subset; ND is the reference. Positive log2 fold change means higher modeled expression in T2D.

- Apply `edgeR::filterByExpr` with each analysis design.
- Use the same retained genes for DESeq2 and edgeR within that analysis.
- DESeq2: size-factor normalization, dispersion estimation, Wald tests; disable automatic count replacement. Default Cook’s handling remains.
- Save native DESeq2 results and a comparison result with independent filtering disabled.
- edgeR: TMM normalization, robust dispersion estimation, robust quasi-likelihood fitting, and disease-coefficient testing.
- Use each method’s BH-adjusted p-value below 0.05. No minimum fold-change threshold or LFC shrinkage was applied.
- Define the core as significant in both methods in all three analyses, with consistent directions across all six fits.

Filters may produce different tested-gene sets across analyses. BMI adjustment changes the conditional estimand; V3 restriction also reduces sample size. Neither sensitivity analysis constitutes independent replication.

## Quality control and UMAP

QC includes all 48 donors, including PD. The exploratory count filter (at least 10 counts in at least 10 donors) retained 16,829 genes.

DESeq2 VST with `blind=TRUE` supports PCA and donor distances. PCA uses the top 500 variable genes; PC1 and PC2 explain 21.53% and 15.19%. A median-distance threshold of median plus three MADs flagged no donor. No donors were excluded.

UMAP inputs are the first ten PCs from the top 500 variable VST genes, recomputed separately for all donors and the 37 V3 donors. The V3 UMAP subset includes PD and differs from the 23-donor DE subset.

| Setting | Values |
|---|---|
| Neighbors | 5, 10, 15 |
| Seeds | 11, 42, 2026 |
| Minimum distance | 0.3 |
| Metric | Euclidean |
| Threads | One for neighbor processing and optimization |
| Embeddings | 18 |

Visual review shows chemistry-associated grouping and overlapping disease groups. Some five-neighbor groups change at larger neighborhood settings. UMAP is exploratory; the plots do not supply a differential-expression test, classifier, or quantified stability result.

## Current results

Milestone results from the October 4, 2026 analysis:

| Analysis | Tested genes | DESeq2 FDR <0.05 | edgeR FDR <0.05 | Both methods, same direction | LFC Spearman |
|---|---:|---:|---:|---:|---:|
| Primary | 17,352 | 1,144 | 983 | 913 | 0.9871 |
| BMI sensitivity | 17,352 | 892 | 691 | 656 | 0.9838 |
| V3 only | 17,346 | 662 | 378 | 364 | 0.9845 |

- Primary consensus overlaps BMI consensus in **611 genes** and V3 consensus in **322 genes**.
- **285 genes** are supported across both methods and all three analyses, with consistent effect directions.
- The full primary results remain the main analysis; the core is an additional sensitivity-supported subset.
- `maximum_FDR_all_six` is a descriptive sorting field, not a combined FDR.
- Preliminary donor summaries show concordant median decreases for PHLDB2 and MPP1 within both chemistries; low control abundance warrants further review for some strongly increasing genes.
- Genetic causality, independent replication, formal influence assessment, and functional validation remain unestablished.

See [core genes](metadata/GSE221156_expression/core_supported_genes.tsv), [robustness counts](metadata/GSE221156_expression/robustness_summary.tsv), and [audit documentation](GSE221156_audit/docs/GSE221156_AUDIT.md).

## Tools and methods

| Component | Tools |
|---|---|
| Version control and execution | Git, GitHub CLI, Bash |
| H5AD audit and sparse aggregation | Python, h5py, NumPy |
| Differential expression | R, Bioconductor, DESeq2, edgeR, statmod |
| Exploratory embeddings | Base R PCA; uwot UMAP |
| Independent table checks | Python pandas, NumPy, SciPy |
| Visual outputs | Base R PDF graphics |

Exact package versions are recorded in generated `sessionInfo.txt` files. No environment lockfile, CI workflow, Docker image, Nextflow pipeline, or Snakemake workflow is currently implemented.

## Repository structure

### Tracked files

The following inventory reflects the verified repository snapshot before this README update. Paths are relative to the repository root.

| Path | Purpose |
|---|---|
| `README.md` | Research rationale, status, methods, execution, and navigation |
| `.gitignore` | Excludes large data, generated results, credentials, and working files |
| `config/analysis.yaml` | Initial study configuration; current R scripts contain their own operational settings |
| `metadata/data_manifest.tsv` | Initial manifest; detailed GSE records are in the audit directory |
| `GSE221156_audit/docs/GSE221156_AUDIT.md` | Dataset audit and provenance notes |
| `GSE221156_audit/source_records/series.txt` | Archived GEO series record |
| `GSE221156_audit/source_records/samples.txt` | Archived GEO sample records |
| `GSE221156_audit/source_records/cellxgene.json` | Archived CELLxGENE asset metadata |
| `GSE221156_audit/source_records/zenodo.json` | Archived Zenodo record |
| `GSE221156_audit/scripts/audit_gse221156.py` | Parse records into donor/library/file inventories |
| `GSE221156_audit/scripts/inspect_beta_h5ad.py` | Inspect H5AD metadata and matrix headers |
| `GSE221156_audit/metadata/GSE221156_audit_summary.json` | Cohort and library summary with source hashes |
| `GSE221156_audit/metadata/GSE221156_cellxgene_assets.tsv` | CELLxGENE download asset inventory |
| `GSE221156_audit/metadata/GSE221156_donors.tsv` | Donor inventory |
| `GSE221156_audit/metadata/GSE221156_libraries.tsv` | GEO library inventory |
| `GSE221156_audit/metadata/GSE221156_library_donor_links.tsv` | Library-to-donor links, including pooled libraries |
| `GSE221156_audit/metadata/GSE221156_geo_file_inventory.tsv` | GEO supplementary-file inventory |
| `GSE221156_audit/metadata/GSE221156_matrix_header_checks.json` | Matrix header inspection record |
| `scripts/validate_beta_h5ad.py` | Validate metadata mapping and numerical count properties |
| `scripts/build_beta_pseudobulk.py` | Aggregate beta-cell counts by donor and check conservation |
| `scripts/check_beta_design.R` | Covariate balance and model-rank checks |
| `scripts/qc_beta_pseudobulk.R` | VST, PCA, donor QC, clustering, and review flags |
| `scripts/umap_beta_donors.R` | Repeated donor UMAP across seeds and neighborhoods |
| `scripts/compare_beta_DE_methods.R` | DESeq2 and robust edgeR QL fits for three designs |
| `scripts/summarize_beta_robustness.R` | Cross-model support and core-gene tables |
| `metadata/GSE221156_beta_validation.json` | H5AD validation summary |
| `metadata/GSE221156_beta_library_donor_cells.tsv` | Cell counts by library/donor |
| `metadata/GSE221156_beta_pseudobulk_samples.tsv` | Donor metadata aligned to pseudobulk columns |
| `metadata/GSE221156_beta_pseudobulk_summary.json` | Aggregation provenance, checksums, and count checks |
| `metadata/GSE221156_beta_gene_annotations.tsv` | Gene identifiers and symbols |
| `metadata/GSE221156_beta_design_checks.tsv` | Rank and residual degrees of freedom |
| `metadata/GSE221156_expression/core_supported_genes.tsv` | 285-gene sensitivity-supported subset |
| `metadata/GSE221156_expression/robustness_summary.tsv` | Per-analysis consensus counts |

### Local input and generated outputs

These paths are local workflow artifacts. `data/raw/`, `data/processed/`, and `results/` are ignored by Git and are not included when cloning.

| Directory/path | Contents |
|---|---|
| `data/raw/GSE221156/beta_cells.h5ad` | Downloaded beta-cell source object |
| `data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz` | Donor count matrix |
| `results/qc/GSE221156/` | QC PDF, PCA_scores.tsv, donor_QC.tsv, beta_QC_vst.rds, sessionInfo.txt |
| `results/qc/GSE221156/` | donor_UMAP_sensitivity.pdf, donor_UMAP_coordinates.tsv, UMAP_sessionInfo.txt |
| `results/expression/GSE221156/` | For each of primary, BMI_sensitivity, V3_only: *_tested_genes.tsv, *_samples.tsv, *_DESeq2_edgeR.tsv |
| `results/expression/GSE221156/` | method_comparison_summary.tsv, DESeq2_edgeR_effect_comparison.pdf, sessionInfo.txt |
| `results/expression/GSE221156/robustness/` | all_genes_robustness.tsv, core_supported_genes.tsv, robustness_summary.tsv, sessionInfo.txt |
| `logs/` | Execution logs written using tee; not present in the verified tracked inventory |
| `.venv-audit/` | Local Python environment; do not commit it |

The current ignore file lists `.venv/` but does not specifically list `.venv-audit/` or `logs/`. Use explicit staging and review `git status` before committing.

### Planned extensions

The following are proposed future organization, not existing implemented modules.

| Proposed directory | Intended role |
|---|---|
| `docs/` | Research plan, decisions, dated reports, manuscript methods |
| `config/genetics/` | GWAS/eQTL definitions, genome builds, LD references, priors |
| `scripts/genetics/` | Harmonization, fine-mapping/colocalization, sensitivity checks |
| `scripts/pathways/` | Versioned ranked and overrepresentation analyses |
| `metadata/genetics/` | Source manifests, schemas, download terms, checksums |
| `results/genetics/` | Local genetic integration outputs |
| `results/validation/` | Independent replication and influence analyses |
| `tests/` | Meaningful checks for future harmonization and mapping logic |

## Running the current pipeline

### 1. Clone and prepare dependencies

Access to the private repository is required.

```bash
git clone https://github.com/mtariqi/t2d-islet-regulatory-genomics.git
cd t2d-islet-regulatory-genomics
python3 -m venv .venv-audit
source .venv-audit/bin/activate
python -m pip install numpy h5py
mkdir -p data/raw/GSE221156 logs
```

Install R and Bioconductor packages in your R environment:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(c("DESeq2", "edgeR"), ask = FALSE, update = FALSE)
install.packages(c("uwot", "statmod"), repos = "https://cloud.r-project.org")
```

These commands install available packages; they do not recreate a pinned environment. Preserve the resulting session information.

### 2. Obtain and verify the source object

Download the linked CELLxGENE asset into `data/raw/GSE221156/beta_cells.h5ad`, following source access and reuse terms. Verify its checksum:

```bash
sha256sum data/raw/GSE221156/beta_cells.h5ad
```

A checksum mismatch requires version/provenance review before reproducing the reported milestone.

### 3. Audit and run in dependency order

```bash
set -euo pipefail
source .venv-audit/bin/activate
mkdir -p logs

python GSE221156_audit/scripts/audit_gse221156.py \
  --source-dir GSE221156_audit/source_records \
  --outdir GSE221156_audit/metadata

python GSE221156_audit/scripts/inspect_beta_h5ad.py \
  data/raw/GSE221156/beta_cells.h5ad

python scripts/validate_beta_h5ad.py
python scripts/build_beta_pseudobulk.py \
  2>&1 | tee logs/GSE221156_beta_pseudobulk.log
Rscript scripts/check_beta_design.R \
  2>&1 | tee logs/GSE221156_beta_design_checks.log
Rscript scripts/qc_beta_pseudobulk.R \
  2>&1 | tee logs/GSE221156_beta_expression_QC.log
Rscript scripts/umap_beta_donors.R \
  2>&1 | tee logs/GSE221156_beta_UMAP.log
Rscript scripts/compare_beta_DE_methods.R \
  2>&1 | tee logs/GSE221156_beta_DE_method_comparison.log
Rscript scripts/summarize_beta_robustness.R \
  2>&1 | tee logs/GSE221156_beta_robustness.log
```

Scripts currently use repository-relative paths and settings defined in their source. `config/analysis.yaml` is an initial planning record, not a universal runtime configuration read by all scripts.

## Interpretation and limitations

- These are secondary observational expression analyses, not evidence that a gene mediates inherited T2D risk.
- Method agreement reuses the same data; V3 and BMI comparisons reuse overlapping donors.
- Reduced V3 significance can reflect smaller sample size and changed donor composition.
- Chemistry adjustment does not establish that all residual technical confounding has been removed.
- Recovered beta-cell numbers are sampling quantities, not direct measurements of beta-cell abundance in vivo.
- Large LFCs with low control abundance require donor-level review and influence assessment.
- The original study must be cited, and novelty must be established against its published analyses.
- Self-reported race is not genetically inferred ancestry.
- No clinical prediction performance, wet-laboratory validation, or causal genetic conclusion is reported.

## Reproducibility and data policy

Track source records, permitted metadata, checksums, sample mappings, settings, seeds, software versions, and the Git commit used for each run. Preserve complete tested-gene results and analysis-specific sample lists.

Large data, restricted participant information, credentials, and analysis intermediates should remain outside version control. Public availability of a dataset does not automatically remove its reuse obligations. Respect the source dataset’s terms and cite the authors and curators.

There is currently no repository LICENSE, CITATION.cff, published release, or archived software DOI in the verified tree. Add these after choosing appropriate licensing and citation details; source-data terms remain separate.

## References

- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology, 15*, Article 550. https://doi.org/10.1186/s13059-014-0550-8
- Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics, 26*(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616
- McInnes, L., Healy, J., & Melville, J. (2018). *UMAP: Uniform manifold approximation and projection for dimension reduction* [Preprint]. arXiv. https://doi.org/10.48550/arXiv.1802.03426

Source study: Bandesh et al. (2026), *Deep single-cell decoding of human pancreatic islets reveals T2D β-cell gene expression defects*, *The EMBO Journal, 45*, 3978–4005. https://doi.org/10.1038/s44318-026-00744-w. The author list is abbreviated here; export the complete publisher citation for manuscript submission.

Implementation: [uwot documentation](https://jlmelville.github.io/uwot/) · [CELLxGENE schema 7.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md)
