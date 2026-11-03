# Independent human islet beta-cell replication plan

## Predefined hypothesis
The discovery-derived 29-gene proteasome/ubiquitin module has a coordinated negative T2D-versus-ND transcriptional shift in an independent human beta-cell donor cohort. This is an *untested* replication hypothesis, not a validated result.

## Candidate datasets and feasibility

1. **GSE153855 (first choice)** — 6 control and 5 T2D human islet donors; GEO lists `GSE153855_Expression_counts_HQ_allsamples.txt.gz` and `GSE153855_Cell_annotation.txt.gz`. Verify which cells are beta cells, donor IDs, per-donor beta-cell counts, and whether counts are genuinely raw integer counts before pseudobulk. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153855
2. **GSE81608 (orthogonal technology)** — human islet single-cell RNA-seq in ND and T2D; GEO processed expression is RPKM. Do **not** use DESeq2 on RPKM. Consider FASTQ reprocessing or an appropriately specified donor-level analysis of transformed expression. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81608
3. **E-MTAB-5061 (additional independent study)** — human pancreatic islet single-cell RNA-seq with ND/T2D comparisons; check raw files and donor independence. https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061

## Analysis protocol

1. Confirm no sample/donor overlap with GSE221156, study provenance, informed use, disease labels, beta-cell annotations and available covariates.
2. Create a donor-level metadata manifest with donor, disease, sex, age, BMI, platform, beta-cell counts, and library information; document missingness.
3. Quality-control beta-cell assignments, gene detection and cell counts per donor. Aggregate raw counts by donor when appropriate. Do not treat cells as independent replicates.
4. Harmonize Ensembl IDs and gene symbols; document detection of each of the 29 prespecified genes and the background universe. Avoid imputing absent genes as zero effect.
5. Estimate T2D-versus-ND donor-level effects with an appropriate design. Keep a simpler prespecified design if sample size cannot support multiple covariates. Preserve direction convention: negative = lower in T2D.
6. Evaluate the fixed 29-gene signature using a donor-aware gene-set test (or prespecified competitive/rank test); additionally report per-gene direction concordance and effect-size correlations with discovery. Account for within-set correlation and avoid treating 29 genes as independent replicates.
7. Test a preregistered proteasome gene set and a small number of nonredundant discovery pathways. Control FDR across the prespecified tests and report effect sizes, uncertainty, coverage, and negative findings.
8. If data permit, repeat after excluding low-beta-cell-count donors and with BMI adjustment; assess platform heterogeneity and donor influence.

## Proposed decision rules (prospective)

- Primary replication endpoint: prespecified 29-gene module shows negative aggregate T2D-associated direction and a statistically supported gene-set test in at least one truly independent donor cohort.
- Secondary: majority of detectable genes have negative effects; proteasome pathway NES is negative; donor-level effect-size agreement is reported with confidence intervals.
- A nonsignificant test does not prove absence of effect. Opposite-sign supported enrichment must be reported prominently. Do not claim independent validation until executed.

## Deliverables

`metadata/validation/<accession>_donors.tsv`, `results/validation/<accession>_beta_DE.tsv`, `results/validation/<accession>_proteasome29.tsv`, `results/validation/replication_summary.tsv`, provenance logs, session info, and a donor-aware methods note.
