# GSE221156 data audit

Audit date: 2026-10-05 UTC (2026-10-04 evening America/Toronto).

## Decision
Suitable for a beta-cell reanalysis, conditional on verifying count values and cell-to-donor assignments in the downloaded processed object. Metadata and file inventory audit is complete; expression payload audit is not complete. No new biological findings have been generated.

## Verified cohort structure
The downloaded GEO text records contain 54 libraries, not 54 donors. Parsing all sample records resolves 48 unique donors: 17 ND, 14 PD, and 17 T2D. There are 42 single-donor libraries and 12 pooled libraries. Library chemistry is V2 for 12 libraries and V3 for 42 libraries. Age, sex, BMI, HbA1c, and a source field named race are present for all donors after resolving paired metadata. Source race is not genetically inferred ancestry. Islet centers include missing labels and inconsistent spellings. Library IDs are identifiers, not independently verified sequencing batch assignments. Preserve diagnosed T2D labels even when HbA1c is below 6.5; do not reclassify diagnosed cases from HbA1c alone.

## Pooled donor libraries
- Islet47 and Islet48: MS19006, MS19007
- Islet57 and Islet58: MS19016, MS19017
- Islet59 and Islet60: MS19020, MS19021
- Islet70 and Islet71: MS19034, MS19035
- Islet84 and Islet85: MS19042, MS19043
- Islet118 and Islet119: MS20086, MS20087

Library membership is not a cell-to-donor assignment. Do not assign an entire pool to one donor, split counts equally, or treat pooled libraries as extra biological replicates. Four paired-only sets account for eight donors without single-donor libraries. A single-donor-only fallback changes cohort coverage and must be reported separately.

## GEO inventory
All 54 sample records list barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz: 162 sample-level files total. The series lists GSE221156_RAW.tar. RAW in the archive filename does not establish whether matrices are unfiltered, uncorrected integer counts. GEO lists HG38 and describes STARsolo, SoupX, genetic demultiplexing, doublet removal, QC, and Harmony. This description does not by itself establish the exact processing state of each deposited matrix. Inspect file contents before choosing analysis inputs. SRA links are captured in the library table. FASTQ downloads are not necessary for the first analysis stage.

## Processed resources
The paper links a CELLxGENE collection and Zenodo analytical resource. Collection API inventory retrieved successfully:
- All annotated cells: 245878 cells, H5AD 5868428962 bytes.
- Alpha cells: 74812 cells, H5AD 1611117039 bytes.
- Delta cells: 10770 cells, H5AD 226349221 bytes.
- Beta cells: 99029 cells, H5AD 2477777480 bytes.

The beta-only object is the preferred next payload for a focused beta-cell study. Downloaded API metadata establish assets and advertised sizes, not the contents of obs, raw, X, or layers. Donor annotations and count semantics remain unverified until payload inspection. A raw group name is not proof of raw integer counts.

Zenodo record 14656366 lists Rmarkdown.zip, 1183611003 bytes, MD5 9445803af9a707eab4fafae2917c7752. Its contents have not been downloaded or inspected. The paper describes Pipeline_html.Rmd and Pipeline_html.html. A 2025 archive may precede the final 2026 publication; check versions and methods before replication.

## Required next checks
1. Download the beta H5AD from the API-listed asset URL, recording download date, source and SHA256.
2. Inspect obs columns for donor identity, disease status, chemistry and library identity. Confirm all intended donors and explicitly check pooled-library assignments.
3. Inspect X, raw.X and layers. Identify a documented count matrix, inspect all stored values for finiteness, nonnegativity and integrality using chunks; fractional ambient-corrected values require an appropriate modeling decision rather than silent rounding.
4. Check sparse format, feature identifiers, mitochondrial feature symbols and uniqueness of cell IDs. Prefix barcodes by library when combining library matrices.
5. Confirm final QC state, remove ambiguous donor calls and doublets where required, and compare donor counts with the published dataset.
6. Aggregate counts by donor and cell type, not GSM. Inspect chemistry, isolation center and disease cross-tabulations before fixing model covariates. Do not include unique library IDs as ordinary donor-level batch covariates.
7. If no validated count layer or donor mapping exists, inspect the authors' Zenodo archive and supplementary data. Resolve the mapping before using pooled cells.

## Novelty review
The associated 2026 paper already reports beta-cell pseudobulk differential expression, genetic-expression integration, networks and experimental follow-up. The 2025 Zenodo description also includes eQTL overlap. Merely repeating differential expression and adding eQTL overlap is not a new discovery paper. Use the published results as a benchmark, then evaluate a distinct contribution such as updated signal-aware colocalization with matched LD and prior sensitivity, independent donor-cohort replication, and allele-specific validation. None of these is established as novel without a dedicated literature review.

## Files in this audit package
- scripts/audit_gse221156.py: Python standard-library parser for saved authoritative records; validates cohort counts, count-triplet inventory and repeated donor metadata consistency.
- metadata/GSE221156_donors.tsv: 48 donor rows.
- metadata/GSE221156_libraries.tsv: 54 library rows.
- metadata/GSE221156_library_donor_links.tsv: membership edges only.
- metadata/GSE221156_geo_file_inventory.tsv: 162 URLs; byte sizes and payload checksums remain blank until download.
- metadata/GSE221156_cellxgene_assets.tsv: advertised asset sizes and URLs.
- metadata/GSE221156_audit_summary.json: audit totals and source checksums.
- metadata/GSE221156_matrix_header_checks.json: representative matrix header probes, when successful; not whole-matrix validation.
- source_records: retrieved GEO records and resource API snapshots.

## Reproduce the saved-record audit
From the audit directory:

```bash
python3 scripts/audit_gse221156.py --source-dir source_records --outdir metadata
```

## Primary sources
GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156
GEO series text: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156&targ=self&form=text&view=full
GEO sample text: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156&targ=gsm&form=text&view=full
Publication: https://doi.org/10.1038/s44318-026-00744-w
Zenodo: https://zenodo.org/records/14656366
CELLxGENE: https://cellxgene.cziscience.com/collections/58e85c2f-d52e-4c19-8393-b854b84d516e

## Reference in APA format
Bandesh, K., Motakis, E., Nargund, S., Kursawe, R., Selvam, V., Ansarullah, Bhuiyan, R. M., Eryilmaz, G. N., Willett, A. M., White, J. K., Krishnan, S. N., Spracklen, C. N., Ucar, D., & Stitzel, M. L. (2026). Deep single-cell decoding of human pancreatic islets reveals T2D beta-cell gene expression defects. The EMBO Journal, 45, 3978–4005. https://doi.org/10.1038/s44318-026-00744-w

## Representative matrix header checks
Partial gzip stream reads succeeded for GSM6846488/MS17001 and GSM6846505/MS19006. Both declare Matrix Market coordinate integer format. The single-donor header reports 36601 features by 5381 barcodes and 14713536 nonzero entries. The pooled header reports 36601 features by 5901 barcodes and 16065554 nonzero entries. Only the header and first entries were inspected. These checks do not validate every value, identify final cells, or supply pooled donor assignments.
