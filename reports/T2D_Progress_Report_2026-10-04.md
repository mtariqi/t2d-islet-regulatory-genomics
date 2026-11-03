# Daily research progress report
**Reporting date:** October 4, 2026 (America/Toronto)  
**Researcher:** Md. Tariqul Islam  
**Project:** T2D islet regulatory genomics  
**Repository:** https://github.com/mtariqi/t2d-islet-regulatory-genomics  
**Verified repository snapshot:** fa33e5f07f08722f3df9c9a3394e985fed03a656  
**Scope:** Work completed in the current research session, from repository initialization through donor-level expression review. This is an interim computational progress report, not a final manuscript.

## 1. Research objective
The overall objective is to connect T2D susceptibility signals to plausible regulatory mechanisms and beta-cell responses by integrating genetics, islet functional genomics, and donor-level transcriptomics. The immediate objective completed today was to establish an auditable beta-cell expression analysis using GSE221156, compare DESeq2 with edgeR, and assess sensitivity to BMI adjustment and sequencing chemistry.

The working question is: Which beta-cell expression changes associated with T2D are consistent across statistical methods and model specifications, and which can subsequently be supported by independent genetic and regulatory evidence? The current disease contrast is T2D versus non-diabetic (ND) donors. Prediabetes (PD) donors are retained in exploratory QC but have not yet undergone formal differential testing.

## 2. Research gap and intended contribution
The Research Topic calls for mechanistic interpretation of diabetes-associated genetic variation. A disease-associated transcript is not necessarily a causal risk gene: expression can reflect disease consequences, medication, cellular state, technical processing, or other donor differences. GWAS loci also often contain several plausible genes and regulatory elements. Integrating evidence therefore requires more than overlapping gene lists.

Our proposed contribution is an explicitly audited chain of evidence: variant association, signal-aware islet eQTL colocalization, compatible chromatin regulation, reproducible beta-cell expression, and eventual functional validation. Chemistry imbalance and donor-level replication are important methodological issues in the present dataset.

The original GSE221156 publication already investigates beta-cell expression defects and genetic/regulatory interpretation. Repeating differential expression or a simple eQTL overlap would not establish novelty. The new contribution remains provisional and must be defined against the original study using updated genetic analyses, explicit sensitivity assessment, and independent validation. Bandesh et al. (2026) is the source study; see the references and its complete methods before drafting a novelty claim.

## 3. Research strategy and experimental setup
This phase is a secondary computational analysis of an existing human donor study. No new human samples, perturbations, genotyping, or wet-laboratory experiments were performed today.

1. Audit GEO, CELLxGENE, and associated study records.
2. Establish library-to-donor mappings, including pooled libraries.
3. Inspect H5AD matrix encodings and identify the appropriate count matrix.
4. Sum assigned beta cells by donor across libraries; treat donors as biological replicates.
5. Check covariate balance and design-matrix estimability.
6. Perform PCA, clustering, count-depth QC, and parameter-sensitive donor UMAP.
7. Fit primary DESeq2 and robust edgeR quasi-likelihood models on identical retained genes within each analysis.
8. Compare effects, significance overlap, and model sensitivity.
9. Review candidate expression at donor level before genetic prioritization.
10. Subsequently audit GWAS/eQTL sources, harmonize alleles and genome builds, perform colocalization, examine chromatin evidence, and seek independent replication.

Future laboratory validation could examine selected regulatory elements or genes by perturbation in an appropriate beta-cell model, followed by expression and beta-cell functional readouts. Candidate selection, biological replication, controls, and power calculations remain to be specified after genetic prioritization; no validated laboratory protocol is claimed here.

## 4. Dataset audit and provenance
| Item | Verified result |
|---|---:|
| GEO libraries | 54 |
| Unique donors | 48 |
| ND / PD / T2D donors | 17 / 14 / 17 |
| Single-donor / pooled libraries | 42 / 12 |
| Listed GEO supplementary files | 162 |
| V2 / V3 libraries | 12 / 42 |
| Curated beta cells | 99,029 |
| Genes in beta-cell object | 26,936 |
| Minimum beta cells per donor | 582 |

Library counts and donor counts are different quantities. Pooled library pairs were resolved using the curated assignments and audited mappings. Donors were not multiplied into independent replicates by counting repeated libraries.

The selected object is “Reintegrated transcriptomes of Beta pancreatic islet cells,” CELLxGENE schema 7.1.0, distributed at:
https://datasets.cellxgene.cziscience.com/faacd98e-af54-459d-8a7a-74427766bc81.h5ad

Input SHA256:
`c0d7e00dc1bf91fb7682413730f7b6cec43dfb5a5d977c7bc3b517197adda39a`

Both X and raw/X have shape 99,029 × 26,936 and 337,547,507 stored CSR values. X contains fractional transformed values and is unsuitable as the input to count-based differential testing. raw/X is finite, nonnegative, and integer-valued despite float32 storage; its stored-value sum is 2,081,765,364. Its use is documented as author-curated unnormalized counts. The upstream correction history was not independently reconstructed. The H5AD key “raw” alone does not prove an untouched experimental count matrix.

Metadata checks found no duplicate cell identifiers, no duplicate feature identifiers, no multiple-donor mappings per Islet, no multiple-Islet mappings per donor ID, no inconsistent disease assignments, and no missing expected Islets.

## 5. Pseudobulk generation
Counts from raw/X were summed by Islet across assigned beta cells and libraries using chunked sparse-matrix reads. Accumulation used integer counts, and the final matrix contains 26,936 genes × 48 donors.

The output is `data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz`. Its total exactly matches the input: 2,081,765,364. All output values are nonnegative integers. The proposed 50-cell threshold excluded no donor; the least represented donor had 582 beta cells. Count conservation was also independently checked against the uploaded matrix and sample metadata.

Beta-cell recovery ranged from 582 to 5,294 cells per donor, with median 2,126. Differences in recovered cells cannot be interpreted directly as differences in beta-cell abundance in vivo.

## 6. Design checks
| Disease | V2 donors | V3 donors | Female | Male |
|---|---:|---:|---:|---:|
| ND | 4 | 13 | 5 | 12 |
| PD | 1 | 13 | 4 | 10 |
| T2D | 7 | 10 | 7 | 10 |

| Analysis | Formula used in differential testing | Donors | Parameters | Residual df |
|---|---|---:|---:|---:|
| Primary | ~ chemistry + sex + age_z + disease_status | 34 | 5 | 29 |
| BMI sensitivity | ~ chemistry + sex + age_z + BMI_z + disease_status | 34 | 6 | 28 |
| V3 only | ~ sex + age_z + disease_status | 23 | 4 | 19 |

All designs were full rank. Age and BMI were centered and scaled within each differential-analysis subset. ND is the disease reference, so positive log2 fold change indicates higher expression in T2D. Full rank establishes estimability; it does not establish absence of confounding.

BMI adjustment is a sensitivity analysis because it can change the estimand and potentially adjust for a variable related to disease pathways. HbA1c was not included as an adjustment variable in these disease-status models. PD was excluded from these three DE contrasts.

## 7. Tools and methods
| Tool | Role |
|---|---|
| Bash, Git, GitHub CLI | Repository setup, scripted execution, logging, version control |
| Python, h5py, NumPy | H5AD inspection, matrix validation, sparse chunk aggregation |
| R, Bioconductor, DESeq2 | Count modeling, size factors, dispersion estimation, Wald testing, VST QC |
| edgeR, statmod | TMM normalization, robust dispersion estimation, robust quasi-likelihood fitting/testing |
| uwot | Repeated donor UMAP embeddings |
| Base R | Model-rank checks, PCA, clustering, tables, PDF plots |
| Python pandas, NumPy, SciPy | Independent uploaded-table alignment, overlap/effect checks, descriptive CPM review |
| Poppler and image inspection | Visual review of uploaded PDFs |

The user ran the pipeline on Ubuntu in ~/t2d-islet-regulatory-genomics, using the Python .venv-audit environment and R. Reported installed packages include h5py 3.16.0 and NumPy 2.5.3; PDF producer metadata reports R 4.6.1. Exact versions of other packages should be taken from the generated sessionInfo files rather than inferred. No cloud or university HPC execution occurred in this session.

QC retained genes with at least 10 counts in at least 10 of the 48 donors. This exploratory filter retained 16,829 genes. DESeq2 variance stabilization used blind=TRUE; PCA used the top 500 variable genes without feature scaling. Donor distances used all QC-retained VST genes. A donor was flagged for review if its median distance exceeded the across-donor median by three MADs. This heuristic flagged no donor and is not a guarantee that all influential observations are absent.

Differential testing used edgeR filterByExpr with the analysis-specific model matrix. Each method used the same retained genes within each analysis, but the tested universe can differ across analyses. DESeq2 automatic count replacement was disabled; default Cook’s handling remained. Native DESeq2 adjusted results were saved alongside results with independent filtering disabled for method comparison. No log-fold-change shrinkage or minimum effect threshold was applied in these results.

edgeR used TMM, estimateDisp(robust=TRUE), glmQLFit(robust=TRUE), and glmQLFTest for the disease coefficient. Significance was defined by each method’s adjusted p-value below 0.05. Agreement statistics are descriptive and do not define a new calibrated joint FDR.

## 8. QC and UMAP findings
PCA explained 21.53% and 15.19% of variance on PC1 and PC2. Chemistry-associated structure was evident, supporting chemistry adjustment and the V3-only sensitivity analysis. No donors were excluded.

UMAP used the first ten PCs from the top 500 variable VST genes, with variable genes and PCA recomputed within each subset. The two subsets were all 48 donors and all 37 V3 donors, including PD. This differs from the 23-donor V3 DE subset, which includes ND and T2D only.

UMAP settings: neighbors 5, 10, and 15; seeds 11, 42, and 2026; min_dist 0.3; Euclidean distance; one thread for both neighbor processing and stochastic optimization. Eighteen embeddings were saved. V2 donors remained grouped apart in the all-donor views, while disease groups overlapped. V3-only disconnected groups at five neighbors were sensitive to neighborhood size and did not remain comparably disconnected at larger settings. These are visual observations, not quantitative embedding-stability tests or disease classification results.

## 9. Differential-expression results
| Analysis | Tested genes | DESeq2 FDR<0.05 | edgeR FDR<0.05 | Both methods | Jaccard | LFC Spearman |
|---|---:|---:|---:|---:|---:|---:|
| Primary | 17,352 | 1,144 | 983 | 913 | 0.7521 | 0.9871 |
| BMI sensitivity | 17,352 | 892 | 691 | 656 | 0.7077 | 0.9838 |
| V3 only | 17,346 | 662 | 378 | 364 | 0.5385 | 0.9845 |

The native and comparison DESeq2 significant-gene counts were identical in these runs. All genes significant in both methods had matching effect directions within each analysis. Effect-direction agreement across all finite tested results was 98.06%, 97.85%, and 97.25%; p-value rank Spearman correlations were 0.9807, 0.9776, and 0.9779.

The primary comparison contains 231 genes significant only in DESeq2 and 70 only in edgeR. Method concordance supports internal consistency, but both methods analyze the same samples and related count models.

## 10. Cross-model robustness
The primary consensus overlaps the BMI consensus in 611 genes and the V3-only consensus in 322 genes. The intersection across all three consensus sets contains 285 genes. Every one retains a consistent direction across all six method-model combinations.

| Cross-model comparison | Common tested genes | DESeq2 LFC Spearman | edgeR LFC Spearman |
|---|---:|---:|---:|
| Primary vs BMI-adjusted | 17,352 | 0.9561 | 0.9535 |
| Primary vs V3 only | 17,063 | 0.8951 | 0.9049 |

The local robustness script independently reproduced counts of 913, 656, 364, and 285. The core table’s maximum_FDR_all_six is a descriptive sorting field, not a combined p-value or formally calibrated FDR. The 285 genes are a sensitivity-supported subset; they do not replace the full primary results or constitute independent replication.

## 11. Preliminary donor-level review
Twelve candidates ranked by the maximum of the six adjusted p-values were examined descriptively. The following illustrate cross-chemistry patterns; CPM is counts divided by total donor counts times one million, without TMM or covariate adjustment.

| Gene | Direction | V2 ND median CPM | V2 T2D median CPM | V3 ND median CPM | V3 T2D median CPM |
|---|---|---:|---:|---:|---:|
| PHLDB2 | Lower in T2D | 15.73 | 11.45 | 17.15 | 10.68 |
| MPP1 | Lower in T2D | 42.80 | 28.20 | 39.88 | 24.08 |
| TAFA4 | Lower in T2D | 1.13 | 0.46 | 1.03 | 0.30 |
| DACT2 | Lower in T2D | 1.49 | 0.43 | 1.20 | 0.22 |
| FAIM2 | Higher in T2D | 0.05 | 3.61 | 0.15 | 1.89 |
| TBX2 | Higher in T2D | 0.03 | 1.03 | 0.02 | 1.36 |

PHLDB2 and MPP1 were detected in every ND and T2D donor. TBX2, TBX2-AS1, ELFN1, and BARX1 warrant low-abundance review because their control expression is small. This review does not demonstrate a causal function, exclude residual cell-mixture effects, or establish leave-one-out stability. Formal candidate plots and refitting were not completed.

## 12. Reproducibility and repository progress
Repository initialization, metadata audits, pseudobulk generation, design checks, QC, UMAP, DE comparison, and robustness summarization are now represented in the verified repository snapshot. Nine exact Python/R scripts are reproduced in Appendix B. README and analysis.yaml are also included as context.

Known session milestones: da61092 (initialization), 2bfce40 (audit), f60b30a (validation), and 9482978 (pseudobulk/design checks). The current verified main snapshot is fa33e5f07f08722f3df9c9a3394e985fed03a656; subsequent script/summary commits are included in that snapshot.

Important generated outputs:
- metadata/GSE221156_beta_validation.json
- metadata/GSE221156_beta_pseudobulk_samples.tsv
- metadata/GSE221156_beta_pseudobulk_summary.json
- metadata/GSE221156_beta_gene_annotations.tsv
- metadata/GSE221156_beta_design_checks.tsv
- metadata/GSE221156_expression/core_supported_genes.tsv
- metadata/GSE221156_expression/robustness_summary.tsv
- results/qc/GSE221156/beta_pseudobulk_QC.pdf
- results/qc/GSE221156/donor_UMAP_sensitivity.pdf
- results/expression/GSE221156/{primary,BMI_sensitivity,V3_only}_DESeq2_edgeR.tsv
- results/expression/GSE221156/DESeq2_edgeR_effect_comparison.pdf
- results/expression/GSE221156/robustness/all_genes_robustness.tsv

Large data and generated results may remain outside Git; source data hashes and sessionInfo files support traceability. No newly generated analysis has been pushed by the report-writing assistant.

## 13. Limitations and remaining work
This is observational expression evidence from 48 donors, with 34 in the main contrast and 23 in V3-only testing. Smaller sample size and changed composition can reduce significance in the V3 analysis. Chemistry adjustment cannot automatically remove every technical or donor-level difference. The absence of a heuristic QC flag does not imply absence of influential donors.

Genotyping, GWAS fine-mapping, eQTL colocalization, chromatin linkage, ancestry analyses, independent replication, pathway enrichment, PD differential testing, predictive modeling, and functional validation have not been completed. Self-reported race must not be substituted for genetically inferred ancestry. No causal genetic mechanism or clinical biomarker is established.

Next priorities:
1. Generate candidate donor-expression plots with disease and chemistry annotations.
2. Run leave-one-donor-out influence assessment using the full modeling gene universe and documented normalization/filtering decisions.
3. Conduct primary ranked pathway enrichment with a documented gene-set version and multiple-testing control.
4. Audit T2D GWAS and islet eQTL datasets, then harmonize builds and alleles.
5. Perform signal-aware colocalization with appropriate LD and assess chromatin support.
6. Define independent replication and a defensible novelty claim before manuscript preparation.

## 14. Research summary
Today’s work established a reproducible donor-level beta-cell analysis from an audited public single-cell dataset. The count matrix conserved all 2,081,765,364 input counts. Primary DESeq2 and edgeR analyses identified 913 jointly significant genes; 285 retained significance and consistent direction after BMI adjustment and V3 restriction. PCA and UMAP revealed chemistry-associated structure, supporting the chosen sensitivity strategy. Preliminary donor summaries show interpretable cross-chemistry changes for selected candidates, while highlighting low-abundance genes for further review.

Suggested interim Results text:
“Donor-level beta-cell pseudobulk analysis identified 913 genes significant in both DESeq2 and edgeR in the primary T2D-versus-ND comparison. Of these, 285 remained significant in both methods following BMI adjustment and restriction to V3 libraries, with consistent effect directions across all analyses.”

## 15. References and source records
Method references (APA style):
- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology, 15*, Article 550. https://doi.org/10.1186/s13059-014-0550-8
- Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics, 26*(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616
- McInnes, L., Healy, J., & Melville, J. (2018). *UMAP: Uniform manifold approximation and projection for dimension reduction* [Preprint]. arXiv. https://doi.org/10.48550/arXiv.1802.03426

Source-study bibliographic record: Bandesh et al. (2026). Deep single-cell decoding of human pancreatic islets reveals T2D β-cell gene expression defects. *The EMBO Journal, 45*, 3978–4005. https://doi.org/10.1038/s44318-026-00744-w
(The author list is abbreviated here; export the complete publisher citation for an APA manuscript reference.)

Data and software records:
- GEO GSE221156: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156
- CELLxGENE collection: https://cellxgene.cziscience.com/collections/58e85c2f-d52e-4c19-8393-b854b84d516e
- Schema: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md
- uwot: https://jlmelville.github.io/uwot/
- Exact source snapshot: https://github.com/mtariqi/t2d-islet-regulatory-genomics/tree/fa33e5f07f08722f3df9c9a3394e985fed03a656

Primary evidence for numerical results consists of user execution logs, uploaded DE tables, uploaded sample/count files, and the reviewed QC/UMAP/effect-comparison PDFs. Repository source code was fetched at the fixed snapshot above. Numerical conclusions were independently checked from the supplied tables and count matrix.

## Appendix A. Execution commands
Run commands from the repository root. These reproduce the workflow assuming the audited source records and H5AD input have already been downloaded. Repository initialization is a historical step and should not be repeated in an existing repository.

```bash
# Historical initialization (already completed)
git add README.md .gitignore config/analysis.yaml metadata/data_manifest.tsv
git commit -m "Initialize T2D islet regulatory genomics study"
gh repo create mtariqi/t2d-islet-regulatory-genomics \
  --private --source=. --remote=origin --push \
  --description "Integrating T2D genetics and human islet functional genomics"

# Analysis environment and logs
cd ~/t2d-islet-regulatory-genomics
source .venv-audit/bin/activate
mkdir -p logs
set -o pipefail

python GSE221156_audit/scripts/audit_gse221156.py \
  --source-dir GSE221156_audit/source_records \
  --outdir GSE221156_audit/metadata

python GSE221156_audit/scripts/inspect_beta_h5ad.py \
  data/raw/GSE221156/beta_cells.h5ad

# Consult each script's argparse defaults before supplying optional paths.
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

Package installation used in this session:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc <- c("DESeq2", "edgeR")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
cran <- c("uwot", "statmod")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
```

## Appendix B. Exact repository source code
These are the final committed versions at the verified snapshot, including the age-scaling and PDF-title corrections. Transient earlier variants, terminal authentication prompts, and credentials are intentionally not duplicated.

### GSE221156_audit/scripts/audit_gse221156.py
```python
#!/usr/bin/env python3
"""Audit saved GEO records. Does not download expression data or infer cell donors."""
import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def write_tsv(path, rows, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def audit(source, out):
    out.mkdir(parents=True, exist_ok=True)
    text = (source / 'samples.txt').read_text()
    libraries, links, files, donors = [], [], [], {}
    stable_fields = {'disease_status_source': 'disease state', 'sex': 'sex',
                     'age': 'age', 'BMI': 'bmi', 'HbA1c': 'hba1c',
                     'reported_race_source': 'race'}
    for block in text.split('^SAMPLE = ')[1:]:
        lines = block.splitlines()
        gsm = lines[0].strip()
        values, char = {}, {}
        for line in lines[1:]:
            if ' = ' not in line:
                continue
            key, value = line.split(' = ', 1)
            values.setdefault(key, []).append(value)
            if key == '!Sample_characteristics_ch1':
                ck, cv = value.split(':', 1)
                char[ck.lower().strip()] = cv.strip()
        title = values['!Sample_title'][0]
        if not re.fullmatch(r'Islet\d+(?:/\d+)?(?:_MS\d+)?', title):
            raise ValueError(f'Unexpected donor title: {title}')
        donor_ids = ['Islet' + n for n in re.findall(r'\d+', title.split('_')[0])]
        pooled = len(donor_ids) > 1
        library = char['library id']
        supplementary = [v for k, vs in values.items()
                         if k.startswith('!Sample_supplementary_file') for v in vs]
        types = []
        for url in supplementary:
            kind = next((k for k in ('barcodes', 'features', 'matrix')
                         if ('_' + k + '.') in url), 'other')
            types.append(kind)
            files.append({'gsm': gsm, 'library_id': library, 'file_type': kind,
                          'file_name': url.rsplit('/', 1)[-1],
                          'url': url.replace('ftp://', 'https://'),
                          'byte_size': '', 'sha256': '',
                          'validation_status': 'listed_only_not_downloaded'})
        if sorted(types) != ['barcodes', 'features', 'matrix']:
            raise ValueError(f'Incomplete triplet: {gsm}')
        libraries.append({'gsm': gsm, 'library_id': library, 'title': title,
                          'donor_ids': ';'.join(donor_ids), 'pooled': str(pooled).lower(),
                          'chemistry': char['chemistry'],
                          'disease_status_source': char['disease state'],
                          'center_source': char.get('center', ''),
                          'assembly_source': ';'.join(v.removeprefix('Assembly: ')
                              for v in values.get('!Sample_data_processing', [])
                              if v.startswith('Assembly: ')),
                          'sra_relation': ';'.join(v for v in values.get('!Sample_relation', [])
                                                  if v.startswith('SRA: '))})
        for i, donor in enumerate(donor_ids):
            row = {'donor_id': donor}
            for field, key in stable_fields.items():
                parts = char[key].split('/')
                if len(parts) != len(donor_ids):
                    raise ValueError(f'Cannot resolve {key} for {gsm}')
                row[field] = parts[i]
            row['disease_status'] = {'NonDiabetic': 'ND', 'PreDiabetic': 'PD',
                                    'Type2Diabetic': 'T2D'}[row['disease_status_source']]
            if donor in donors and donors[donor] != row:
                raise ValueError(f'Conflicting donor metadata: {donor}')
            donors[donor] = row
            links.append({'gsm': gsm, 'library_id': library, 'donor_id': donor,
                          'mapping_level': 'library_membership_only',
                          'cell_assignment_status': 'requires_demultiplexing_map' if pooled
                          else 'single_donor_library_pending_cell_QC'})
    groups = Counter(d['disease_status'] for d in donors.values())
    if len(libraries) != 54 or len(donors) != 48 or groups != Counter(ND=17, PD=14, T2D=17):
        raise ValueError(f'Unexpected cohort: {len(libraries)}, {len(donors)}, {groups}')
    for name, rows in [('GSE221156_libraries.tsv', libraries),
                       ('GSE221156_donors.tsv', sorted(donors.values(), key=lambda d: int(d['donor_id'][5:]))),
                       ('GSE221156_library_donor_links.tsv', links),
                       ('GSE221156_geo_file_inventory.tsv', files)]:
        write_tsv(out / name, rows, list(rows[0]))
    assets = []
    collection = json.loads((source / 'cellxgene.json').read_text())
    for dataset in collection['datasets']:
        for asset in dataset.get('assets', []):
            assets.append({'dataset_id': dataset['dataset_id'], 'title': dataset['title'],
                           'cell_count': dataset['cell_count'], 'file_type': asset['filetype'],
                           'byte_size': asset['filesize'], 'url': asset['url'],
                           'payload_status': 'not_downloaded_or_inspected'})
    write_tsv(out / 'GSE221156_cellxgene_assets.tsv', assets, list(assets[0]))
    summary = {'geo_libraries': len(libraries), 'unique_donors': len(donors),
               'donor_groups': dict(groups), 'pooled_libraries': sum(l['pooled']=='true' for l in libraries),
               'single_donor_libraries': sum(l['pooled']=='false' for l in libraries),
               'listed_geo_files': len(files), 'chemistry_library_counts': dict(Counter(l['chemistry'] for l in libraries)),
               'expression_payloads_inspected': False,
               'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.iterdir()) if p.is_file()}}
    (out / 'GSE221156_audit_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    audit(args.source_dir, args.outdir)

```

### GSE221156_audit/scripts/inspect_beta_h5ad.py
```python
#!/usr/bin/env python3
"""Inspect H5AD metadata and matrix layout without loading expression into memory.
Requires h5py. Does not certify integer counts or infer donor assignments.
"""
import argparse
import json
import h5py


def scalar(value):
    if isinstance(value, bytes):
        return value.decode()
    if hasattr(value, 'tolist'):
        return scalar(value.tolist())
    if isinstance(value, (list, tuple)):
        return [scalar(v) for v in value]
    return value


def layout(node):
    result = {'encoding': scalar(node.attrs.get('encoding-type', ''))}
    if isinstance(node, h5py.Dataset):
        result.update(shape=list(node.shape), dtype=str(node.dtype))
    else:
        result.update(keys=list(node.keys()), shape=scalar(node.attrs.get('shape', [])))
        if 'data' in node:
            result['data_dtype'] = str(node['data'].dtype)
            result['stored_values'] = node['data'].shape[0]
    return result


p = argparse.ArgumentParser()
p.add_argument('input')
a = p.parse_args()
with h5py.File(a.input, 'r') as f:
    result = {'top_level_keys': list(f.keys()), 'obs_columns': list(f['obs'].keys()),
              'var_columns': list(f['var'].keys()), 'matrices': {}}
    for key in ('X', 'raw/X'):
        if key in f:
            result['matrices'][key] = layout(f[key])
    if 'layers' in f:
        for key in f['layers']:
            result['matrices']['layers/'+key] = layout(f['layers'][key])
    result['candidate_obs_fields'] = {}
    for key, node in f['obs'].items():
        if any(t in key.lower() for t in ('donor', 'sample', 'library', 'disease', 'chemistry', 'batch', 'cell_type')):
            if isinstance(node, h5py.Group) and 'categories' in node:
                result['candidate_obs_fields'][key] = {'categories': scalar(node['categories'][:])}
            elif isinstance(node, h5py.Dataset):
                result['candidate_obs_fields'][key] = {'first_10': scalar(node[:10])}
    print(json.dumps(result, indent=2))

```

### README.md
```markdown
# T2D Islet Regulatory Genomics

## Purpose
Integrate T2D genetics and human islet functional genomics
to prioritize beta-cell regulatory mechanisms.

## Status
Study design and repository initialization.
No research findings have been generated.

## Primary comparison
T2D versus non-diabetic donor-level beta-cell pseudobulk expression.

## Planned analyses
- Data and metadata audit
- Single-cell quality control and annotation
- Donor-level pseudobulk differential expression
- GWAS and eQTL harmonization
- Colocalization and sensitivity analysis
- Chromatin annotation
- Candidate prioritization
- Independent validation

## Data policy
Large data files and restricted participant information are excluded.
Metadata, provenance, checksums, and permitted derived outputs
are tracked when source terms allow.

## Reproducibility
Record software versions, configuration, random seeds,
input checksums, and Git commit for each analysis run.

## Interpretation
Expression changes and computational predictions do not alone
establish genetic causality.

```

### config/analysis.yaml
```yaml
project: t2d-islet-regulatory-genomics
status: planning
species: Homo_sapiens

expression:
  primary_dataset: GSE221156
  cell_type: beta
  primary_comparison: T2D_vs_ND
  fdr: 0.05
  proposed_min_cells_per_donor: 50

genetics:
  harmonized_genome_build: null
  coloc_h4_threshold: 0.80
  require_prior_sensitivity: true

reproducibility:
  seed: 20261004

```

### scripts/build_beta_pseudobulk.py
```python
#!/usr/bin/env python3
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE221156/beta_cells.h5ad"
OUT = ROOT / "data/processed/GSE221156"
META = ROOT / "metadata"
OUT.mkdir(parents=True, exist_ok=True)
META.mkdir(exist_ok=True)


def decode(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def column(group, name):
    node = group[name]
    if isinstance(node, h5py.Dataset):
        return [decode(x) for x in node[:]]
    if "categories" in node and "codes" in node:
        categories = [decode(x) for x in node["categories"][:]]
        return [
            categories[int(code)] if code >= 0 else ""
            for code in node["codes"][:]
        ]
    raise ValueError(f"Unsupported encoding: {name}")


def normalize_islet(value):
    value = value.strip()
    if value.startswith("Islet"):
        return value
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"Invalid Islet value: {value}")
    return f"Islet{int(number)}"


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


validation = json.loads(
    (META / "GSE221156_beta_validation.json").read_text()
)
if not validation["metadata_pass"]:
    raise ValueError("Previous metadata validation failed")
if not validation["matrices"]["raw/X"]["passes_numeric_count_checks"]:
    raise ValueError("Previous count validation failed")

with (
    ROOT / "GSE221156_audit/metadata/GSE221156_donors.tsv"
).open() as handle:
    donor_metadata = {
        row["donor_id"]: row
        for row in csv.DictReader(handle, delimiter="\t")
    }

with (
    ROOT / "GSE221156_audit/metadata/GSE221156_libraries.tsv"
).open() as handle:
    chemistry = {
        row["library_id"]: row["chemistry"]
        for row in csv.DictReader(handle, delimiter="\t")
    }

with (
    ROOT / "GSE221156_audit/metadata/GSE221156_library_donor_links.tsv"
).open() as handle:
    allowed = {
        (row["library_id"], row["donor_id"])
        for row in csv.DictReader(handle, delimiter="\t")
    }

disease_map = {
    "normal": "ND",
    "prediabetes syndrome": "PD",
    "type 2 diabetes mellitus": "T2D",
}

with h5py.File(INPUT, "r") as f:
    obs = f["obs"]
    islets = [normalize_islet(x) for x in column(obs, "Islet")]
    donor_ids = column(obs, "donor_id")
    libraries = column(obs, "LibraryID")
    diseases = column(obs, "disease")

    if set(islets) != set(donor_metadata):
        raise ValueError("Observed donors differ from the audited cohort")

    order = sorted(donor_metadata, key=lambda x: int(x[5:]))
    donor_position = {donor: i for i, donor in enumerate(order)}
    cell_donors = np.array(
        [donor_position[x] for x in islets], dtype=np.int64
    )

    donor_id_sets = defaultdict(set)
    library_sets = defaultdict(set)
    chemistry_sets = defaultdict(set)
    n_cells = Counter(islets)

    for islet, donor_id, library, disease in zip(
        islets, donor_ids, libraries, diseases
    ):
        base_library = library.split("_", 1)[0]
        if (base_library, islet) not in allowed:
            raise ValueError(f"Invalid membership: {library}, {islet}")
        if disease_map.get(disease) != donor_metadata[islet]["disease_status"]:
            raise ValueError(f"Disease mismatch: {islet}")
        donor_id_sets[islet].add(donor_id)
        library_sets[islet].add(base_library)
        chemistry_sets[islet].add(chemistry[base_library])

    if any(len(x) != 1 for x in donor_id_sets.values()):
        raise ValueError("Donor identifier mapping is inconsistent")

    var = f["raw/var"]
    index_name = decode(var.attrs.get("_index", "_index"))
    genes = column(var, index_name)
    symbols = (
        column(var, "feature_name")
        if "feature_name" in var else [""] * len(genes)
    )
    if len(set(genes)) != len(genes):
        raise ValueError("Duplicate gene identifiers")

    matrix = f["raw/X"]
    shape = tuple(int(x) for x in matrix.attrs["shape"])
    if shape != (len(islets), len(genes)):
        raise ValueError("Matrix and metadata dimensions differ")

    indptr = matrix["indptr"][:]
    if (
        len(indptr) != shape[0] + 1
        or indptr[0] != 0
        or indptr[-1] != len(matrix["data"])
        or np.any(indptr[1:] < indptr[:-1])
    ):
        raise ValueError("Invalid CSR row pointers")

    # Approximately 10 MB for the final donor-by-gene matrix.
    pseudobulk = np.zeros((len(order), len(genes)), dtype=np.int64)
    input_total = 0

    for start in range(0, shape[0], 1000):
        stop = min(start + 1000, shape[0])
        left, right = int(indptr[start]), int(indptr[stop])
        values = matrix["data"][left:right]
        indices = matrix["indices"][left:right]

        if (
            not np.isfinite(values).all()
            or np.any(values < 0)
            or np.any(values != np.floor(values))
            or np.any(indices < 0)
            or np.any(indices >= shape[1])
        ):
            raise ValueError(f"Invalid counts or indices near cell {start}")

        integer_values = values.astype(np.int64)
        entry_donors = np.repeat(
            cell_donors[start:stop],
            np.diff(indptr[start:stop + 1]),
        )
        np.add.at(pseudobulk, (entry_donors, indices), integer_values)
        input_total += int(integer_values.sum())

        if start % 10000 == 0:
            print(f"Aggregated {stop:,}/{shape[0]:,} cells", flush=True)

    if int(pseudobulk.sum()) != input_total:
        raise ValueError("Aggregation did not conserve counts")
    if input_total != int(
        validation["matrices"]["raw/X"]["sum_stored_values"]
    ):
        raise ValueError("Total differs from the validated input")

    provenance_fields = {}
    for key in ("schema_version", "schema_reference", "title", "citation"):
        if key in f["uns"] and isinstance(f["uns"][key], h5py.Dataset):
            provenance_fields[key] = decode(f["uns"][key][()])

counts_path = OUT / "beta_pseudobulk_counts.tsv.gz"
with gzip.open(counts_path, "wt", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["gene_id", *order])
    for j, gene in enumerate(genes):
        writer.writerow([gene, *pseudobulk[:, j].tolist()])

write_tsv(
    META / "GSE221156_beta_gene_annotations.tsv",
    [{"gene_id": g, "gene_symbol": s} for g, s in zip(genes, symbols)],
    ["gene_id", "gene_symbol"],
)

rows = []
for donor in order:
    row = dict(donor_metadata[donor])
    row.update(
        cellxgene_donor_id=next(iter(donor_id_sets[donor])),
        n_beta_cells=n_cells[donor],
        total_counts=int(pseudobulk[donor_position[donor]].sum()),
        library_ids=";".join(sorted(library_sets[donor])),
        chemistry=";".join(sorted(chemistry_sets[donor])),
        meets_proposed_50_cell_threshold=n_cells[donor] >= 50,
    )
    rows.append(row)

write_tsv(
    META / "GSE221156_beta_pseudobulk_samples.tsv",
    rows,
    list(rows[0]),
)

sha256 = hashlib.sha256()
with INPUT.open("rb") as handle:
    for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
        sha256.update(block)

summary = {
    "input_sha256": sha256.hexdigest(),
    "matrix_used": "raw/X",
    "input_description": "Author-curated unnormalized counts",
    "upstream_correction_history": "Not independently reconstructed",
    "aggregation": "Sum assigned beta cells by Islet across libraries",
    "n_cells": len(islets),
    "n_donors": len(order),
    "n_genes": len(genes),
    "donor_groups": dict(Counter(r["disease_status"] for r in rows)),
    "minimum_beta_cells_per_donor": min(n_cells.values()),
    "donors_below_proposed_50_cell_threshold": [
        donor for donor in order if n_cells[donor] < 50
    ],
    "total_counts": input_total,
    "count_conservation_pass": True,
    "nonnegative_integer_output": True,
    "counts_output": str(counts_path.relative_to(ROOT)),
    "object_provenance": provenance_fields,
}
(META / "GSE221156_beta_pseudobulk_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n"
)
print(json.dumps(summary, indent=2))

```

### scripts/check_beta_design.R
```r
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

```

### scripts/compare_beta_DE_methods.R
```r
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))

out <- "results/expression/GSE221156"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(
  "metadata/GSE221156_beta_pseudobulk_samples.tsv",
  stringsAsFactors = FALSE, check.names = FALSE
)
counts <- as.matrix(read.delim(
  gzfile("data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz"),
  row.names = 1, check.names = FALSE
))
annotations <- read.delim(
  "metadata/GSE221156_beta_gene_annotations.tsv",
  stringsAsFactors = FALSE
)

stopifnot(
  !anyDuplicated(samples$donor_id),
  !anyDuplicated(rownames(counts)),
  setequal(colnames(counts), samples$donor_id),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == floor(counts)),
  max(counts) <= .Machine$integer.max
)
samples <- samples[match(colnames(counts), samples$donor_id), ]
rownames(samples) <- samples$donor_id
stopifnot(all(colSums(counts) == samples$total_counts))
storage.mode(counts) <- "integer"

configs <- list(
  primary = list(
    selected = samples$disease_status %in% c("ND", "T2D"),
    formula = ~ chemistry + sex + age_z + disease_status
  ),
  BMI_sensitivity = list(
    selected = samples$disease_status %in% c("ND", "T2D"),
    formula = ~ chemistry + sex + age_z + BMI_z + disease_status
  ),
  V3_only = list(
    selected = samples$disease_status %in% c("ND", "T2D") &
      samples$chemistry == "V3",
    formula = ~ sex + age_z + disease_status
  )
)

summaries <- list()
pdf(file.path(out, "DESeq2_edgeR_effect_comparison.pdf"),
    width = 8, height = 7)

for (name in names(configs)) {
  config <- configs[[name]]
  m <- samples[config$selected, , drop = FALSE]
  m$disease_status <- factor(m$disease_status, levels = c("ND", "T2D"))
  m$chemistry <- factor(m$chemistry)
  m$sex <- factor(m$sex)
  m$age_z <- as.numeric(scale(as.numeric(m$age)))
  m$BMI_z <- as.numeric(scale(as.numeric(m$BMI)))

  design <- model.matrix(config$formula, data = m)
  stopifnot(
    !anyNA(design),
    nrow(design) == nrow(m),
    qr(design)$rank == ncol(design)
  )

  cts <- counts[, rownames(m), drop = FALSE]

  # Same expression filter for both methods within this analysis.
  keep <- edgeR::filterByExpr(
    edgeR::DGEList(counts = cts), design = design
  )
  cts <- cts[keep, , drop = FALSE]
  if (!nrow(cts)) stop("No genes retained")

  write.table(
    data.frame(gene_id = rownames(cts)),
    file.path(out, paste0(name, "_tested_genes.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    m, file.path(out, paste0(name, "_samples.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # DESeq2: primary method; disable automatic count replacement.
  dds <- DESeqDataSetFromMatrix(
    countData = cts, colData = m, design = config$formula
  )
  dds <- DESeq(dds, minReplicatesForReplace = Inf)

  native <- results(
    dds, contrast = c("disease_status", "T2D", "ND"), alpha = 0.05
  )
  comparison <- results(
    dds,
    contrast = c("disease_status", "T2D", "ND"),
    alpha = 0.05,
    independentFiltering = FALSE
  )

  # edgeR: TMM normalization and robust quasi-likelihood model.
  y <- edgeR::DGEList(counts = cts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  coefficient <- match("disease_statusT2D", colnames(design))
  stopifnot(!is.na(coefficient))
  test <- edgeR::glmQLFTest(fit, coef = coefficient)
  er <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  er <- er[rownames(cts), , drop = FALSE]

  tab <- data.frame(
    gene_id = rownames(cts),
    gene_symbol = annotations$gene_symbol[
      match(rownames(cts), annotations$gene_id)
    ],
    DESeq2_baseMean = native$baseMean,
    DESeq2_log2FC = native$log2FoldChange,
    DESeq2_pvalue = native$pvalue,
    DESeq2_FDR_native = native$padj,
    DESeq2_FDR_no_independent_filter = comparison$padj,
    edgeR_log2FC = er$logFC,
    edgeR_logCPM = er$logCPM,
    edgeR_pvalue = er$PValue,
    edgeR_FDR = er$FDR
  )

  valid <- is.finite(tab$DESeq2_pvalue) &
    is.finite(tab$edgeR_pvalue) &
    is.finite(tab$DESeq2_log2FC) &
    is.finite(tab$edgeR_log2FC)

  ds_sig <- !is.na(tab$DESeq2_FDR_no_independent_filter) &
    tab$DESeq2_FDR_no_independent_filter < 0.05
  er_sig <- !is.na(tab$edgeR_FDR) & tab$edgeR_FDR < 0.05
  both <- ds_sig & er_sig
  either <- ds_sig | er_sig

  tab$significant_both <- both
  tab$same_effect_direction <-
    sign(tab$DESeq2_log2FC) == sign(tab$edgeR_log2FC)

  write.table(
    tab, file.path(out, paste0(name, "_DESeq2_edgeR.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  summaries[[name]] <- data.frame(
    analysis = name,
    n_donors = nrow(m),
    genes_tested = nrow(tab),
    genes_with_finite_results_both = sum(valid),
    DESeq2_native_FDR05 = sum(
      !is.na(tab$DESeq2_FDR_native) & tab$DESeq2_FDR_native < 0.05
    ),
    DESeq2_comparison_FDR05 = sum(ds_sig),
    edgeR_FDR05 = sum(er_sig),
    overlap_FDR05 = sum(both),
    overlap_same_direction = sum(both & tab$same_effect_direction),
    significance_Jaccard = if (sum(either)) sum(both) / sum(either) else NA,
    log2FC_Spearman = cor(
      tab$DESeq2_log2FC[valid], tab$edgeR_log2FC[valid],
      method = "spearman"
    ),
    effect_direction_agreement = mean(tab$same_effect_direction[valid]),
    pvalue_rank_Spearman = cor(
      tab$DESeq2_pvalue[valid], tab$edgeR_pvalue[valid],
      method = "spearman"
    )
  )

  plot(
    tab$DESeq2_log2FC[valid], tab$edgeR_log2FC[valid],
    pch = 16, cex = 0.35,
    col = ifelse(both[valid], "#C83E4D", "#77777755"),
    xlab = "DESeq2 log2 fold change",
    ylab = "edgeR log2 fold change",
    main = paste(name, "- T2D versus ND")
  )
  abline(0, 1, col = "#2878B5", lty = 2)
}
dev.off()

summary <- do.call(rbind, summaries)
write.table(
  summary, file.path(out, "method_comparison_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()), file.path(out, "sessionInfo.txt")
)
print(summary, row.names = FALSE)

```

### scripts/qc_beta_pseudobulk.R
```r
suppressPackageStartupMessages(library(DESeq2))

out <- "results/qc/GSE221156"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(
  "metadata/GSE221156_beta_pseudobulk_samples.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts <- as.matrix(read.delim(
  gzfile("data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz"),
  row.names = 1,
  check.names = FALSE
))

stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(samples$donor_id),
  setequal(colnames(counts), samples$donor_id),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == floor(counts)),
  max(counts) <= .Machine$integer.max
)

samples <- samples[match(colnames(counts), samples$donor_id), ]
rownames(samples) <- samples$donor_id
stopifnot(
  identical(colnames(counts), rownames(samples)),
  all(colSums(counts) == samples$total_counts),
  sum(counts) == 2081765364
)

storage.mode(counts) <- "integer"

samples$disease_status <- factor(
  samples$disease_status, levels = c("ND", "PD", "T2D")
)
samples$chemistry <- factor(samples$chemistry)
samples$sex <- factor(samples$sex)
samples$age <- as.numeric(scale(as.numeric(samples$age)))

# Exploratory QC filter, applied across all 48 donors.
# Differential-testing filters will be documented separately.
keep <- rowSums(counts >= 10) >= 10
if (sum(keep) < 2) stop("Too few expressed genes for QC")

dds <- DESeqDataSetFromMatrix(
  countData = counts[keep, ],
  colData = samples,
  design = ~ chemistry + sex + age + disease_status
)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
transformed <- assay(vsd)

variances <- apply(transformed, 1, var)
top <- head(order(variances, decreasing = TRUE), 500)
pca <- prcomp(t(transformed[top, ]), center = TRUE, scale. = FALSE)
variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)

scores <- data.frame(
  donor_id = rownames(pca$x),
  disease_status = samples$disease_status,
  chemistry = samples$chemistry,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
write.table(
  scores, file.path(out, "PCA_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

distances <- as.matrix(dist(t(transformed)))
median_distances <- apply(
  distances, 1,
  function(x) median(x[x > 0])
)
distance_mad <- mad(median_distances)
review_flag <- if (distance_mad > 0) {
  median_distances >
    median(median_distances) + 3 * distance_mad
} else {
  rep(FALSE, nrow(samples))
}

qc <- data.frame(
  donor_id = samples$donor_id,
  disease_status = samples$disease_status,
  chemistry = samples$chemistry,
  n_beta_cells = samples$n_beta_cells,
  total_counts = colSums(counts),
  detected_genes = colSums(counts > 0),
  counts_per_beta_cell = colSums(counts) / samples$n_beta_cells,
  size_factor = sizeFactors(dds),
  median_VST_distance = median_distances,
  review_flag = review_flag
)
write.table(
  qc, file.path(out, "donor_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

colors <- c(ND = "#2878B5", PD = "#E6A000", T2D = "#C83E4D")
point_colors <- colors[as.character(samples$disease_status)]
point_shapes <- ifelse(samples$chemistry == "V2", 17, 16)

pdf(file.path(out, "beta_pseudobulk_QC.pdf"), width = 11, height = 8.5)

par(mar = c(5, 5, 3, 2))
plot(
  pca$x[, 1], pca$x[, 2],
  col = point_colors, pch = point_shapes,
  xlab = sprintf("PC1 (%.1f%%)", variance_percent[1]),
  ylab = sprintf("PC2 (%.1f%%)", variance_percent[2]),
  main = "Beta-cell donor PCA - top 500 variable genes"
)
text(
  pca$x[, 1], pca$x[, 2],
  labels = samples$donor_id, pos = 3, cex = 0.65
)
legend(
  "topright", legend = names(colors),
  col = colors, pch = 16, bty = "n"
)
legend(
  "bottomright", legend = c("V2", "V3"),
  pch = c(17, 16), bty = "n"
)

par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
boxplot(
  total_counts ~ disease_status, data = qc,
  log = "y", ylab = "Total counts", xlab = "Disease"
)
boxplot(
  n_beta_cells ~ disease_status, data = qc,
  ylab = "Beta cells", xlab = "Disease"
)
plot(
  qc$n_beta_cells, qc$total_counts,
  log = "y", col = point_colors, pch = point_shapes,
  xlab = "Beta cells", ylab = "Total counts"
)
plot(
  qc$n_beta_cells, qc$counts_per_beta_cell,
  col = point_colors, pch = point_shapes,
  xlab = "Beta cells", ylab = "Counts per beta cell"
)

par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
plot(
  hclust(as.dist(distances)),
  main = "Donor clustering using VST expression",
  xlab = "", sub = "", cex = 0.7
)
dev.off()

saveRDS(vsd, file.path(out, "beta_QC_vst.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(out, "sessionInfo.txt")
)

cat("QC donors:", ncol(counts), "\n")
cat("Input genes:", nrow(counts), "\n")
cat("Genes retained for QC:", sum(keep), "\n")
cat("PC1 variance:", round(variance_percent[1], 2), "%\n")
cat("PC2 variance:", round(variance_percent[2], 2), "%\n")
cat("Donors flagged for review:",
    paste(qc$donor_id[qc$review_flag], collapse = ", "), "\n")
cat("No donors were excluded.\n")
cat("Saved:", file.path(out, "beta_pseudobulk_QC.pdf"), "\n")

```

### scripts/summarize_beta_robustness.R
```r
input <- "results/expression/GSE221156"
output <- file.path(input, "robustness")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

analyses <- c("primary", "BMI_sensitivity", "V3_only")

tables <- setNames(lapply(analyses, function(name) {
  path <- file.path(input, paste0(name, "_DESeq2_edgeR.tsv"))
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

  required <- c(
    "gene_id", "gene_symbol",
    "DESeq2_log2FC", "DESeq2_FDR_no_independent_filter",
    "edgeR_log2FC", "edgeR_FDR"
  )
  stopifnot(all(required %in% names(x)))
  stopifnot(!anyNA(x$gene_id), !anyDuplicated(x$gene_id))

  x$supported <- with(
    x,
    is.finite(DESeq2_FDR_no_independent_filter) &
    is.finite(edgeR_FDR) &
    is.finite(DESeq2_log2FC) &
    is.finite(edgeR_log2FC) &
    DESeq2_FDR_no_independent_filter < 0.05 &
    edgeR_FDR < 0.05 &
    sign(DESeq2_log2FC) == sign(edgeR_log2FC)
  )

  rownames(x) <- x$gene_id
  x
}), analyses)

# Keep all genes tested in any analysis; untested genes remain NA.
genes <- sort(unique(unlist(lapply(tables, function(x) x$gene_id))))
result <- data.frame(gene_id = genes, stringsAsFactors = FALSE)

for (name in analyses) {
  x <- tables[[name]]
  idx <- match(genes, x$gene_id)

  if (name == "primary") {
    result$gene_symbol <- x$gene_symbol[idx]
  }

  for (column in c(
    "DESeq2_log2FC", "DESeq2_FDR_no_independent_filter",
    "edgeR_log2FC", "edgeR_FDR", "supported"
  )) {
    result[[paste(name, column, sep = "_")]] <- x[[column]][idx]
  }
}

support_columns <- paste0(analyses, "_supported")
support <- as.matrix(result[, support_columns])
result$n_analyses_supported <- rowSums(support, na.rm = TRUE)

fc_columns <- unlist(lapply(analyses, function(name) {
  paste(name, c("DESeq2_log2FC", "edgeR_log2FC"), sep = "_")
}))
fc <- as.matrix(result[, fc_columns])

result$consistent_direction_all_six <- apply(fc, 1, function(values) {
  all(is.finite(values)) &&
    (all(values > 0) || all(values < 0))
})

result$core_supported <- (
  result$n_analyses_supported == 3 &
  result$consistent_direction_all_six
)

fdr_columns <- unlist(lapply(analyses, function(name) {
  paste(name, c("DESeq2_FDR_no_independent_filter", "edgeR_FDR"), sep = "_")
}))
fdr <- as.matrix(result[, fdr_columns])

result$maximum_FDR_all_six <- apply(fdr, 1, function(values) {
  if (all(is.finite(values))) max(values) else NA_real_
})

core <- result[result$core_supported, ]
core <- core[order(core$maximum_FDR_all_six, core$gene_id), ]

write.table(
  result, file.path(output, "all_genes_robustness.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  core, file.path(output, "core_supported_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary <- data.frame(
  analysis = analyses,
  genes_tested = vapply(tables, nrow, integer(1)),
  supported_both_methods = vapply(
    tables, function(x) sum(x$supported), integer(1)
  )
)
write.table(
  summary, file.path(output, "robustness_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output, "sessionInfo.txt")
)

print(summary, row.names = FALSE)
cat("\nCore genes supported across all six tests:", nrow(core), "\n")
cat("Saved outputs to:", output, "\n")

```

### scripts/umap_beta_donors.R
```r
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(uwot))

out <- "results/qc/GSE221156"
vsd <- readRDS(file.path(out, "beta_QC_vst.rds"))
expression <- assay(vsd)
metadata <- as.data.frame(colData(vsd))

colors <- c(ND = "#2878B5", PD = "#E6A000", T2D = "#C83E4D")
coordinates <- list()

pdf(file.path(out, "donor_UMAP_sensitivity.pdf"), width = 9, height = 7)

for (subset_name in c("all_donors", "V3_only")) {
  selected <- if (subset_name == "all_donors") {
    rep(TRUE, ncol(expression))
  } else {
    metadata$chemistry == "V3"
  }

  x <- expression[, selected, drop = FALSE]
  m <- metadata[selected, , drop = FALSE]
  variable_genes <- head(
    order(apply(x, 1, var), decreasing = TRUE), 500
  )
  pcs <- prcomp(
    t(x[variable_genes, , drop = FALSE]),
    center = TRUE, scale. = FALSE
  )$x[, 1:10, drop = FALSE]

  for (neighbors in c(5, 10, 15)) {
    for (seed in c(11, 42, 2026)) {
      set.seed(seed)
      embedding <- uwot::umap(
        pcs,
        n_neighbors = neighbors,
        min_dist = 0.3,
        metric = "euclidean",
        n_threads = 1,
        n_sgd_threads = 1,
        verbose = FALSE
      )

      key <- paste(subset_name, neighbors, seed, sep = "_")
      coordinates[[key]] <- data.frame(
        subset = subset_name,
        n_neighbors = neighbors,
        seed = seed,
        donor_id = rownames(m),
        disease_status = m$disease_status,
        chemistry = m$chemistry,
        UMAP1 = embedding[, 1],
        UMAP2 = embedding[, 2]
      )

      plot(
        embedding,
        col = colors[as.character(m$disease_status)],
        pch = ifelse(m$chemistry == "V2", 17, 16),
        xlab = "UMAP1", ylab = "UMAP2",
        main = sprintf(
          "%s | neighbors=%d | seed=%d",
          subset_name, neighbors, seed
        )
      )
      legend(
        "topright", legend = names(colors),
        col = colors, pch = 16, bty = "n"
      )
      if (subset_name == "all_donors") {
        legend(
          "bottomright", legend = c("V2", "V3"),
          pch = c(17, 16), bty = "n"
        )
      }
    }
  }
}
dev.off()

write.table(
  do.call(rbind, coordinates),
  file.path(out, "donor_UMAP_coordinates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out, "UMAP_sessionInfo.txt")
)
cat("Saved donor_UMAP_sensitivity.pdf and coordinates.\n")

```

### scripts/validate_beta_h5ad.py
```python
#!/usr/bin/env python3
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE221156/beta_cells.h5ad"
OUT = ROOT / "metadata"
OUT.mkdir(exist_ok=True)


def decode(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def read_column(group, name):
    node = group[name]
    if isinstance(node, h5py.Dataset):
        return [decode(x) for x in node[:]]
    if "categories" in node and "codes" in node:
        categories = [decode(x) for x in node["categories"][:]]
        return [
            categories[int(code)] if code >= 0 else ""
            for code in node["codes"][:]
        ]
    raise ValueError(f"Unsupported column encoding: {name}")


def write_tsv(name, rows, fields):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def normalize_islet(value):
    value = value.strip()
    if value.startswith("Islet"):
        return value
    try:
        number = float(value)
        if number.is_integer():
            return f"Islet{int(number)}"
    except ValueError:
        pass
    return value


def scan_matrix(node, name, chunk=1_000_000):
    encoding = decode(node.attrs.get("encoding-type", ""))
    if encoding != "csr_matrix":
        raise ValueError(f"{name}: expected CSR, found {encoding}")

    data = node["data"]
    shape = [int(x) for x in node.attrs["shape"]]
    indptr = node["indptr"][:]

    if len(indptr) != shape[0] + 1:
        raise ValueError(f"{name}: invalid indptr length")
    if indptr[0] != 0 or indptr[-1] != len(data):
        raise ValueError(f"{name}: invalid indptr endpoints")
    if np.any(indptr[1:] < indptr[:-1]):
        raise ValueError(f"{name}: nonmonotonic indptr")
    if len(node["indices"]) != len(data):
        raise ValueError(f"{name}: indices/data length mismatch")

    nonfinite = negative = fractional = invalid_indices = 0
    minimum, maximum, total = float("inf"), float("-inf"), 0.0

    for start in range(0, len(data), chunk):
        values = data[start:start + chunk]
        indices = node["indices"][start:start + chunk]
        finite = np.isfinite(values)
        valid = values[finite]

        nonfinite += int((~finite).sum())
        negative += int((valid < 0).sum())
        # Exact integrality: no silent rounding.
        fractional += int((valid != np.floor(valid)).sum())
        invalid_indices += int(
            ((indices < 0) | (indices >= shape[1])).sum()
        )

        if valid.size:
            minimum = min(minimum, float(valid.min()))
            maximum = max(maximum, float(valid.max()))
            total += float(valid.sum(dtype=np.float64))

    result = {
        "shape": shape,
        "stored_values": len(data),
        "dtype": str(data.dtype),
        "nonfinite_values": nonfinite,
        "negative_values": negative,
        "fractional_values": fractional,
        "invalid_column_indices": invalid_indices,
        "minimum_stored_value": minimum if len(data) else None,
        "maximum_stored_value": maximum if len(data) else None,
        "sum_stored_values": total,
        "passes_numeric_count_checks": not (
            nonfinite or negative or fractional or invalid_indices
        ),
        "count_provenance_confirmed": False,
    }
    print(f"Finished scanning {name}", flush=True)
    return result


with (ROOT / "GSE221156_audit/metadata/GSE221156_donors.tsv").open() as f:
    expected = {
        row["donor_id"]: row["disease_status"]
        for row in csv.DictReader(f, delimiter="\t")
    }

with (
    ROOT / "GSE221156_audit/metadata/"
    "GSE221156_library_donor_links.tsv"
).open() as f:
    allowed = {
        (row["library_id"], row["donor_id"])
        for row in csv.DictReader(f, delimiter="\t")
    }

disease_map = {
    "normal": "ND",
    "prediabetes syndrome": "PD",
    "type 2 diabetes mellitus": "T2D",
}

with h5py.File(INPUT, "r") as f:
    obs = f["obs"]
    islets = [normalize_islet(x) for x in read_column(obs, "Islet")]
    donors = read_column(obs, "donor_id")
    libraries = read_column(obs, "LibraryID")
    diseases = read_column(obs, "disease")
    cell_types = read_column(obs, "cell_type")

    index_name = decode(obs.attrs.get("_index", "_index"))
    cell_ids = read_column(obs, index_name)

    n = len(islets)
    if any(len(x) != n for x in (
        donors, libraries, diseases, cell_types, cell_ids
    )):
        raise ValueError("Observation column lengths differ")

    counts = Counter()
    donor_sets = defaultdict(set)
    disease_sets = defaultdict(set)
    reverse_donor_sets = defaultdict(set)
    problems = Counter()

    for islet, donor, library, disease in zip(
        islets, donors, libraries, diseases
    ):
        base_library = library.split("_", 1)[0]
        status = disease_map.get(disease, "UNKNOWN")
        counts[(islet, donor, base_library, library, status)] += 1
        donor_sets[islet].add(donor)
        disease_sets[islet].add(status)
        reverse_donor_sets[donor].add(islet)

        if not donor or donor.lower() in {"unknown", "nan", "none"}:
            problems["missing_or_unknown_donor"] += 1
        if islet not in expected:
            problems["unexpected_islet"] += 1
        elif status != expected[islet]:
            problems["disease_mismatch"] += 1
        if (base_library, islet) not in allowed:
            problems["invalid_library_donor_membership"] += 1
        if "_" in library and library.split("_", 1)[1] != islet:
            problems["library_suffix_islet_mismatch"] += 1

    problems["islets_with_multiple_donor_ids"] = sum(
        len(x) != 1 for x in donor_sets.values()
    )
    problems["donor_ids_with_multiple_islets"] = sum(
        len(x) != 1 for x in reverse_donor_sets.values()
    )
    problems["islets_with_multiple_disease_states"] = sum(
        len(x) != 1 for x in disease_sets.values()
    )
    problems["duplicate_cell_ids"] = n - len(set(cell_ids))

    rows = [
        dict(zip(
            ["islet", "donor_id", "base_library_id",
             "annotated_library_id", "disease_status", "n_cells"],
            (*key, count),
        ))
        for key, count in sorted(counts.items())
    ]
    write_tsv(
        "GSE221156_beta_library_donor_cells.tsv",
        rows,
        ["islet", "donor_id", "base_library_id",
         "annotated_library_id", "disease_status", "n_cells"],
    )

    summary = {
        "n_cells": n,
        "observed_islets": len(donor_sets),
        "observed_donor_ids": len(reverse_donor_sets),
        "missing_expected_islets": sorted(set(expected) - set(islets)),
        "observed_cell_types": dict(Counter(cell_types)),
        "metadata_problems": dict(problems),
        "matrices": {},
    }

    for key in ("X", "raw/X"):
        if key in f:
            summary["matrices"][key] = scan_matrix(f[key], key)

            var = f["raw/var"] if key == "raw/X" else f["var"]
            var_index = decode(var.attrs.get("_index", "_index"))
            genes = read_column(var, var_index)
            summary["matrices"][key]["n_feature_ids"] = len(genes)
            summary["matrices"][key]["duplicate_feature_ids"] = (
                len(genes) - len(set(genes))
            )
            summary["matrices"][key]["shape_matches_metadata"] = (
                summary["matrices"][key]["shape"] == [n, len(genes)]
            )

summary["metadata_pass"] = (
    not any(problems.values())
    and not summary["missing_expected_islets"]
    and summary["observed_islets"] == 48
    and summary["observed_donor_ids"] == 48
)

output = OUT / "GSE221156_beta_validation.json"
output.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
print(f"Saved: {output}")

```

