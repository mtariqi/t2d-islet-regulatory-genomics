# Phase 7 — Leave-One-Donor-Out Robustness Analysis

## Research objective

Assess the sensitivity of T2D-associated beta-cell
differential expression to individual donor removal
in GSE221156.

## Study design

- Comparison: non-diabetic versus T2D
- Donors: 34
- Tested genes: 17,352
- Statistical methods: DESeq2 and edgeR
- Model: ~ chemistry + sex + age_z + disease_status
- Sensitivity approach: leave-one-donor-out refitting
- Donor-deletion refits: 34 per statistical method

## Reproducibility

The full-cohort log2 fold changes reproduced the
original primary analysis within tolerance 1e-5.

The complete LODO analysis exited successfully.

## Results

| Metric | DESeq2 | edgeR |
|---|---:|---:|
| Genes with complete log2FC estimates | 17,352 | 17,352 |
| Genes with 100% direction consistency | 12,849 | 12,780 |
| Genes with >=90% FDR retention | 792 | 662 |
| Median maximum absolute log2FC change | 0.06777 | 0.06793 |
| Maximum absolute log2FC change | 8.43381 | 3.56607 |

Total donor-gene deletion results: 589,968.

## Interpretation

Most tested genes exhibited relatively small
maximum absolute changes in estimated log2FC.

A subset showed substantial donor sensitivity,
requiring further investigation of expression
abundance, individual donor contributions, and
statistical significance stability.

The direction-consistency and FDR-retention
statistics are descriptive sensitivity metrics.

## Limitations

- LODO is an internal robustness assessment,
  not independent biological replication.
- Removing donors reduces sample size and power.
- Statistical significance retention is not
  equivalent to effect-size stability.
- Genes with low expression or large estimated
  effects require additional scrutiny.
- The 285 sensitivity-supported genes have not
  yet been separately evaluated using LODO.
- These findings do not establish causal
  genetic or regulatory mechanisms.

## Next steps

1. Evaluate the 285 sensitivity-supported genes.
2. Identify influential donor-gene combinations.
3. Investigate low-expression outliers.
4. Generate publication-quality robustness figures.
5. Perform ranked pathway enrichment.
6. Integrate T2D GWAS and islet molecular QTL evidence.
