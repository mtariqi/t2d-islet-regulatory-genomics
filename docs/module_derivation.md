# Derivation of the 29-gene proteasome/ubiquitin module

**Frozen file:** config/gene_sets/proteasome29.tsv (SHA256 in proteasome29.sha256)
**Source file:** results/enrichment/GSE221156/Reactome_negative18_shared29_gene_characterization.tsv
**Generating script:** scripts/enrichment/characterize_shared_reactome_genes.R
**Sensitivity audit:** scripts/enrichment/audit_proteasome_sensitivity.R
**Discovery data:** GSE221156 beta-cell donor pseudobulk, primary model (17 ND vs 17 T2D)

## Selection rule
1. Ranked Reactome enrichment (fgsea) on the primary-model DE statistic (DESeq2 and edgeR).
2. Negatively enriched pathways (lower in T2D) were grouped into redundancy clusters
   (results/enrichment/GSE221156/pathway_redundancy_clusters_v2.tsv).
3. Cluster "Reactome.negative_18" contains 38 Reactome pathways.
4. The module is the 29 genes present in the fgsea leading edge of all 38 pathways.

Reactome release: not recorded

## Discovery behaviour (descriptive; not independent evidence)
| Model | Negative in both methods | Significant in both | Median DESeq2 log2FC |
|---|---:|---:|---:|
| Primary | 29/29 | 2/29 | -0.182 |
| BMI-adjusted | 29/29 | 0/29 | -0.181 |
| V3 only | 28/29 | 0/29 | -0.082 |

PSMA3 was not negative in all models.

## Interpretation constraints
- Genes were selected from discovery leading edges, so discovery p-values for this set are optimistically biased and are not evidence for the module.
- The effect is a small coordinated shift; module-level tests are primary, per-gene results are descriptive.
- The effect is attenuated in V3-only libraries; a disease x chemistry interaction test is required before claiming robustness.
- Inference rests on independent replication (docs/Independent_Validation_Plan.md).

## Pre-specified replication comparator
In each replication cohort, the full Reactome proteasome pathway (all annotated members detected in that cohort) is tested alongside the 29-gene module. Its definition does not depend on discovery ranks.
