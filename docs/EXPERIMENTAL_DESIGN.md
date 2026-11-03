# Experimental Design by Aim (v2)

**Project:** t2d-islet-regulatory-genomics
**Discovery data:** GSE221156 beta-cell pseudobulk (48 donors; primary contrast 17 ND vs 17 T2D)
**Status basis:** Progress report of 2026-10-04 (snapshot `fa33e5f`) and Phase 8 deliverables. Nothing below reports a new result; items marked **DONE** summarize work already in the repository, everything else is planned.

---

## 0. Central question and novelty framing

**Working title:** *Separating cause from consequence in the type 2 diabetes beta-cell transcriptome: a replicated proteasome/ubiquitin module and its genetic anchoring*

**Central question:** Which T2D-associated beta-cell expression changes are (a) robust to method and model, (b) reproducible in independent donors, and (c) supported by human genetics as plausibly *upstream* of disease rather than a *consequence* of hyperglycaemia, treatment, or donor/islet-processing factors?

**Why this framing:** The source study (Bandesh et al., EMBO J 2026) already reports beta-cell expression defects in these donors and integrates genetics. Re-running DE plus an eQTL overlap would not be novel. The defensible contributions are:

1. A pre-specified, independently replicated pathway-level module (proteasome/ubiquitin), tested with donor-aware statistics.
2. A disease-stage gradient using the 14 prediabetes (PD) donors, which the primary contrast does not use.
3. An explicit cause-versus-consequence classification of DE genes using genetic evidence with direction of effect.

Read the source study's complete methods and supplementary tables before finalizing this claim. If it already reports a proteasome finding, novelty must rest on points 2–3 and on replication.

---

## 1. Global rules (apply to all aims)

1. **Donors are the unit of replication.** Never test cells as replicates.
2. **Freeze before you look.** The 29-gene module file, its SHA256, the validation plan, and the decision rules are committed and tagged (`prereg-validation-v1`) **before** any replication count matrix is downloaded. Any later change goes in `docs/deviations.md`.
3. **Direction convention:** ND is the reference. Negative log2FC means lower in T2D. Use this convention in every table and figure.
4. **Thresholds live in `config/analysis.yaml`.** Scripts never hard-code them.
5. **Every result traces to a script plus a log.** Use the existing `Rscript … | tee logs/…` pattern, and add `sessionInfo()` / `pip freeze` to each log.

---

## Aim 1 — Discovery robustness in GSE221156

**Hypothesis 1:** T2D-associated beta-cell expression changes are stable to statistical method, BMI adjustment, sequencing chemistry, and individual donors.

### Completed (DONE)
- Metadata audit, library-to-donor resolution, and count-conserving pseudobulk (2,081,765,364 counts preserved).
- Design rank checks for three models: primary, BMI-adjusted, and V3-only.
- PCA/UMAP QC, which shows chemistry-associated structure.
- DESeq2 vs edgeR robust QL in each model. This gives 913 jointly significant genes in the primary model and 285 that hold up across all six method-model combinations.

### Remaining steps
| Step | Analysis | Why reviewers will ask | Output |
|---|---|---|---|
| 1.1 | **Leave-one-donor-out (LODO)** refit for the 285 core genes and the module score; report sign flips and max |ΔLFC| | 17 vs 17 is small; one donor can drive a result | `results/expression/GSE221156/robustness/lodo.tsv` |
| 1.2 | **Disease × chemistry interaction** test (primary model plus interaction term) | V2 is unbalanced (4 ND vs 7 T2D); additive chemistry adjustment assumes the same disease effect in V2 and V3 | `…/chemistry_interaction.tsv` |
| 1.3 | **Islet-processing covariates** where available: culture time, islet purity, cold ischaemia time, cause of death. Run as sensitivity models, one covariate at a time | Common non-disease drivers of islet transcriptomes | `…/processing_sensitivity.tsv` |
| 1.4 | **Recovered beta-cell count** as a covariate sensitivity (log cells per donor) | Pseudobulk depth and composition vary 582–5,294 cells | `…/ncell_sensitivity.tsv` |
| 1.5 | **Low-abundance review**: flag genes with control median CPM < 1 (e.g., TBX2, FAIM2, BARX1) and exclude them from headline claims unless they pass LODO | Fold changes from near-zero baselines are fragile | `…/low_abundance_flags.tsv` |
| 1.6 | **Donor-level candidate plots** (disease × chemistry, each donor a point) | Already in your next priorities | `figures/expression/candidates/*.pdf` |

**Success criterion:** A "core" gene list where each gene keeps the same sign in every LODO fit and shows no supported disease × chemistry interaction (interaction FDR ≥ 0.05). Report the list size, whatever it is.

---

## Aim 2 — Disease-stage gradient (ND → PD → T2D)

**Hypothesis 2:** If a change is an early contributor rather than a late consequence of overt hyperglycaemia, it is partly present in prediabetes. That means a monotonic ND → PD → T2D trend, or an association with HbA1c within non-T2D donors.

This aim uses the 14 PD donors that are currently unused. It is the cheapest new evidence you can add.

| Step | Analysis | Output |
|---|---|---|
| 2.1 | Ordinal trend model: `~ chemistry + sex + age_z + stage_numeric` (ND=0, PD=1, T2D=2), all 48 donors | `results/expression/GSE221156/stage_trend.tsv` |
| 2.2 | PD vs ND contrast (exploratory; expect low power) | `…/PD_vs_ND.tsv` |
| 2.3 | HbA1c as a continuous predictor within ND + PD only (excludes T2D treatment effects) | `…/HbA1c_nonT2D.tsv` |
| 2.4 | Module-level versions of 2.1–2.3 (see Aim 3 for the module test) | `…/module_stage.tsv` |

**Interpretation rule (pre-specify):** A gene or module is labelled *early* if the trend coefficient and the PD-vs-ND point estimate are in the same direction as T2D-vs-ND. It is labelled *late/consequence-compatible* if PD ≈ ND (the confidence interval includes 0 and the point estimate is less than 25% of the T2D effect). A null result in PD is reported as "not detected," not as "absent."

---

## Aim 3 — Pathway structure and the proteasome/ubiquitin module (discovery)

**Hypothesis 3:** Proteasome and ubiquitin genes show a coordinated decrease in T2D beta cells beyond what is expected from inter-gene correlation.

**Key risk: selection circularity.** If the 29 genes were chosen as the leading edge of an enrichment test on GSE221156, their discovery p-value is optimistically biased. In that case, discovery evidence for the 29-gene set is descriptive only, and inference rests on Aim 4.

| Step | Analysis | Output |
|---|---|---|
| 3.1 | Document exactly how the 29 genes were derived: gene-set source and version (MSigDB/Reactome release), ranking statistic, leading-edge rule. Freeze as `config/gene_sets/proteasome29.tsv` with SHA256 | `docs/module_derivation.md` |
| 3.2 | Use the **full curated proteasome gene set** (not the 29) as the primary discovery pathway test, with a correlation-aware competitive test (`limma::camera` on voom/edgeR-normalized pseudobulk, plus `fry` as a self-contained test) | `results/pathway/GSE221156/camera_primary.tsv` |
| 3.3 | Ranked GSEA (fgsea) on the primary-model statistic as a secondary analysis; record gene-set version and seed | `…/fgsea_primary.tsv` |
| 3.4 | Redundancy reduction of significant pathways (leading-edge overlap clustering) | `…/pathways_nonredundant.tsv` |
| 3.5 | Repeat 3.2 for BMI and V3-only models | `…/camera_sensitivity.tsv` |
| 3.6 | Per-donor module score (mean z-scored log-CPM of the 29 genes) for plotting and for Aims 2 and 4 | `…/module_scores_by_donor.tsv` |

**Specificity control:** Compare against size-matched random gene sets with expression levels matched to the module (1,000 draws). Also test a beta-cell stress/UPR set to see whether the proteasome signal is separable from general ER stress.

---

## Aim 4 — Independent replication (prospective)

This builds on your `Independent_Validation_Plan.md`, which is already well specified. The additions below cover power, overlap, pooling, and in vivo relevance.

### 4.1 Cohort ladder
| Tier | Dataset | Notes |
|---|---|---|
| A | GSE153855 (6 ND / 5 T2D) | First choice; verify raw counts and beta-cell labels |
| A | E-MTAB-5061 (single-cell, ND/T2D) | Count-based; check donor independence |
| B | GSE81608 (RPKM) | Use limma on log-transformed donor means, or reprocess FASTQ; never DESeq2 on RPKM |
| C | Human Pancreas Analysis Program (HPAP / PANC-DB) single-cell data | Larger donor numbers; **mandatory donor-overlap check against GSE221156** |
| C | Bulk islet RNA-seq cohorts with HbA1c (e.g., Lund/Fadista-type cohorts) | Not beta-specific; adjust for estimated beta-cell proportion by deconvolution |
| D | Islets from living surgical donors (partial pancreatectomy cohorts) | Tests whether the signal persists without brain-death and isolation artefacts |

### 4.2 Added safeguards
1. **Overlap audit first.** Several public islet cohorts draw on the same organ-procurement networks. Match on donor ID, age, sex, BMI, and HbA1c, and exclude any probable duplicates. Record the result in `metadata/validation/overlap_audit.tsv`.
2. **Power statement.** With about 11 donors, single-gene replication is underpowered. The **module-level** test is primary. Report the minimum detectable standardized module effect for each cohort.
3. **Pooling.** Pool donor-level module effects across Tier A/B cohorts with a random-effects meta-analysis (`metafor`), reporting I² and leave-one-cohort-out estimates. Pre-specify that the pooled estimate, not any single small cohort, is the primary replication endpoint.
4. **Per-gene concordance.** Report sign agreement with an exact binomial test as descriptive only (the genes are correlated), plus Spearman correlation of effects with a bootstrap CI over donors.

### 4.3 Decision matrix (pre-registered)
| Pooled module effect | Gene-set test | Conclusion wording |
|---|---|---|
| Negative, CI excludes 0 | Supported in ≥ 1 independent cohort | "Replicated in independent donors" |
| Negative, CI includes 0 | Not supported | "Directionally consistent; not statistically replicated" |
| Positive, CI excludes 0 | Supported, opposite sign | "Failed to replicate; opposite direction" (report prominently) |

---

## Aim 5 — Genetic anchoring: cause versus consequence

**Hypothesis 5:** A subset of robust DE genes, possibly including module members, lies under T2D genetic signals. For that subset, the T2D risk allele shifts islet expression in the **same direction** as the T2D-vs-ND expression change, which supports an upstream role. Genes without such evidence are classified as consequence-compatible or unresolved.

| Step | Analysis | Output |
|---|---|---|
| 5.1 | Audit and harmonize T2D GWAS (multi-ancestry T2DGGI, DIAMANTE ancestry-specific) and islet eQTL (InsPIRE, TIGER) to one genome build; record versions in `metadata/data_manifest.tsv`; fill `harmonized_genome_build` in config | `data/processed/gwas/…`, `docs/genetics_sources.md` |
| 5.2 | **Gene-set heritability test:** MAGMA gene-set analysis and S-LDSC on ±100 kb windows around (a) core DE genes, (b) up vs down genes separately, (c) the proteasome set, against expression-matched control genes | `results/genetics/geneset_h2.tsv` |
| 5.3 | **Signal-aware colocalization** of T2D GWAS with islet eQTL for every core gene with a nearby (±1 Mb) T2D signal (coloc-SuSiE; PP.H4 > 0.8 from config; ancestry-matched LD) | `results/genetics/coloc_islet.tsv` |
| 5.4 | **Direction-of-effect test** at colocalized genes: does the risk allele raise or lower islet expression, and does that match the DE sign? | `results/genetics/direction_concordance.tsv` |
| 5.5 | **TWAS** (islet expression weights, e.g., FUSION/S-PrediXcan) as a complementary gene-level test | `results/genetics/twas_islet.tsv` |
| 5.6 | **Mendelian randomization** using cis-eQTL instruments with Steiger filtering; SMR/HEIDI to flag linkage | `results/genetics/mr_smr.tsv` |
| 5.7 | Prior-sensitivity analysis for coloc (config already has `require_prior_sensitivity: true`) | `results/genetics/coloc_prior_sensitivity.tsv` |

### Gene classification (final deliverable of this aim)
| Class | Evidence required |
|---|---|
| **Genetically anchored, concordant** | Coloc PP.H4 > 0.8 **and** risk allele direction matches DE direction **and** passes Aim 1 robustness |
| **Genetically anchored, discordant** | Coloc supported but direction opposite (consistent with compensation; report separately) |
| **Consequence-compatible** | Robust DE, no genetic support, and late pattern in Aim 2 |
| **Unresolved** | Everything else, including underpowered loci |

**Expectation management:** Most DE genes will not be genetically anchored, and that is itself a finding. Heritability enrichment for downregulated genes but not upregulated genes (or the reverse) would be a meaningful result about cause versus consequence.

---

## Aim 6 — Regulatory mechanism at anchored genes

Run only for genes in the "genetically anchored, concordant" class.

| Step | Analysis | Output |
|---|---|---|
| 6.1 | Fine-map the colocalized signal (SuSiE-RSS); intersect credible variants with beta-cell snATAC peaks (Chiou et al. 2021) | `results/regulatory/cs_beta_peaks.tsv` |
| 6.2 | Link peaks to the gene (co-accessibility, islet promoter-capture Hi-C / enhancer hubs) | `results/regulatory/peak_gene_links.tsv` |
| 6.3 | Transcription-factor motif disruption for credible variants (motifbreakR, JASPAR 2024) | `results/regulatory/motif_disruption.tsv` |
| 6.4 | Allelic imbalance in islet chromatin or expression data, where available | `results/regulatory/allelic_imbalance.tsv` |

---

## Aim 7 — Functional validation (proposed; for collaborators or the Discussion)

Only for one or two top anchored genes and/or the module as a whole.

- **Model:** EndoC-βH1 or EndoC-βH5 human beta-cell lines; primary human islets if available.
- **Perturbation:** siRNA/CRISPRi knockdown of anchored genes. For the module, use partial proteasome inhibition with a titrated dose that does not cause overt cytotoxicity, alongside a viability control.
- **Readouts:** glucose-stimulated insulin secretion (2.8 vs 16.7 mM glucose), insulin content, proteasome activity (chymotrypsin-like assay), UPR markers (spliced XBP1, ATF4, CHOP), and apoptosis.
- **Design:** at least three independent biological replicates (separate passages or donors), non-targeting controls, a pre-specified primary readout (GSIS stimulation index), and a power calculation from pilot variance.
- **Enhancer test (if Aim 6 yields a variant):** allele-specific luciferase assay, then CRISPRi of the enhancer with target-gene readout.

---

## 2. Figure plan

| Fig. | Content | Aim | Source |
|---|---|---|---|
| 1 | Study design and donor/cohort flow, with the evidence chain | all | schematic |
| 2 | Discovery DE: method concordance, robustness funnel (913 → 285 → LODO-stable) | 1 | existing robustness outputs + 1.1 |
| 3 | Pathway landscape and proteasome module with donor-level scores by chemistry | 3 | 3.2–3.6 |
| 4 | ND → PD → T2D gradient for module and core genes | 2 | 2.1–2.4 |
| 5 | Replication forest plot (cohort and pooled module effects) | 4 | 4.2 |
| 6 | Genetic anchoring: heritability enrichment and cause/consequence classification | 5 | 5.2–5.4 |
| 7 | Example anchored locus: GWAS, credible set, beta-cell chromatin, eQTL | 6 | 6.1–6.2 |

---

## 3. Repository actions arising from this design

1. **Update README "Status."** It still says no findings have been generated, but Phases 1–8 are complete. Add a phase table with the status of each phase.
2. **Update `config/analysis.yaml`:** change `status: planning` to `analysis`, set `harmonized_genome_build`, and add sections for `pathway`, `validation`, `stage_gradient`, and `genetics.coloc_priors`.
3. **Freeze and tag the validation plan now:** commit `Independent_Validation_Plan.md`, the decision matrix (Aim 4.3), and `proteasome29.tsv` with its hash, then `git tag -a prereg-validation-v1`. This must happen before any replication data are downloaded.
4. **Lock environments:** `renv::snapshot()` → `renv.lock`; `pip freeze > requirements-audit.txt`.
5. **Add a single entry point:** a `Makefile` or Snakemake workflow that chains the existing scripts in Appendix A order, so a reviewer can run one command.
6. **Add `LICENSE`, `CITATION.cff`, and `CHANGELOG.md`.** Connect the repo to Zenodo, and release `v1.0.0` at submission.
7. **Add `docs/deviations.md`** and `docs/module_derivation.md` (Aim 3.1).
