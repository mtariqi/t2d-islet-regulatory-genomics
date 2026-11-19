#!/usr/bin/env Rscript
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       scripts/influence/leave_one_donor_out.R
# Description:  Leave-one-donor-out refits to test the stability of core-gene effects.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================

# Phase 7: donor-deletion sensitivity using fixed primary gene universe.
suppressPackageStartupMessages({library(DESeq2); library(edgeR)})
output <- 'results/influence/GSE221156'
dir.create(output, recursive=TRUE, showWarnings=FALSE)
samples <- read.delim('metadata/GSE221156_beta_pseudobulk_samples.tsv', check.names=FALSE)
counts <- as.matrix(read.delim(gzfile('data/processed/GSE221156/beta_pseudobulk_counts.tsv.gz'), row.names=1, check.names=FALSE))
stopifnot(!anyDuplicated(samples$donor_id), !anyDuplicated(rownames(counts)),
          setequal(colnames(counts), samples$donor_id), all(is.finite(counts)),
          all(counts >= 0), all(counts == floor(counts)), max(counts) <= .Machine$integer.max)
samples <- samples[match(colnames(counts), samples$donor_id), ,drop=FALSE]
stopifnot(all(colSums(counts) == samples$total_counts))
storage.mode(counts) <- 'integer'
m <- samples[samples$disease_status %in% c('ND','T2D'), ,drop=FALSE]
rownames(m) <- m$donor_id
m$disease_status <- factor(m$disease_status, levels=c('ND','T2D'))
m$chemistry <- factor(m$chemistry)
m$sex <- factor(m$sex)
m$age_z <- as.numeric(scale(as.numeric(m$age)))
formula <- ~ chemistry + sex + age_z + disease_status
full_design <- model.matrix(formula, data=m)
stopifnot(qr(full_design)$rank == ncol(full_design), !anyNA(full_design), nrow(m)==34)
cts <- counts[, rownames(m),drop=FALSE]
# Keep the original primary tested-gene universe fixed in every deletion.
primary <- read.delim('results/expression/GSE221156/primary_DESeq2_edgeR.tsv', check.names=FALSE)
stopifnot(!anyDuplicated(primary$gene_id), all(primary$gene_id %in% rownames(cts)))
cts <- cts[primary$gene_id,,drop=FALSE]
fit_one <- function(drop_id=NA_character_) {
  sub <- if(is.na(drop_id)) m else m[rownames(m)!=drop_id,,drop=FALSE]
  sub$chemistry <- droplevels(sub$chemistry)
  sub$sex <- droplevels(sub$sex)
  sub$disease_status <- factor(as.character(sub$disease_status), levels=c('ND','T2D'))
  # Retain full-cohort age scaling for comparable conditional effects.
  design <- model.matrix(formula, data=sub)
  if(qr(design)$rank != ncol(design)) stop('Rank-deficient design after removing: ',drop_id)
  x <- cts[,rownames(sub),drop=FALSE]
  dds <- DESeqDataSetFromMatrix(countData=x, colData=sub, design=formula)
  dds <- DESeq(dds, minReplicatesForReplace=Inf, quiet=TRUE)
  ds <- results(dds, contrast=c('disease_status','T2D','ND'), alpha=.05, independentFiltering=FALSE)
  y <- DGEList(counts=x)
  y <- calcNormFactors(y, method='TMM')
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  er <- topTags(glmQLFTest(fit, coef=match('disease_statusT2D',colnames(design))), n=Inf,sort.by='none')$table
  er <- er[rownames(x),,drop=FALSE]
  data.frame(gene_id=rownames(x), omitted_donor=if(is.na(drop_id)) 'FULL' else drop_id,
             DESeq2_log2FC=ds$log2FoldChange, DESeq2_FDR=ds$padj,
             edgeR_log2FC=er$logFC, edgeR_FDR=er$FDR, row.names=NULL)
}
full <- fit_one()
# Guard against silent drift from original analysis.
orig <- primary[match(full$gene_id, primary$gene_id),]
stopifnot(isTRUE(all.equal(full$DESeq2_log2FC,orig$DESeq2_log2FC,tolerance=1e-5,check.attributes=FALSE)),
          isTRUE(all.equal(full$edgeR_log2FC,orig$edgeR_log2FC,tolerance=1e-5,check.attributes=FALSE)))
write.table(full,file.path(output,'full_refit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
results_list <- lapply(rownames(m),function(id){message('Omitting ',id); fit_one(id)})
all_deletions <- do.call(rbind,results_list)
write.table(all_deletions,file.path(output,'all_donor_deletions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
summary <- do.call(rbind,lapply(split(all_deletions,all_deletions$gene_id),function(z){
  base <- full[match(z$gene_id[1],full$gene_id),]
  data.frame(gene_id=z$gene_id[1], full_DESeq2_log2FC=base$DESeq2_log2FC,
    full_edgeR_log2FC=base$edgeR_log2FC,
    max_abs_delta_DESeq2=max(abs(z$DESeq2_log2FC-base$DESeq2_log2FC),na.rm=TRUE),
    max_abs_delta_edgeR=max(abs(z$edgeR_log2FC-base$edgeR_log2FC),na.rm=TRUE),
    direction_consistency_DESeq2=mean(sign(z$DESeq2_log2FC)==sign(base$DESeq2_log2FC),na.rm=TRUE),
    direction_consistency_edgeR=mean(sign(z$edgeR_log2FC)==sign(base$edgeR_log2FC),na.rm=TRUE),
    fraction_FDR05_DESeq2=mean(z$DESeq2_FDR<.05,na.rm=TRUE),
    fraction_FDR05_edgeR=mean(z$edgeR_FDR<.05,na.rm=TRUE),
   valid_logFC_DESeq2=sum(is.finite(z$DESeq2_log2FC)),
   valid_logFC_edgeR=sum(is.finite(z$edgeR_log2FC)),
   valid_FDR_DESeq2=sum(!is.na(z$DESeq2_FDR)),
   valid_FDR_edgeR=sum(!is.na(z$edgeR_FDR)))
}))
write.table(summary,file.path(output,'gene_influence_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
message('Saved ',nrow(all_deletions),' deletion-gene results and ',nrow(summary),' gene summaries to ',output)
