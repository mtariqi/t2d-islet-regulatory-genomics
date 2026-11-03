#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
base <- 'results/enrichment/GSE221156'
expr <- 'results/expression/GSE221156'
out <- 'figures/phase8/publication'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
readtsv <- function(f) {if(!file.exists(f)) stop('Missing: ',f);read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)}
models <- c(primary='primary',BMI='BMI_sensitivity',V3='V3_only')
cluster <- readtsv(file.path(base,'proteasome38_cross_model_robustness.tsv'))
stopifnot(nrow(cluster)==38)
# FIGURE 1: pathway-level paired NES values across models, with descriptive medians.
nes <- do.call(rbind,lapply(names(models),function(label){
 m<-models[[label]]
 do.call(rbind,lapply(c('DESeq2','edgeR'),function(method){
  data.frame(model=label,method=method,pathway=cluster$pathway,
             NES=cluster[[paste0(m,'_NES_',method)]])
 }))
}))
nes$model<-factor(nes$model,levels=names(models))
nes$method<-factor(nes$method,levels=c('DESeq2','edgeR'))
p1<-ggplot(nes,aes(model,NES,group=pathway))+geom_hline(yintercept=0,linetype=2,color='grey60')+
 geom_line(alpha=.10,color='grey50')+geom_point(alpha=.28,size=.85,color='#486a87')+
 stat_summary(aes(group=1),fun=median,geom='point',size=3.2,color='#bd4b3a')+
 stat_summary(aes(group=1),fun=median,geom='line',linewidth=.9,color='#bd4b3a')+
 facet_wrap(~method,nrow=1)+theme_bw(base_size=12)+
 labs(x=NULL,y='Normalized enrichment score (NES)',title='Proteasome-associated Reactome cluster across models',
 subtitle='38 discovery-selected overlapping pathways; red = median NES')
ggsave(file.path(out,'Figure1_Proteasome38_NES.pdf'),p1,width=10,height=5.5)
ggsave(file.path(out,'Figure1_Proteasome38_NES.png'),p1,width=10,height=5.5,dpi=300)
# FIGURE 2: exact 29 shared genes from characterized primary file.
characterization <- readtsv(file.path(base,'Reactome_negative18_shared29_gene_characterization.tsv'))
stopifnot('gene_id'%in%names(characterization),nrow(characterization)==29)
ids<-unique(characterization$gene_id)
heat <- do.call(rbind,lapply(names(models),function(label){
 m<-models[[label]]
 x<-readtsv(file.path(expr,paste0(m,'_DESeq2_edgeR.tsv')))
 stopifnot(all(c('gene_id','gene_symbol','DESeq2_log2FC','edgeR_log2FC')%in%names(x)))
 x<-x[match(ids,x$gene_id),]
 if(anyNA(x$gene_id)) stop('Missing 29 gene IDs in ',m)
 rbind(data.frame(gene_id=x$gene_id,gene=x$gene_symbol,model=label,method='DESeq2',log2FC=x$DESeq2_log2FC),
       data.frame(gene_id=x$gene_id,gene=x$gene_symbol,model=label,method='edgeR',log2FC=x$edgeR_log2FC))
}))
heat$gene[is.na(heat$gene)|!nzchar(heat$gene)]<-heat$gene_id[is.na(heat$gene)|!nzchar(heat$gene)]
heat$gene<-factor(heat$gene,levels=rev(unique(heat$gene)))
heat$specification<-factor(paste(heat$model,heat$method,sep=' / '),
 levels=as.vector(t(outer(names(models),c('DESeq2','edgeR'),paste,sep=' / '))))
lim<-max(abs(heat$log2FC),na.rm=TRUE)
p2<-ggplot(heat,aes(specification,gene,fill=log2FC))+geom_tile(color='white',linewidth=.25)+
 scale_fill_gradient2(low='#3468a0',mid='white',high='#c54c40',midpoint=0,limits=c(-lim,lim))+
 theme_minimal(base_size=11)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+
 labs(x=NULL,y=NULL,fill='log2FC',title='Shared proteasome/ubiquitin genes across models',
 subtitle='Negative values indicate lower expression in T2D; color is not significance')
ggsave(file.path(out,'Figure2_Proteasome29_log2FC_heatmap.pdf'),p2,width=10,height=10)
ggsave(file.path(out,'Figure2_Proteasome29_log2FC_heatmap.png'),p2,width=10,height=10,dpi=300)
# FIGURE 3: nonredundant cluster representatives, sourced from audited v2 representative table.
repfile<-file.path(base,'representative_pathways_v2.tsv')
r<-readtsv(repfile)
cat('Representative table columns:',paste(names(r),collapse=', '),'\n')
if(!('pathway'%in%names(r)) || !('collection'%in%names(r))) stop('Inspect representative file columns before figure 3')
primary<-readtsv(file.path(base,'pathway_method_consensus.tsv'))
if(!all(c('pathway','collection','NES_DESeq2','NES_edgeR','consensus_supported','padj_DESeq2','padj_edgeR')%in%names(primary))) stop('Unexpected primary consensus schema')
x<-merge(unique(r[,c('pathway','collection')]),primary,by=c('pathway','collection'))
x<-x[!is.na(x$consensus_supported)&x$consensus_supported,]
x$worst_fdr<-pmax(x$padj_DESeq2,x$padj_edgeR)
x<-x[order(x$worst_fdr),]
x<-head(x,20)
if(nrow(x)==0) stop('No representative consensus pathways')
x$mean_NES<-(x$NES_DESeq2+x$NES_edgeR)/2
x$short<-gsub('^(REACTOME_|GOBP_)','',x$pathway)
x$short<-gsub('_',' ',x$short)
x$short<-substr(x$short,1,68)
x$short<-make.unique(x$short)
x$short<-factor(x$short,levels=rev(x$short))
p3<-ggplot(x,aes(mean_NES,short,fill=collection))+geom_col(width=.7)+geom_vline(xintercept=0,linetype=2)+
 theme_bw(base_size=11)+theme(legend.position='bottom')+
 labs(x='Mean NES (DESeq2 and edgeR; descriptive)',y=NULL,
 title='Nonredundant primary pathway representatives',
 subtitle='Top 20 complete-linkage cluster representatives by worse method FDR')
ggsave(file.path(out,'Figure3_Nonredundant_pathways.pdf'),p3,width=11,height=8.5)
ggsave(file.path(out,'Figure3_Nonredundant_pathways.png'),p3,width=11,height=8.5,dpi=300)
write.table(x,file.path(out,'Figure3_representative_source.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Three figure pairs created under ',out,'\n',sep='')
