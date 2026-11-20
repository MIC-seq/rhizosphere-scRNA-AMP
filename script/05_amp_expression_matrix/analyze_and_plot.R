#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat); library(ggplot2)
})
set.seed(20260719)
root <- "/public3/labmember/qianqh/soil_amp_analysis/10_scAMP_8samples"
cfg <- fread(file.path(root,"00_config/samples.tsv"))
ref <- fread(file.path(root,"01_reference/amp_reference_metadata.tsv"))
stats_dir <- file.path(root,"04_stats"); fig_dir <- file.path(root,"05_figures")
dir.create(stats_dir,recursive=TRUE,showWarnings=FALSE); dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
sample_levels <- cfg$sample; crop_levels <- c("soybean","rice","wheat")
cols_source <- c(APD6="#2166AC",AMPSphere="#FDAE61",Novel="#1A9850")
cols_crop <- c(soybean="#4DAF4A",rice="#377EB8",wheat="#E6AB02")

read_first_col <- function(path) fread(path,header=FALSE,select=1)[[1]]
read_annotation <- function(path, name) {
  z <- fread(path,header=FALSE,select=1:2); setnames(z,c("cell_barcode",name)); unique(z,by="cell_barcode")
}
cells <- rbindlist(lapply(seq_len(nrow(cfg)),function(i){
  s <- cfg$sample[i]; w <- data.table(cell_barcode=read_first_col(cfg$whitelist[i]))
  w <- merge(w,read_annotation(cfg$species_annotation[i],"species"),by="cell_barcode",all.x=TRUE)
  w <- merge(w,read_annotation(cfg$genus_annotation[i],"genus"),by="cell_barcode",all.x=TRUE)
  w[,`:=`(sample=s,crop=cfg$crop[i],cell_id=paste(s,cell_barcode,sep="__"))]
  w
}),fill=TRUE)
cells[is.na(species)|species=="",species:="Unannotated"]
cells[is.na(genus)|genus=="",genus:="Unannotated"]

counts <- rbindlist(lapply(sample_levels,function(s){
  p <- file.path(root,"03_matrix",s,paste0(s,".primary.cell_amp_umi.tsv.gz"))
  if(!file.exists(p)) stop("Missing count file: ",p)
  x <- fread(cmd=paste("gzip -cd",shQuote(p)))
  x[,cell_id:=paste(sample,cell_barcode,sep="__")]; x
}),use.names=TRUE)
counts <- merge(counts,ref[,.(amp_id,primary_source,source_pattern,length)],by="amp_id",all.x=TRUE)
if(anyNA(counts$primary_source)) stop("AMP IDs missing from reference metadata")

cell_expr <- counts[,.(amp_umi=sum(umi_count),n_amp=uniqueN(amp_id)),by=.(cell_id,sample,cell_barcode)]
cells <- merge(cells,cell_expr,by=c("cell_id","sample","cell_barcode"),all.x=TRUE)
cells[is.na(amp_umi),`:=`(amp_umi=0,n_amp=0)]
cells[,`:=`(amp_positive=amp_umi>=1,amp_positive_2umi=amp_umi>=2)]
setcolorder(cells,c("cell_id","sample","crop","cell_barcode","species","genus","amp_positive","amp_positive_2umi","n_amp","amp_umi"))
fwrite(cells,file.path(stats_dir,"cell_amp_metrics.tsv"),sep="\t")

sample_summary <- cells[,.(total_cells=.N,amp_positive_cells=sum(amp_positive),amp_positive_cells_2umi=sum(amp_positive_2umi),
  amp_positive_fraction=mean(amp_positive),amp_positive_fraction_2umi=mean(amp_positive_2umi),
  median_amp_per_positive=median(n_amp[amp_positive]),median_umi_per_positive=median(amp_umi[amp_positive])),by=.(sample,crop)]
sample_summary[,sample:=factor(sample,levels=sample_levels)]; setorder(sample)
fwrite(sample_summary,file.path(stats_dir,"sample_amp_summary.tsv"),sep="\t")

amp_by_sample <- counts[,.(umi=sum(umi_count),positive_cells=uniqueN(cell_id)),by=.(sample,amp_id,primary_source,source_pattern)]
den <- cells[,.(total_cells=.N),by=sample]
amp_by_sample <- merge(amp_by_sample,den,by="sample"); amp_by_sample[,cell_prevalence:=positive_cells/total_cells]
fwrite(amp_by_sample,file.path(stats_dir,"amp_by_sample.tsv.gz"),sep="\t",compress="gzip")
amp_global <- amp_by_sample[,.(total_umi=sum(umi),total_positive_cells=sum(positive_cells),n_samples=uniqueN(sample),
  n_crops=uniqueN(cfg$crop[match(sample,cfg$sample)]),max_cell_prevalence=max(cell_prevalence),
  min_positive_cells_across_detected=min(positive_cells)),by=.(amp_id,primary_source,source_pattern)]
core_ids <- amp_by_sample[positive_cells>=3,.(qualified_samples=uniqueN(sample)),by=amp_id][qualified_samples==8,amp_id]
amp_global[,`:=`(
  core_amp=amp_id %in% core_ids,
  high_frequency_amp=max_cell_prevalence>=0.01 & total_positive_cells>=10,
  sample_restricted_amp=n_samples==1 & total_positive_cells>=5 & total_umi>=10
)]
amp_global[,classification:=fifelse(core_amp,"Core",fifelse(sample_restricted_amp,"Sample-restricted",fifelse(high_frequency_amp,"High-frequency","Other")))]
fwrite(amp_global,file.path(stats_dir,"amp_classification.tsv"),sep="\t")

source_comp <- counts[,.(umi=sum(umi_count),positive_cells=uniqueN(cell_id)),by=.(sample,primary_source)]
source_comp <- merge(source_comp,cfg[,.(sample,crop)],by="sample"); source_comp[,umi_fraction:=umi/sum(umi),by=sample]
fwrite(source_comp,file.path(stats_dir,"source_composition_by_sample.tsv"),sep="\t")
crop_comp <- source_comp[,.(umi=sum(umi),positive_cells=sum(positive_cells)),by=.(crop,primary_source)]; crop_comp[,umi_fraction:=umi/sum(umi),by=crop]
fwrite(crop_comp,file.path(stats_dir,"source_composition_by_crop.tsv"),sep="\t")

tax_counts <- merge(counts,cells[,.(cell_id,species,genus)],by="cell_id")
species_amp <- tax_counts[species!="Unannotated",.(umi=sum(umi_count),positive_cells=uniqueN(cell_id)),by=.(species,amp_id,primary_source)]
fwrite(species_amp,file.path(stats_dir,"species_amp_association.tsv.gz"),sep="\t",compress="gzip")
species_heterogeneity <- cells[,.(cells=.N,positive_cells=sum(amp_positive),positive_fraction=mean(amp_positive),
 mean_amp_umi=mean(amp_umi),variance_amp_umi=var(amp_umi),cv_amp_umi=ifelse(mean(amp_umi)>0,sd(amp_umi)/mean(amp_umi),NA_real_),
 mean_n_amp=mean(n_amp)),by=.(sample,species)]
fwrite(species_heterogeneity,file.path(stats_dir,"species_amp_heterogeneity.tsv"),sep="\t")

# Figure 4A: reference composition before and after 90/90 clustering.
input <- fread(file.path(root,"01_reference/all_amp.input_records.tsv"))[,.(sequences=.N),by=source]
post <- ref[,.(sequences=.N),by=.(source=primary_source)]
dbcomp <- rbind(input[,stage:="Input"],post[,stage:="NR 90/90"])
fwrite(dbcomp,file.path(stats_dir,"reference_composition.tsv"),sep="\t")
pA1 <- ggplot(dbcomp,aes(factor(source,levels=names(cols_source)),sequences,fill=source))+geom_col()+facet_wrap(~stage,scales="free_y")+
 scale_fill_manual(values=cols_source)+scale_y_continuous(labels=scales::comma)+labs(x=NULL,y="AMP sequences",title="Reference database composition")+theme_bw()+theme(legend.position="none",axis.text.x=element_text(angle=30,hjust=1))
patterns <- ref[,.(clusters=.N),by=source_pattern][order(-clusters)][1:min(10,.N)]
pA2 <- ggplot(patterns,aes(reorder(source_pattern,clusters),clusters))+geom_col(fill="#636363")+coord_flip()+scale_y_continuous(labels=scales::comma)+labs(x=NULL,y="NR clusters",title="Top source-membership patterns")+theme_bw()
pdf(file.path(fig_dir,"Figure4A_reference_database.pdf"),11,5); print(pA1); print(pA2); dev.off()
ggsave(file.path(fig_dir,"Figure4A_reference_database.png"),pA1,width=8,height=4,dpi=300)

# Figure 4B.
pB1 <- ggplot(sample_summary,aes(sample,total_cells,fill=crop))+geom_col()+geom_text(aes(label=scales::comma(total_cells)),vjust=-.25,size=3)+
 scale_fill_manual(values=cols_crop)+labs(x=NULL,y="Cells",title="Cells retained by whitelist")+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
pB2 <- ggplot(sample_summary,aes(sample,amp_positive_fraction*100,fill=crop))+geom_col()+geom_text(aes(label=sprintf("%.1f%%",amp_positive_fraction*100)),vjust=-.25,size=3)+
 scale_fill_manual(values=cols_crop)+labs(x=NULL,y="AMP-positive cells (%)",title="AMP-positive fraction (>=1 assigned UMI)")+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="none")
pdf(file.path(fig_dir,"Figure4B_cells_positive_fraction.pdf"),11,5); print(pB1); print(pB2); dev.off()
ggsave(file.path(fig_dir,"Figure4B_cells_positive_fraction.png"),pB2,width=8,height=4,dpi=300)

# Figure 4C.
pos <- cells[amp_positive]; pos[,sample:=factor(sample,levels=sample_levels)]
pC1 <- ggplot(pos,aes(sample,n_amp,fill=crop))+geom_violin(scale="width",trim=TRUE)+geom_boxplot(width=.12,outlier.shape=NA)+scale_y_log10()+scale_fill_manual(values=cols_crop)+labs(x=NULL,y="AMPs per cell (log10)")+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
pC2 <- ggplot(pos,aes(sample,amp_umi,fill=crop))+geom_violin(scale="width",trim=TRUE)+geom_boxplot(width=.12,outlier.shape=NA)+scale_y_log10()+scale_fill_manual(values=cols_crop)+labs(x=NULL,y="AMP UMIs per cell (log10)")+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="none")
pdf(file.path(fig_dir,"Figure4C_single_cell_amp_distributions.pdf"),11,5); print(pC1); print(pC2); dev.off()
ggsave(file.path(fig_dir,"Figure4C_single_cell_amp_distributions.png"),pC2,width=8,height=4,dpi=300)

# Figure 4D.
source_comp[,sample:=factor(sample,levels=sample_levels)]
pD1 <- ggplot(source_comp,aes(sample,umi_fraction,fill=primary_source))+geom_col()+scale_fill_manual(values=cols_source)+scale_y_continuous(labels=scales::percent)+labs(x=NULL,y="AMP UMI composition",fill="Source")+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))
pD2 <- ggplot(crop_comp,aes(factor(crop,levels=crop_levels),umi_fraction,fill=primary_source))+geom_col()+scale_fill_manual(values=cols_source)+scale_y_continuous(labels=scales::percent)+labs(x=NULL,y="AMP UMI composition",fill="Source")+theme_bw()
pdf(file.path(fig_dir,"Figure4D_source_composition.pdf"),11,5); print(pD1); print(pD2); dev.off()
ggsave(file.path(fig_dir,"Figure4D_source_composition.png"),pD1,width=8,height=4,dpi=300)

# Figure 4E: independent AMP UMAP only when information content is adequate.
features <- amp_by_sample[,.(cells=sum(positive_cells)),by=amp_id][cells>=5,amp_id]
umap_ok <- nrow(pos)>=200 && length(features)>=20 && median(pos$n_amp)>=2
decision <- data.table(positive_cells=nrow(pos),features_detected_in_5_cells=length(features),median_amp_per_positive=median(pos$n_amp),independent_umap_supported=umap_ok,
 reason=if(umap_ok) "positive_cells>=200; features>=20; median_nAMP>=2" else "AMP matrix too sparse; full-transcriptome coordinates required")
fwrite(decision,file.path(stats_dir,"umap_decision.tsv"),sep="\t")
if(umap_ok){
  cc <- counts[amp_id %in% features]
  ci <- match(cc$cell_id,pos$cell_id); fi <- match(cc$amp_id,features)
  m <- sparseMatrix(i=fi,j=ci,x=cc$umi_count,dims=c(length(features),nrow(pos)),dimnames=list(features,pos$cell_id))
  md <- as.data.frame(pos[match(colnames(m),cell_id)])
  rownames(md) <- md$cell_id
  seu <- CreateSeuratObject(counts=m,meta.data=md)
  seu <- NormalizeData(seu,verbose=FALSE); seu <- FindVariableFeatures(seu,nfeatures=min(2000,nrow(m)),verbose=FALSE)
  seu <- ScaleData(seu,verbose=FALSE); npcs <- min(30,nrow(m)-1,ncol(m)-1)
  seu <- RunPCA(seu,npcs=npcs,verbose=FALSE); dims_use <- 1:min(20,npcs)
  seu <- FindNeighbors(seu,dims=dims_use,verbose=FALSE); seu <- FindClusters(seu,resolution=.5,verbose=FALSE); seu <- RunUMAP(seu,dims=dims_use,seed.use=20260719,verbose=FALSE)
  emb <- as.data.table(Embeddings(seu,"umap"),keep.rownames="cell_id"); setnames(emb,c("cell_id","UMAP_1","UMAP_2"))
  em <- merge(emb,pos,by="cell_id"); em[,cluster:=as.character(Idents(seu)[em$cell_id])]
  fwrite(em,file.path(stats_dir,"amp_umap_coordinates.tsv"),sep="\t")
  pE1 <- ggplot(em,aes(UMAP_1,UMAP_2,color=sample))+geom_point(size=.25,alpha=.7)+labs(title="AMP-expression UMAP by sample")+theme_void()
  pE2 <- ggplot(em,aes(UMAP_1,UMAP_2,color=log1p(amp_umi)))+geom_point(size=.25,alpha=.8)+scale_color_viridis_c()+labs(color="log1p AMP UMI",title="AMP expression intensity")+theme_void()
  pdf(file.path(fig_dir,"Figure4E_amp_umap.pdf"),11,5); print(pE1); print(pE2); dev.off()
  ggsave(file.path(fig_dir,"Figure4E_amp_umap.png"),pE2,width=6,height=5,dpi=300)
  cluster_heterogeneity <- em[,.(cells=.N,mean_amp_umi=mean(amp_umi),median_amp_umi=median(amp_umi),mean_n_amp=mean(n_amp)),by=.(cluster,sample)]
  fwrite(cluster_heterogeneity,file.path(stats_dir,"amp_cluster_heterogeneity.tsv"),sep="\t")
}

# Figure 4F: sample-level logCPM heatmap for classified AMPs.
selected <- amp_global[classification!="Other"][order(-core_amp,-high_frequency_amp,-total_umi),head(amp_id,60)]
if(length(selected)<20) selected <- amp_global[order(-total_umi),head(amp_id,60)]
hm <- amp_by_sample[amp_id %in% selected,.(umi=sum(umi)),by=.(amp_id,sample)]
hm <- merge(CJ(amp_id=selected,sample=sample_levels),hm,all.x=TRUE); hm[is.na(umi),umi:=0]
hm[,logCPM:=log1p(1e6*umi/sum(umi)),by=sample]
hm <- merge(hm,amp_global[,.(amp_id,classification,primary_source)],by="amp_id")
hm[,amp_label:=paste(amp_id,classification,primary_source,sep=" | ")]
pF <- ggplot(hm,aes(factor(sample,levels=sample_levels),factor(amp_label,levels=rev(unique(amp_label))),fill=logCPM))+geom_tile()+scale_fill_viridis_c()+labs(x=NULL,y=NULL,fill="log1p CPM",title="Core, high-frequency and sample-restricted AMPs")+theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1),axis.text.y=element_text(size=5))
ggsave(file.path(fig_dir,"Figure4F_amp_heatmap.pdf"),pF,width=9,height=12)
ggsave(file.path(fig_dir,"Figure4F_amp_heatmap.png"),pF,width=9,height=12,dpi=300)

# Supplemental species-AMP association dot plot.
top_species <- species_amp[,.(umi=sum(umi)),by=species][order(-umi),head(species,15)]
top_amp <- species_amp[,.(umi=sum(umi)),by=amp_id][order(-umi),head(amp_id,20)]
assoc <- species_amp[species %in% top_species & amp_id %in% top_amp]
pS <- ggplot(assoc,aes(amp_id,species,size=positive_cells,color=log1p(umi)))+geom_point(alpha=.8)+scale_color_viridis_c()+labs(x="AMP",y="Potential source species",size="Positive cells",color="log1p UMI",title="AMP–species associations")+theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1,size=7))
ggsave(file.path(fig_dir,"FigureS_AMP_species_association.pdf"),pS,width=11,height=7)
ggsave(file.path(fig_dir,"FigureS_AMP_species_association.png"),pS,width=11,height=7,dpi=300)

writeLines(capture.output(sessionInfo()),file.path(root,"logs/04_R_sessionInfo.txt"))
