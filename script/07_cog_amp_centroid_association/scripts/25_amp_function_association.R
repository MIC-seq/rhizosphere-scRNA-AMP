#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(Matrix); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop("config matrix_root embedding_root output_root")
cfg <- fromJSON(args[1], simplifyVector=TRUE)
matrix_root <- args[2]; embedding_root <- args[3]; output_root <- args[4]
dir.create(output_root, recursive=TRUE, showWarnings=FALSE)
slug <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
read_sparse <- function(path) { con <- gzfile(path,"rt"); on.exit(close(con)); as(readMM(con),"dgCMatrix") }

all_continuous <- list(); all_state_amp <- list(); all_profiles <- list()
for (species in cfg$species) {
  s <- slug(species); input <- file.path(matrix_root,s); embed <- file.path(embedding_root,s)
  meta <- fread(file.path(input,"cells.tsv")); cogs <- fread(file.path(input,"cogs.tsv"))$COG_category
  mat <- read_sparse(file.path(input,"cog_by_cell.mtx.gz")); rownames(mat) <- cogs; colnames(mat) <- meta$cell_id
  state <- fread(file.path(embed,"cog_umap_cells.tsv"))[,.(cell_id,seurat_clusters)]
  meta <- merge(meta,state,by="cell_id",all=FALSE); mat <- mat[,meta$cell_id,drop=FALSE]
  totals <- Matrix::colSums(mat); norm <- log1p(t(t(mat)/pmax(totals,1))*10000)
  scores <- as.data.table(t(as.matrix(norm)),keep.rownames="cell_id")
  dat <- merge(meta,scores,by="cell_id",all=FALSE)
  dat[,`:=`(amp_load=as.numeric(gene_amp_load), log_umi=log1p(as.numeric(gene_total_non_rrna_umi)),
             detected=as.numeric(gene_detected_non_rrna_genes))]

  per_sample <- list()
  for (sample_id in unique(dat$sample)) {
    d <- dat[sample==sample_id & is.finite(amp_load)]
    if (nrow(d) < cfg$association$min_cells_per_sample_species || sd(d$amp_load)==0) next
    for (cog in cogs) {
      if (sd(d[[cog]])==0) next
      fit <- tryCatch(lm(scale(d[[cog]]) ~ scale(amp_load) + scale(log_umi) + scale(detected),data=d),error=function(e)NULL)
      if (is.null(fit)) next
      co <- summary(fit)$coefficients
      if (!("scale(amp_load)" %in% rownames(co)) || !is.finite(co["scale(amp_load)","Std. Error"])) next
      per_sample[[length(per_sample)+1]] <- data.table(species=species,sample=sample_id,COG_category=cog,
        beta=co["scale(amp_load)","Estimate"],se=co["scale(amp_load)","Std. Error"],
        pvalue=co["scale(amp_load)","Pr(>|t|)"],cells=nrow(d))
    }
  }
  per_sample <- if(length(per_sample)) rbindlist(per_sample) else data.table()
  if (nrow(per_sample)) {
    meta_rows <- per_sample[,{
      ok <- is.finite(beta)&is.finite(se)&se>0
      b <- beta[ok]; ss <- se[ok]; w <- 1/(ss^2)
      if(length(b)<cfg$association$min_samples) list(samples=length(b),cells=sum(cells[ok]),beta=NA_real_,se=NA_real_,z=NA_real_,pvalue=NA_real_,positive_samples=sum(b>0),negative_samples=sum(b<0))
      else { bb<-sum(w*b)/sum(w); see<-sqrt(1/sum(w)); list(samples=length(b),cells=sum(cells[ok]),beta=bb,se=see,z=bb/see,pvalue=2*pnorm(-abs(bb/see)),positive_samples=sum(b>0),negative_samples=sum(b<0)) }
    },by=.(species,COG_category)]
    meta_rows[,FDR:=p.adjust(pvalue,"BH"),by=species]
    all_continuous[[length(all_continuous)+1]] <- meta_rows
    fwrite(per_sample,file.path(output_root,paste0(s,"_continuous_amp_cog_per_sample.tsv")),sep="\t")
  }

  profiles <- dat[,lapply(.SD,mean),by=seurat_clusters,.SDcols=cogs]
  profiles[,species:=species]; all_profiles[[length(all_profiles)+1]] <- profiles

  clusters <- sort(unique(dat$AMP_cluster)); states <- sort(unique(dat$seurat_clusters))
  assoc <- list()
  for (amp_cluster in clusters) for (functional_state in states) {
    per <- list()
    for (sample_id in unique(dat$sample)) {
      d <- dat[sample==sample_id]
      if (nrow(d)<cfg$association$min_cells_per_sample_species) next
      d[,`:=`(cluster_bin=as.integer(AMP_cluster==amp_cluster),state_bin=as.integer(seurat_clusters==functional_state))]
      if(sum(d$cluster_bin)<cfg$association$minimum_amp_cluster_cells || uniqueN(d$cluster_bin)<2 || uniqueN(d$state_bin)<2) next
      fit<-tryCatch(glm(cluster_bin~state_bin+scale(log_umi)+scale(detected),family=binomial(),data=d),error=function(e)NULL)
      if(is.null(fit)||!"state_bin"%in%rownames(coef(summary(fit)))) next
      z<-coef(summary(fit))["state_bin",];if(!all(is.finite(z[c("Estimate","Std. Error")]))||z["Std. Error"]<=0)next
      per[[length(per)+1]] <- c(lor=unname(z["Estimate"]),se=unname(z["Std. Error"]),cells=nrow(d))
    }
    if(length(per)<cfg$association$min_samples) next
    x<-do.call(rbind,per); w<-1/x[,"se"]^2; bb<-sum(w*x[,"lor"])/sum(w); see<-sqrt(1/sum(w))
    assoc[[length(assoc)+1]] <- data.table(species=species,AMP_cluster=amp_cluster,functional_state=functional_state,
      samples=nrow(x),cells=sum(x[,"cells"]),log_odds_ratio=bb,se=see,pvalue=2*pnorm(-abs(bb/see)))
  }
  if(length(assoc)) { aa<-rbindlist(assoc); aa[,FDR:=p.adjust(pvalue,"BH")]; all_state_amp[[length(all_state_amp)+1]]<-aa }
}

continuous <- if(length(all_continuous)) rbindlist(all_continuous,fill=TRUE) else data.table()
state_amp <- if(length(all_state_amp)) rbindlist(all_state_amp,fill=TRUE) else data.table()
profiles <- rbindlist(all_profiles,fill=TRUE)
fwrite(continuous,file.path(output_root,"continuous_amp_cog_meta_analysis.tsv"),sep="\t")
fwrite(state_amp,file.path(output_root,"amp_cluster_functional_state_meta_analysis.tsv"),sep="\t")
fwrite(profiles,file.path(output_root,"functional_state_cog_profiles.tsv"),sep="\t")
cat("continuous rows",nrow(continuous),"significant",sum(continuous$FDR<cfg$association$fdr,na.rm=TRUE),"\n")
cat("state-AMP rows",nrow(state_amp),"significant",sum(state_amp$FDR<cfg$association$fdr,na.rm=TRUE),"\n")
