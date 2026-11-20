#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Matrix);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=3)stop("root config outdir")
root<-normalizePath(args[1]);cfg<-fromJSON(args[2]);out<-args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(as.integer(cfg$embedding$random_seed))
amp_root<-cfg$inputs$amp_matrix_root
cells<-fread(cmd=paste("gzip -cd",shQuote(file.path(root,"results/12_integrated_amp_function/figure5_optimized/data/figure5_cells_with_state_program.tsv.gz"))))
cells[,cell_barcode:=sub("^[^_]+__","",cell_id)]
targets<-cells[,.(sample,cell_barcode,cell_id)]
files<-list.files(amp_root,pattern="primary.cell_amp_umi.tsv.gz$",recursive=TRUE,full.names=TRUE)
long<-rbindlist(lapply(files,function(f){
 d<-fread(cmd=paste("gzip -cd",shQuote(f)),select=c("sample","cell_barcode","amp_id","umi_count"))
 merge(d,targets,by=c("sample","cell_barcode"),all=FALSE)
}))
long<-long[,.(umi_count=sum(umi_count)),by=.(cell_id,amp_id)]
cell_levels<-cells$cell_id;amp_levels<-sort(unique(long$amp_id))
X<-sparseMatrix(i=match(long$cell_id,cell_levels),j=match(long$amp_id,amp_levels),x=long$umi_count,dims=c(length(cell_levels),length(amp_levels)),dimnames=list(cell_levels,amp_levels))
tot<-rowSums(X);Xn<-X;Xn@x<-log1p(X@x*1e4/pmax(tot[X@i+1L],1))

predict_loso<-function(meta,Xmat,labels,nfeat){
 pred<-rep(NA_character_,nrow(meta));score_out<-rep(NA_real_,nrow(meta));selected<-list()
 for(heldout in unique(meta$sample)){
  tr<-which(meta$sample!=heldout);te<-which(meta$sample==heldout)
  y<-labels[tr];states<-sort(unique(y));if(length(states)<2)next
  keep<-which(Matrix::colSums(Xmat[tr,,drop=FALSE]>0)>=5)
  if(length(keep)<2)next
  cent0<-do.call(rbind,lapply(states,function(st)Matrix::colMeans(Xmat[tr[y==st],keep,drop=FALSE])))
  between<-apply(cent0,2,var);ord<-order(between,decreasing=TRUE,na.last=NA);use<-keep[head(ord,min(nfeat,length(ord)))]
  A<-as.matrix(Xmat[tr,use,drop=FALSE]);B<-as.matrix(Xmat[te,use,drop=FALSE])
  mu<-colMeans(A);ss<-apply(A,2,sd);ss[!is.finite(ss)|ss==0]<-1
  A<-sweep(sweep(A,2,mu,"-"),2,ss,"/");B<-sweep(sweep(B,2,mu,"-"),2,ss,"/")
  cent<-do.call(rbind,lapply(states,function(st)colMeans(A[y==st,,drop=FALSE])));rownames(cent)<-states
  distmat<-sapply(states,function(st)rowSums((B-matrix(cent[st,],nrow(B),ncol(B),byrow=TRUE))^2))
  if(is.null(dim(distmat)))distmat<-matrix(distmat,ncol=length(states))
  best<-max.col(-distmat,ties.method="first");pred[te]<-states[best]
  sorted<-t(apply(distmat,1,sort));score_out[te]<-if(ncol(sorted)>1)sorted[,2]-sorted[,1] else NA_real_
  selected[[heldout]]<-data.table(heldout_sample=heldout,amp_id=colnames(Xmat)[use],feature_rank=seq_along(use),between_state_variance=between[match(use,keep)])
 }
 list(pred=pred,margin=score_out,selected=rbindlist(selected))
}
balanced_accuracy<-function(truth,pred){
 ok<-!is.na(pred);mean(vapply(unique(truth[ok]),function(st)mean(pred[ok][truth[ok]==st]==st),numeric(1)))
}
nfeat<-as.integer(cfg$association$amp_state_prediction_features);B<-as.integer(cfg$association$amp_state_prediction_permutations)
all_pred<-list();all_metric<-list();all_sel<-list();all_conf<-list()
for(sp in unique(cells$species)){
 idx<-which(cells$species==sp);m<-cells[idx];xs<-Xn[idx,,drop=FALSE];truth<-m$functional_state
 obs<-predict_loso(m,xs,truth,nfeat);ba<-balanced_accuracy(truth,obs$pred)
 perm<-rep(NA_real_,B)
 for(b in seq_len(B)){
  yp<-truth;for(ii in split(seq_along(yp),m$sample))yp[ii]<-sample(yp[ii])
  pp<-predict_loso(m,xs,yp,nfeat)$pred;perm[b]<-balanced_accuracy(yp,pp)
 }
 pv<-(1+sum(perm>=ba,na.rm=TRUE))/(1+sum(is.finite(perm)))
 pr<-copy(m);pr[,`:=`(AMP_predicted_state=obs$pred,prediction_margin=obs$margin,prediction_correct=as.integer(functional_state==obs$pred))]
 all_pred[[sp]]<-pr;all_metric[[sp]]<-data.table(species=sp,n_cells=nrow(m),states=uniqueN(truth),balanced_accuracy=ba,permutations=B,permutation_p=pv,permuted_mean=mean(perm,na.rm=TRUE),permuted_sd=sd(perm,na.rm=TRUE),features_per_fold=nfeat)
 cf<-as.data.table(table(truth,obs$pred));setnames(cf,c("observed_state","predicted_state","N"));cf[,species:=sp];all_conf[[sp]]<-cf
 z<-obs$selected;z[,species:=sp];all_sel[[sp]]<-z
}
fwrite(rbindlist(all_pred),file.path(out,"AMP_feature_LOSO_state_predictions.tsv.gz"),sep="\t",compress="gzip")
fwrite(rbindlist(all_metric),file.path(out,"AMP_feature_LOSO_prediction_metrics.tsv"),sep="\t")
fwrite(rbindlist(all_conf),file.path(out,"AMP_feature_LOSO_confusion.tsv"),sep="\t")
fwrite(rbindlist(all_sel),file.path(out,"AMP_feature_LOSO_selected_features.tsv.gz"),sep="\t",compress="gzip")
fwrite(long,file.path(out,"target_cells_individual_AMP_counts.tsv.gz"),sep="\t",compress="gzip")
print(rbindlist(all_metric))
