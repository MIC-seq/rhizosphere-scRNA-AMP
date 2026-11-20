.libPaths(c('/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/software/Rlib',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix)})
b <- '/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions';o <- file.path(b,'fg_MES_20261027');prior<-file.path(b,'six_class_joint_mes_20261026');author<-file.path(b,'imeta_mes_20261025/author_code_only')
for(d in c('scripts','tables','figures','reports','source_data','provenance','logs'))dir.create(file.path(o,d),recursive=TRUE,showWarnings=FALSE)
source(file.path(author,'scripts/author_functions_unmodified.R'))
S<-fread(file.path(prior,'tables/six_class_cell_both_scores.tsv'));K<-fread(file.path(b,'inputs/state_key.tsv'))
S<-merge(S,K[,.(species,functional_state,state_code)],by=c('species','functional_state'),all.x=TRUE,sort=FALSE);stopifnot(nrow(S)==1708*6,!anyNA(S$state_code))
sp<-c('Chloracidobacterium validum','Halomonas sp. JS92-SW72','Pseudomonas aeruginosa');fo<-data.table(species=sp,class_id=c('energy','amino_acid','nucleotide')); classes<-unique(S$class_id)
soil<-states<-audit<-list()
for(i in seq_along(sp)){
 z<-readRDS(file.path(prior,'objects',paste0('species_',i,'_six_classes.rds')));M<-as.matrix(z$raw)
 C<-S[species==sp[i]&class_id==classes[1]];C<-C[match(z$cell_ids,cell_id)];stopifnot(identical(C$cell_id,z$cell_ids),ncol(M)==nrow(C))
 n<-author_normalize(M,z$random_scores);audit[[i]]<-data.table(species=sp[i],max_existing_pooled_MES_difference=max(abs(n$MES[,1]-z$norm$MES[,1])))
 for(st in unique(C$state_code)){
  j<-which(C$state_code==st);v<-author_normalize(M[,j,drop=FALSE],lapply(z$random_scores,function(q)q[,j,drop=FALSE]))
  states[[paste(i,st)]]<-data.table(species=sp[i],state_code=st,class_id=rownames(M),n_cells=length(j),MES=v$MES[,1],P_background=v$P[,1])
 }
 for(sa in unique(C$sample))for(st in unique(C$state_code)){
  j<-which(C$sample==sa & C$state_code==st);if(!length(j))next
  v<-author_normalize(M[,j,drop=FALSE],lapply(z$random_scores,function(q)q[,j,drop=FALSE]))
  for(k in rownames(M)){
   q<-S[species==sp[i]&class_id==k&sample==sa&state_code==st]
   soil[[paste(i,sa,st,k)]]<-data.table(species=sp[i],sample=sa,crop=C$crop[j[1]],state_code=st,class_id=k,class=unique(q$class),n_cells=length(j),MES=v$MES[k,1],P_background=v$P[k,1],raw_mean=mean(as.numeric(M[k,j])),pathway_UMI=sum(q$pathway_UMI),total_UMI=sum(q$retained_gene_umi),detected_cells=sum(q$pathway_UMI>0),UMI_fraction_pp=100*sum(q$pathway_UMI)/sum(q$retained_gene_umi))
  }
 }
}
T<-rbindlist(soil);ST<-rbindlist(states);A<-rbindlist(audit);stopifnot(max(A$max_existing_pooled_MES_difference)<1e-12)
P<-merge(T[state_code=='S1'],T[state_code=='S2'],by=c('species','sample','crop','class_id','class'),suffixes=c('_S1','_S2'))
P[,`:=`(delta_MES=MES_S2-MES_S1,delta_UMI_pp=UMI_fraction_pp_S2-UMI_fraction_pp_S1,detected_both=pathway_UMI_S1>0&pathway_UMI_S2>0)]
flip<-function(d){if(length(d)<3)return(NA_real_);if(all(abs(d)<1e-12))return(1);x<-as.matrix(expand.grid(rep(list(c(-1,1)),length(d))));mean(abs(drop(x%*%d)/length(d))>=abs(mean(d))-1e-12)}
Tests<-list()
for(m in c(3L,5L,10L)){
 Q<-P[n_cells_S1>=m&n_cells_S2>=m]
 v<-Q[,.(n_soils=.N,mean_delta_MES=mean(delta_MES),min_delta_MES=min(delta_MES),max_delta_MES=max(delta_MES),positive_MES=sum(delta_MES>0),negative_MES=sum(delta_MES<0),P_MES=flip(delta_MES),mean_delta_UMI_pp=mean(delta_UMI_pp),positive_UMI=sum(delta_UMI_pp>0),negative_UMI=sum(delta_UMI_pp<0),P_UMI=flip(delta_UMI_pp),direction_agreement=sum(sign(delta_MES)==sign(delta_UMI_pp)),n_zero_detection_pairs=sum(!detected_both)),by=.(species,class_id,class)]
 v[,`:=`(minimum_cells=m,q_MES_BH=p.adjust(P_MES,'BH',n=18),q_UMI_BH=p.adjust(P_UMI,'BH',n=18))];Tests[[as.character(m)]]<-v
}
Tests<-rbindlist(Tests);Q<-P[n_cells_S1>=5&n_cells_S2>=5];F<-merge(Q,fo,by=c('species','class_id'));FS<-merge(S,fo,by=c('species','class_id'))
C<-F[,.(n=.N,mean_delta_MES=mean(delta_MES),lo=min(delta_MES),hi=max(delta_MES),mean_delta_UMI_pp=mean(delta_UMI_pp),positive_MES=sum(delta_MES>0),positive_UMI=sum(delta_UMI_pp>0),zero_detection_pairs=sum(!detected_both)),by=.(species,class_id,class,crop)]
C<-merge(merge(CJ(species=sp,crop=c('soybean','rice','wheat')),fo,by='species'),C,by=c('species','class_id','crop'),all.x=TRUE);C[is.na(n),n:=0L]
# Compare to the existing Fig5g UMI results exactly on the paired-soil scale.
old<-fread(file.path(b,'broad_metabolism_20261018/tables/focus_soil_contrasts.tsv'))
B<-merge(F,old[,.(species,sample,crop,old_UMI_delta=delta_pp)],by=c('species','sample','crop'));stopifnot(nrow(B)==14,max(abs(B$delta_UMI_pp-B$old_UMI_delta))<1e-10)
H<-fread(file.path(author,'tables/state_MES.tsv'));sel<-unique(fread(file.path(author,'tables/Fig_author_MES_state_24pathways_source.tsv'))[,.(pathway_id,pathway,group,row)]);setorder(sel,row);stopifnot(nrow(sel)==24,nrow(H)==315)
H[,q_background_BH:=p.adjust(P,'BH')];H[,`:=`(display_MES=ifelse(pathway_UMI==0,0,MES),detected=pathway_UMI>0)]
HD<-merge(H[unit=='S1'],H[unit=='S2'],by=c('species','pathway_id','pathway','group'),suffixes=c('_S1','_S2'));HD[,`:=`(delta_pooled_MES=MES_S2-MES_S1,detected_both=detected_S1&detected_S2)]
write<-function(x,n,sub='tables')fwrite(x,file.path(o,sub,paste0(n,'.tsv')),sep='\t',na='NA')
write(T,'six_class_soil_state_MES');write(ST,'six_class_state_MES');write(A,'score_identity_QA');write(P,'all_soil_state_pairs');write(Tests,'all18_state_tests_sensitivity');write(F,'g_paired_soils');write(FS,'g_UMAP_cells');write(C,'g_crop_contrasts');write(H,'f_all45_author_MES');write(sel,'f_fixed_24_pathways');write(HD,'f_pooled_state_differences');write(B,'old_UMI_identity_QA')
write(H[pathway_id%in%sel$pathway_id],'f_heatmap','source_data');write(FS,'g_cells','source_data');write(F,'g_soil_pairs','source_data');write(C,'g_crop_points','source_data')
capture.output(sessionInfo(),file=file.path(o,'reports/analysis_sessionInfo.txt'));print(merge(Tests[minimum_cells==5],fo,by=c('species','class_id')));print(C);cat('FG_FIXED_METHOD_ANALYSIS_COMPLETE\n')
