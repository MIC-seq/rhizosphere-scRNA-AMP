.libPaths(c('/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/software/Rlib',.libPaths()))
suppressPackageStartupMessages(library(data.table))
b<-'/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/fg_MES_20261027';o<-file.path(b,'g_statistics_v6');dir.create(o,recursive=TRUE,showWarnings=FALSE)
P<-fread(file.path(b,'tables/all_soil_state_pairs.tsv'));P<-P[n_cells_S1>=5 & n_cells_S2>=5];prev<-fread(file.path(b,'tables/all18_state_tests_sensitivity.tsv'))[minimum_cells==5]
focus<-fread(file.path(b,'tables/replacement_g_selection.tsv'))[,.(species,class_id)];crops<-c('soybean','rice','wheat')
# Frozen plan before new tests: independent soils; all six classes and three taxa.
# Family 1: original 18 overall paired state contrasts, exact sign-flip unchanged.
# Family 2: all within-crop state contrasts with at least 3 paired soils, BH eligible tests.
# Family 3: all crop pairs with at least 2 soils/group, for S1, S2 and delta; BH eligible tests.
flip<-function(d){n<-length(d);if(n<3)return(NA_real_);s<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));mean(abs(drop(s%*%d)/n)>=abs(mean(d))-1e-12)}
perm<-function(a,b){n<-length(a);v<-c(a,b);sets<-combn(length(v),n);obs<-mean(a)-mean(b);stats<-apply(sets,2,function(ix)mean(v[ix])-mean(v[-ix]));list(P=mean(abs(stats)>=abs(obs)-1e-12),delta=obs,permutations=length(stats))}
O<-P[,.(n_soils=.N,mean_delta=mean(delta_MES),P=flip(delta_MES)),by=.(species,class_id)];O[,q_BH:=p.adjust(P,'BH',n=18)];check<-merge(O,prev,by=c('species','class_id'));stopifnot(nrow(check)==18,max(abs(check$P-check$P_MES))<1e-12,max(abs(check$q_BH-check$q_MES_BH))<1e-12)
W<-list();B<-list();keys<-unique(P[,.(species,class_id)])
for(i in seq_len(nrow(keys))){z<-P[species==keys$species[i]&class_id==keys$class_id[i]]
 for(cr in crops){d<-z[crop==cr,delta_MES];W[[length(W)+1]]<-data.table(species=keys$species[i],class_id=keys$class_id[i],crop=cr,n_soils=length(d),mean_delta=if(length(d))mean(d)else NA_real_,P=flip(d),status=if(length(d)>=3)'tested'else'insufficient_paired_soils')}
 for(pair in combn(crops,2,simplify=FALSE))for(ep in c('MES_S1','MES_S2','delta_MES')){
  a<-z[crop==pair[1]][[ep]];bb<-z[crop==pair[2]][[ep]];ok<-length(a)>=2&&length(bb)>=2
  v<-if(ok)perm(a,bb)else list(P=NA_real_,delta=if(length(a)&&length(bb))mean(a)-mean(bb)else NA_real_,permutations=NA_integer_)
  B[[length(B)+1]]<-data.table(species=keys$species[i],class_id=keys$class_id[i],crop_A=pair[1],crop_B=pair[2],endpoint=ep,n_A=length(a),n_B=length(bb),mean_A=if(length(a))mean(a)else NA_real_,mean_B=if(length(bb))mean(bb)else NA_real_,delta_A_minus_B=v$delta,P=v$P,n_permutations=v$permutations,status=if(ok)'tested'else'insufficient_independent_soils')
 }
}
W<-rbindlist(W);B<-rbindlist(B);W[,q_BH:=NA_real_];B[,q_BH:=NA_real_];W[status=='tested',q_BH:=p.adjust(P,'BH')];B[status=='tested',q_BH:=p.adjust(P,'BH')]
fwrite(O,file.path(o,'overall_state_all18.tsv'),sep='\t');fwrite(W,file.path(o,'within_crop_state_all54.tsv'),sep='\t');fwrite(B,file.path(o,'between_crop_all162.tsv'),sep='\t')
for(n in c('O','W','B')){d<-merge(get(n),focus,by=c('species','class_id'));fwrite(d,file.path(o,paste0('focus_',n,'.tsv')),sep='\t')}
summary<-rbind(data.table(family='Overall state',eligible=nrow(O),nominal_significant=sum(O$P<.05),FDR_significant=sum(O$q_BH<.05)),data.table(family='Within-crop state',eligible=sum(W$status=='tested'),nominal_significant=sum(W$P<.05,na.rm=TRUE),FDR_significant=sum(W$q_BH<.05,na.rm=TRUE)),data.table(family='Between crops',eligible=sum(B$status=='tested'),nominal_significant=sum(B$P<.05,na.rm=TRUE),FDR_significant=sum(B$q_BH<.05,na.rm=TRUE)))
fwrite(summary,file.path(o,'test_family_summary.tsv'),sep='\t');capture.output(sessionInfo(),file=file.path(o,'stats_sessionInfo.txt'));print(summary);print(merge(B[status=='tested'],focus,by=c('species','class_id'))[,.(species,endpoint,n_A,n_B,delta_A_minus_B,P,q_BH)]);cat('STATE_CROP_TESTS_COMPLETE\n')