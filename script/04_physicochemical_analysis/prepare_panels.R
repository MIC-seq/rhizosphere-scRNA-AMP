suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid);library(scales)})
setDTthreads(2)
B <- '/public3/labmember/qianqh/soil_amp_analysis/AMP_four_model_consensus_out'
O <- '/public3/labmember/qianqh/soil_amp_analysis/13_model_test'; D <- file.path(O,'source_data')
dir.create(file.path(O,'figures'),showWarnings=FALSE);dir.create(file.path(O,'panels'),showWarnings=FALSE)
Sys.setenv(FONTCONFIG_FILE=file.path(O,'work/fonts.conf'))
FONT <- 'Arial'
publication_safe_fallback <- 'Arial'
blue <- '#387DA6'; teal <- '#23998D'; gold <- '#D9A349'; purple <- '#8D7AB5'; ink <- '#293B4B'; grey <- '#AAB4BC'
model_cols <- c('ESM2-LGBM'=blue,'LightGBM'=gold,'CNN'='#84B5A1','LSTM'=purple,'Transformer'='#8EBFD2','Attention'='#687787','Macrel'='#A9A9A9')
nov_cols <- c('Database-unmatched'=teal,'Distant match'=gold,'Known-like'=ink)
tier_cols <- c('Conservative'=blue,'Expanded-only'=purple,'Single-model only'='#CFD6DB')
theme_set(theme_classic(base_size=6.5,base_family="Arial")+theme(axis.line=element_line(linewidth=.28,colour=ink),axis.ticks=element_line(linewidth=.25),axis.ticks.length=unit(1,'mm'),axis.text=element_text(colour=ink,size=5.7),axis.title=element_text(size=6.2),plot.title=element_text(size=7,face='bold',colour=ink),plot.subtitle=element_text(size=5.7,colour='#576674',margin=margin(b=3)),plot.caption=element_text(size=5.3,hjust=0,colour='#576674'),legend.title=element_blank(),legend.text=element_text(size=5.5),legend.key.size=unit(2.5,'mm'),legend.spacing.x=unit(1,'mm'),strip.background=element_blank(),strip.text=element_text(size=6,face='bold'),plot.margin=margin(4,5,4,5)))
read <- function(n) fread(file.path(D,n))
write <- function(x,n) fwrite(x,file.path(D,n),sep='\t')
savefig <- function(p,n,w,h,where='figures',tif=TRUE) {
 base <- file.path(O,where,n)
 cairo_pdf(paste0(base,'.pdf'),width=w/25.4,height=h/25.4,family=FONT);if(is.function(p))p() else print(p);dev.off()
 if(requireNamespace('svglite',quietly=TRUE)) svglite::svglite(paste0(base,'.svg'),width=w/25.4,height=h/25.4) else svg(paste0(base,'.svg'),width=w/25.4,height=h/25.4,family=FONT)
 if(is.function(p))p() else print(p);dev.off()
 png(paste0(base,'.png'),width=w,height=h,units='mm',res=300,type='cairo');if(is.function(p))p() else print(p);dev.off()
 if(tif){tiff(paste0(base,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw',type='cairo');if(is.function(p))p() else print(p);dev.off()}
}
pred <- read('benchmark_predictions.tsv'); stopifnot(nrow(pred)==5082,sum(pred$label)==2496)
curves <- rbindlist(lapply(names(model_cols),function(m){
 z <- data.table(y=pred$label,p=pred[[m]])[,.(pos=sum(y),neg=sum(1-y)),by=p][order(-p)]
 z[,`:=`(tp=cumsum(pos),fp=cumsum(neg))];z[,`:=`(tpr=tp/2496,fpr=fp/2586,precision=tp/(tp+fp))]
 rbind(data.table(model=m,fpr=0,tpr=0,precision=1),z[,.(model=m,fpr,tpr,precision)])
}))
write(curves,'Fig2bc_ROC_PR_coordinates.tsv')
metrics <- read('fig2_seven_model_metrics.tsv');metrics[,short:=c('ESM2-LGBM','LightGBM','LSTM','CNN','Transformer','Attention','Macrel')]
auc_check <- curves[,.(auc=sum(diff(fpr)*(head(tpr,-1)+tail(tpr,-1))/2),ap=sum(diff(tpr)*tail(precision,-1))),by=model]
check <- merge(auc_check,metrics,by.x='model',by.y='short');stopifnot(max(abs(check$auc-check$AUROC))<1e-6,max(abs(check$ap-check$AUPRC))<1e-6)
write(check[,.(model,auc,ap,AUROC,AUPRC)],'QA_recomputed_benchmark.tsv')
pB <- ggplot(curves,aes(fpr,tpr,colour=model))+geom_abline(slope=1,intercept=0,linewidth=.25,linetype=2,colour='#C3CAD0')+geom_path(linewidth=.43)+geom_path(data=curves[model=='ESM2-LGBM'],linewidth=.7)+scale_colour_manual(values=model_cols)+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+scale_x_continuous(breaks=c(0,.5,1))+scale_y_continuous(breaks=c(0,.5,1))+labs(title='Benchmark discrimination',subtitle='2,496 AMP + 2,586 non-AMP sequences',x='False-positive rate',y='True-positive rate')+theme(legend.position='none')
lg <- data.table(model=names(model_cols),x=.38,y=seq(.51,.09,length.out=7))
pB <- pB+geom_segment(data=lg,aes(x=x,xend=x+.09,y=y,yend=y,colour=model),inherit.aes=FALSE,linewidth=.6)+geom_text(data=lg,aes(x=x+.12,y=y,label=model,colour=model),inherit.aes=FALSE,hjust=0,size=1.75)
pC <- ggplot(curves,aes(tpr,precision,colour=model))+geom_hline(yintercept=2496/5082,linetype=2,linewidth=.25,colour='#B5BEC5')+geom_path(linewidth=.43)+geom_path(data=curves[model=='ESM2-LGBM'],linewidth=.7)+scale_colour_manual(values=model_cols)+scale_x_continuous(breaks=c(0,.5,1))+scale_y_continuous(breaks=c(0,.5,1))+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+labs(title='Precision-recall performance',subtitle='ESM2-LGBM: AUROC 0.9904; AP 0.9903',x='Recall',y='Precision')+annotate('text',x=.05,y=.18,label='ESM2-LGBM\nPrecision 0.950 | Recall 0.953\nF1 0.951 | MCC 0.904',hjust=0,size=1.95,colour=blue)+theme(legend.position='none')
counts <- data.table(model=c('ESM2-LGBM','Macrel','AMPlify','AMP Scanner v2'),n=c(537,343,4806,3819));counts[,fraction:=n/18454];counts[,model:=factor(model,levels=rev(model))];write(counts,'Fig2d_model_yield.tsv')
pD <- ggplot(counts,aes(fraction,model))+geom_segment(aes(x=0,xend=fraction,yend=model),colour='#DFE6EB',linewidth=2.2)+geom_point(size=2,colour=blue)+geom_text(aes(label=sprintf('%s (%.1f%%)',comma(n),100*fraction)),nudge_y=.29,hjust=.5,size=1.9,colour=ink)+scale_x_continuous(labels=percent,limits=c(-.025,.34),breaks=c(0,.1,.2,.3),expand=expansion(mult=c(0,.02)))+labs(title='Complementary prediction yields',subtitle='18,454 expression-supported smORFs',x='Predicted fraction (threshold 0.5)',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=5.5))
inter <- read('exact_intersection_counts.tsv');stopifnot(sum(inter$n_sequences)==18454)
ints <- inter[exact_pattern!='None'];stopifnot(all(ints$n_sequences>0));ints[,support:=lengths(strsplit(exact_pattern,' + ',fixed=TRUE))];ints[,tier:=ifelse(support==1,'Single-model only',ifelse(grepl('ESM2-LGBM|Macrel',exact_pattern),'Conservative','Expanded-only'))];ints[,x:=.31+(seq_len(.N)-1)*.044]
stopifnot(ints[support>=2,sum(n_sequences)]==2028,ints[tier=='Conservative',sum(n_sequences)]==567)
write(ints,'Fig2e_intersections_and_tiers.tsv')
tx <- function(label,x,y,size=6,face='plain',col=ink,just='left')grid.text(label,x,y,just=just,gp=gpar(fontfamily=FONT,fontsize=size,fontface=face,col=col))
line <- function(x,y,xend,yend,col=ink,lwd=.5)grid.lines(c(x,xend),c(y,yend),gp=gpar(col=col,lwd=lwd))
draw_upset <- function(){
 tx('Four-model intersections define the catalogue',.047,.97,7,'bold');tx('2,028 candidates with >=2 predictors; 11,246 with no positive prediction',.015,.912,5.6)
 for(k in c(1,10,100,1000)){ yy<-.40+log10(k)/3.6*.40;line(.29,yy,.96,yy,'#E5E9EC',.4);tx(as.character(k),.275,yy,5,just='right') }
 tx('Count (log10)',.07,.66,5.5)
 for(i in seq_len(nrow(ints))){x<-ints$x[i];height<-log10(ints$n_sequences[i])/3.6*.40;grid.rect(x,.40+height/2,width=.031,height=max(.001,height),gp=gpar(fill=tier_cols[ints$tier[i]],col=NA));tx(as.character(ints$n_sequences[i]),x,.415+height,4.8,just='centre')
  active<-vapply(c('ESM2-LGBM','Macrel','AMPlify','AMP Scanner v2'),function(m)grepl(m,ints$exact_pattern[i],fixed=TRUE),logical(1));ys<-c(.29,.235,.18,.125)
  if(sum(active)>1)line(x,min(ys[active]),x,max(ys[active]),tier_cols[ints$tier[i]],1)
  grid.points(rep(x,4),ys,pch=16,size=unit(1.25,'mm'),gp=gpar(col=ifelse(active,tier_cols[ints$tier[i]],'#E1E6EA')))
 }
 for(k in 1:4)tx(c('ESM2-LGBM','Macrel','AMPlify','AMP Scanner v2')[k],.26,c(.29,.235,.18,.125)[k],5.4,just='right')
 for(k in 1:3){x<-c(.02,.34,.68)[k];grid.points(x,.036,pch=15,size=unit(1.8,'mm'),gp=gpar(col=tier_cols[k]));tx(names(tier_cols)[k],x+.025,.036,5.1)}
}
nov <- read('novelty_class_summary.tsv')
strict4<-read('candidate_evidence.tsv')[support_count_4==4,.N,by=novelty_class];strict4[,`:=`(set='All four',fraction=N/sum(N))]
nov<-rbind(nov,strict4,use.names=TRUE);write(nov,'Fig2f_novelty_with_four_model_subset.tsv')
nov[,class:=factor(c('Known-like'='Known-like','Remote homologous'='Distant match','Potentially novel'='Database-unmatched')[novelty_class],levels=names(nov_cols))];nov[,set:=factor(set,levels=c('Conservative','Expanded-only','All >=2/4','All four'))]
pF <- ggplot(nov,aes(set,fraction,fill=class))+geom_col(width=.65,colour='white',linewidth=.18)+geom_text(aes(label=ifelse(fraction>.04,sprintf('%s\n%.1f%%',comma(N),100*fraction),'')),position=position_stack(vjust=.5),size=1.9,colour='white',lineheight=.9)+scale_fill_manual(values=nov_cols)+scale_x_discrete(labels=c('Conserv.\nn = 567','Expanded\nn = 1,461','Complete\nn = 2,028','All four\nn = 29'))+scale_y_continuous(labels=percent,breaks=c(0,.5,1),expand=expansion(mult=c(0,.015)))+labs(title='Novelty across prediction-support tiers',subtitle='20 of 29 four-model candidates remain unmatched',x=NULL,y='Fraction of candidates')+theme(legend.position='bottom',legend.text=element_text(size=5),legend.key.width=unit(2,'mm'),axis.text.x=element_text(size=5.2))+guides(fill=guide_legend(nrow=2,byrow=TRUE))
best <- read('candidate_best_hits_by_database.tsv.gz');best[,has_hit:=!is.na(sseqid)&nzchar(sseqid)&aln_length>0]; nh<-best[,.(total=.N,aligned=sum(has_hit),no_hit=sum(!has_hit)),by=database];stopifnot(all(nh$total==2028));write(nh,'Fig2g_alignment_accounting.tsv')
best[,class:=factor(c('Known-like'='Known-like','Remote homologous'='Distant match','Potentially novel'='Database-unmatched')[database_class],levels=names(nov_cols))]
db_labels<-setNames(paste0(nh$database,'\n',comma(nh$aligned),' aligned; ',comma(nh$no_hit),' no hit'),nh$database)
best[,db_label:=unname(db_labels[database])]
pG <- ggplot(best[has_hit==TRUE],aes(bidirectional_coverage,pident,colour=class))+geom_vline(xintercept=50,linetype=2,linewidth=.25,colour='#9AA8B3')+geom_vline(xintercept=80,linewidth=.25,colour='#9AA8B3')+geom_hline(yintercept=40,linetype=2,linewidth=.25,colour='#9AA8B3')+geom_hline(yintercept=80,linewidth=.25,colour='#9AA8B3')+geom_point(size=.65,alpha=.45,stroke=0)+facet_wrap(~db_label,nrow=1)+scale_colour_manual(values=nov_cols)+scale_x_continuous(breaks=c(0,50,100),limits=c(0,102))+scale_y_continuous(breaks=c(0,50,100),limits=c(0,102))+labs(title='Sequence identity and alignment coverage',subtitle='One point per candidate with a reported alignment in each database',x='Minimum of query and subject coverage (%)',y='Amino-acid identity (%)',caption='No-hit candidates are counted above, not assigned artificial zero coordinates.')+theme(legend.position='none')
up <- read('SourceData_Fig3b.tsv');cons<-up[candidate_set=='Conservative'];allup<-up[candidate_set=='All expanded'];stopifnot(all(cons$metric==allup$metric));expup<-copy(cons);expup[,`:=`(candidate_set='Expanded-only',n=1461,observed_n=allup$observed_n-cons$observed_n)];up2<-rbind(cons,expup);up2[,fraction:=observed_n/n];z<-qnorm(.975);up2[,`:=`(ci_low=(fraction+z*z/(2*n)-z*sqrt(fraction*(1-fraction)/n+z*z/(4*n*n)))/(1+z*z/n),ci_high=(fraction+z*z/(2*n)+z*sqrt(fraction*(1-fraction)/n+z*z/(4*n*n)))/(1+z*z/n))];up2[,metric_label:=factor(rep(c('Upstream codon','5-prime complete','In-frame stop'),2),levels=rev(c('Upstream codon','5-prime complete','In-frame stop')))];write(up2,'Fig2h_disjoint_coding_context.tsv')
pH<-ggplot(up2,aes(fraction,metric_label,colour=candidate_set))+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',position=position_dodge(width=.5),width=.16,linewidth=.4)+geom_point(position=position_dodge(width=.5),size=1.7)+scale_colour_manual(values=tier_cols)+scale_x_continuous(labels=percent,limits=c(0,.19),breaks=c(0,.05,.1,.15))+labs(title='Recovered coding-context evidence',subtitle='Disjoint sets; points and 95% Wilson intervals',x='Fraction of candidates',y=NULL,caption='Stop observed: 56/567 vs 1/1,461.\nMissing upstream context is not evidence of absence.')+theme(legend.position='bottom',legend.text=element_text(size=5))
sen <- read('novelty_threshold_sensitivity.tsv');sens<-sen[known_identity==80 & class=='Potentially novel'];write(sens,'Fig2i_threshold_sensitivity.tsv')
pI<-ggplot(sens,aes(factor(remote_identity),factor(remote_coverage),fill=fraction))+geom_tile(colour='white',linewidth=1)+geom_text(aes(label=sprintf('%s\n%.1f%%',comma(N),100*fraction)),size=2.2,lineheight=.95,colour=ink)+scale_fill_gradient(low='#EAF3F3',high='#8CCAC2',limits=c(.65,1))+labs(title='Novelty classification across thresholds',subtitle='Database-unmatched candidates / 2,028',x='Minimum identity for a distant match (%)',y='Minimum coverage (%)')+theme(legend.position='none',axis.line=element_blank(),axis.ticks=element_blank())+annotate('rect',xmin=1.5,xmax=2.5,ymin=.5,ymax=1.5,fill=NA,colour=ink,linewidth=.6)
fam<-read('candidate_family_class_counts.tsv');fam[,class:=factor(family_class,levels=rev(family_class))];write(fam,'Fig2j_family_composition.tsv')
pJ<-ggplot(fam,aes(families,class))+geom_col(width=.52,fill=blue)+geom_text(aes(label=paste0(comma(families),' / ',comma(candidates))),hjust=0,nudge_x=20,size=1.9)+scale_x_continuous(limits=c(0,1880),breaks=c(0,500,1000),expand=expansion(mult=c(0,.01)))+scale_y_discrete(labels=c('Candidate-only family'='Unmatched-only','Database-expanded family'='Distant-match associated','Known-associated family'='Known-like associated'))+labs(title='Diversity across sequence families',subtitle='1,466 families; labels show families / peptides',x='Number of sequence families',y=NULL,caption='40% identity; 70% coverage; alignment >=15 aa')+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
box <- function(x,y,w,h,title,body,col){grid.roundrect(x,y,w,h,r=unit(1.7,'mm'),gp=gpar(fill=col,col=NA));tx(title,x,y+h*.22,6.2,'bold',just='centre');tx(body,x,y-h*.12,5.4,just='centre')}
arr <- function(x1,x2,y){grid.lines(c(x1,x2),c(y,y),arrow=arrow(length=unit(1.4,'mm'),type='closed'),gp=gpar(col='#788E9D',lwd=.65))}
draw_a<-function(){
 tx('Language-model representation and consensus discovery',.035,.94,7.7,'bold')
 box(.107,.55,.188,.50,'Labelled training peptides','45,731 total\n4,574 fine-tuning / 41,157 LGBM','#E8F1F5')
 box(.327,.55,.196,.50,'ESM2 650M','33 layers: 29 frozen + 4 tuned\nContextual peptide embeddings','#DBEAF0')
 box(.550,.55,.192,.50,'Sequence representation','CLS / mean / mean-max\nMean-max retained (2,560-D)','#EBE6F3')
 box(.767,.55,.190,.50,'10-fold LightGBM','Average ensemble probability\nDecision threshold = 0.5','#E7EFE4')
 arr(.205,.226,.55);arr(.429,.451,.55);arr(.65,.669,.55)
 grid.roundrect(.944,.55,.11,.50,r=unit(1.7,'mm'),gp=gpar(fill='#F4EADB',col=NA));tx('Candidate',.944,.64,5.8,'bold',just='centre');tx('AMP score',.944,.49,5.6,just='centre');arr(.866,.887,.55)
 tx('Held-out benchmark: 5,082 sequences',.01,.15,5.8,col=blue)
 tx('Discovery: 18,454 smORFs',.39,.15,5.8,col=ink)
 tx('4 predictors -> 2,028 candidates -> 567 conservative',.66,.15,5.6,col=blue)
}
put <- function(p,x,y,w,h,tag){pushViewport(viewport(x=unit(x,'mm'),y=unit(238-y-h,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','bottom')));if(is.function(p))p() else print(p,newpage=FALSE);upViewport();grid.text(tag,x=unit(x+.2,'mm'),y=unit(238-y+.1,'mm'),just=c('left','top'),gp=gpar(fontfamily=FONT,fontface='bold',fontsize=8.5,col=ink))}
draw_main<-function(){grid.newpage();put(draw_a,3,3,177,32,'a');put(pB,3,39,59,44,'b');put(pC,63,39,58,44,'c');put(pD,122,39,59,44,'d');put(draw_upset,3,88,108,49,'e');put(pF,113,88,68,49,'f');put(pG,3,141,103,49,'g');put(pH,108,141,73,49,'h');put(pI,3,196,86,39,'i');put(pJ,92,196,89,39,'j')}
# Supplement S1: all records contribute to kernel densities; no subsampling.
features<-c(length='Length (aa)',basic_residues_pct='Basic residues (%)',isoelectric_point='Isoelectric point',charge_ph7='Charge at pH 7',instability_index='Instability index',boman_index='Boman index',hydrophobic_moment='Hydrophobic moment',gravy='GRAVY')
groups<-c('Conservative','Complete','APD6','AMPSphere');paths<-c('conservative_consensus','expanded_consensus_ge2of4','apd6','ampsphere')
dens<-list();qa<-list();sumrows<-list()
for(i in seq_along(groups)){
 dat<-fread(file.path(O,'inputs/physchem',paste0('physchem_features_',paths[i],'.tsv.gz')),select=names(features))
 for(k in names(features)){v<-dat[[k]];ok<-is.finite(v);dd<-density(v[ok],n=512);dens[[length(dens)+1]]<-data.table(group=groups[i],feature=features[k],x=dd$x,y=dd$y);qa[[length(qa)+1]]<-data.table(group=groups[i],feature=k,input=length(v),finite=sum(ok),excluded=sum(!ok));sumrows[[length(sumrows)+1]]<-data.table(group=groups[i],feature=k,n=sum(ok),median=median(v[ok]),q25=quantile(v[ok],.25),q75=quantile(v[ok],.75)) }
 rm(dat);gc()
}
dd<-rbindlist(dens);dd[,feature:=factor(feature,levels=features)];dd[,group:=factor(group,levels=groups)];write(dd,'S1_density_coordinates.tsv');write(rbindlist(qa),'QA_physchem_record_accounting.tsv');write(rbindlist(sumrows),'Table_S1_physchem_summary.tsv')
pS1<-ggplot(dd,aes(x,y,colour=group))+geom_line(linewidth=.52)+facet_wrap(~feature,scales='free',ncol=4)+scale_colour_manual(values=c(blue,purple,gold,teal),labels=c('Conservative (567)','Complete (2,028)','APD6 (3,306)','AMPSphere (863,498)'))+labs(title='Physicochemical profiles of the candidate and reference collections',subtitle='All available valid records; complete catalogue includes the conservative subset',x=NULL,y='Density')+theme(legend.position='bottom',legend.text=element_text(size=6),panel.spacing=unit(4,'mm'))
# Supplement S2: compact model metrics, pooling uncertainty, length strata and decision trade-off.
ml<-melt(metrics,id.vars=c('Model','short'),measure.vars=c('AUROC','AUPRC','Precision','Recall','F1','MCC'),variable.name='metric');ml[,short:=factor(short,levels=rev(names(model_cols)))]
q1<-ggplot(ml,aes(metric,short,fill=value))+geom_tile(colour='white',linewidth=.5)+geom_text(aes(label=sprintf('%.3f',value)),size=2)+scale_fill_gradient(low='#EEF2F5',high='#9FC7D9',limits=c(.45,1))+labs(title='a  Seven-model benchmark',x=NULL,y=NULL)+theme(legend.position='none',axis.line=element_blank(),axis.ticks=element_blank())
ci<-read('fig2_bootstrap_ci.tsv');q2<-ggplot(ci,aes(estimate,model,colour=metric))+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',width=.2,position=position_dodge(width=.4),linewidth=.4)+geom_point(position=position_dodge(width=.4),size=1.6)+scale_colour_manual(values=c(blue,gold))+labs(title='b  Pooling performance and uncertainty',subtitle='2,000 bootstrap replicates; 95% intervals',x='Metric',y=NULL)+theme(legend.position='bottom')
le<-read('fig2_length_stratified_metrics.tsv');le[,length_group:=factor(length_group,levels=c('<=20','21-30','31-50','51-75','>75'))];q3<-ggplot(le,aes(length_group,AUROC,colour=model,group=model))+geom_line(linewidth=.4)+geom_point(size=1.5)+scale_colour_manual(values=c(purple,gold,blue))+scale_y_continuous(limits=c(.90,1),breaks=c(.90,.95,1))+labs(title='c  Length-stratified discrimination',subtitle='n = 1,353 / 755 / 2,275 / 374 / 325 sequences',x='Peptide length (aa)',y='AUROC')+theme(legend.position='bottom')
tr<-read('fig2_meanmax_threshold_tradeoff.tsv');cat('THRESHOLD COLS',names(tr),'\n')
nm<-names(tr);th<-nm[tolower(nm)=='threshold'];keep<-nm[tolower(nm)%in%c('precision','recall','f1')];tl<-melt(tr,id.vars=th,measure.vars=keep);setnames(tl,th,'threshold')
q4<-ggplot(tl,aes(threshold,value,colour=variable))+geom_line(linewidth=.5)+geom_vline(xintercept=.5,linetype=2,linewidth=.3)+scale_colour_manual(values=c(blue,gold,purple))+labs(title='d  ESM2-LGBM threshold trade-off',x='Probability threshold',y='Metric')+theme(legend.position='bottom')
# Supplement S3: new cross-tabulation of existing support and novelty annotations.
ev<-read('candidate_evidence.tsv');ev[,tier:=ifelse(support_count_4==2,ifelse(is_conservative,'2/4 + ESM2/Macrel','2/4 AMPlify + Scanner'),paste0(support_count_4,'/4'))]
ct<-ev[,.N,by=.(tier,novelty_class)];ct[,fraction:=N/sum(N),by=tier];write(ct,'S3_support_by_novelty.tsv')
ct[,class:=c('Known-like'='Known-like','Remote homologous'='Distant match','Potentially novel'='Database-unmatched')[novelty_class]]
r1<-ggplot(ct,aes(tier,fraction,fill=class))+geom_col(width=.65)+geom_text(aes(label=ifelse(fraction>.04,paste0(N,'\n',sprintf('%.1f%%',fraction*100)),'')),position=position_stack(vjust=.5),size=2.1,colour='white')+scale_fill_manual(values=nov_cols)+scale_y_continuous(labels=percent)+labs(title='a  Sequence novelty by prediction support',x=NULL,y='Fraction within support stratum')+theme(legend.position='bottom',axis.text.x=element_text(angle=15,hjust=1))
fs<-read('candidate_family_summary.tsv');fs[,size_bin:=factor(ifelse(n_candidates>=5,'5+',as.character(n_candidates)),levels=c('1','2','3','4','5+'))];fc<-fs[,.N,by=.(size_bin,family_class)];write(fc,'S3_family_size_distribution.tsv')
r2<-ggplot(fc,aes(size_bin,N,fill=family_class))+geom_col(width=.65,position='stack')+scale_fill_manual(values=c(teal,gold,ink))+labs(title='b  Family size and database association',subtitle='40% identity; 70% bidirectional coverage; >=15-aa alignment',x='Candidates per sequence family',y='Number of families')+theme(legend.position='bottom',legend.text=element_text(size=5))+guides(fill=guide_legend(nrow=3))
us<-read('figure3B_upstream_by_prediction_stratum.tsv');r3<-ggplot(us,aes(prediction_stratum,fraction,colour=metric))+geom_errorbar(aes(ymin=ci_low,ymax=ci_high),position=position_dodge(width=.5),width=.18,linewidth=.4)+geom_point(position=position_dodge(width=.5),size=1.6)+scale_colour_manual(values=c(blue,gold,purple))+scale_y_continuous(labels=percent)+labs(title='c  Coding evidence by prediction support',x=NULL,y='Fraction; 95% Wilson interval')+theme(legend.position='bottom',axis.text.x=element_text(angle=15,hjust=1))
write(ct[tier=='2/4 + ESM2/Macrel' | tier=='3/4' | tier=='4/4'],'QA_conservative_strata.tsv')
# Supplement S4 retains the original 10,000-ORF analysis.
P4<-file.path(O,'inputs/orf10000')
orfs<-fread(file.path(P4,'smorf_pred_lgbm_esm2_mm.csv'));ma<-fread(file.path(P4,'macrel_out/macrel.out.prediction.gz'),skip='Access')
stopifnot(nrow(orfs)==10000);orfs[,id:=paste0('smorf_',seq_len(.N)-1L)];esid<-orfs[esm2_lgbm_pred==1,id];maid<-ma[AMP_probability>=.5,Access]
stopifnot(length(esid)==635,length(maid)==928,length(intersect(esid,maid))==98)
orfs[,set:=ifelse(id%in%esid & id%in%maid,'Shared',ifelse(id%in%esid,'ESM2-LGBM only',ifelse(id%in%maid,'Macrel only','Neither')))];write(orfs[,.(id,Length,esm2_lgbm_prob,esm2_lgbm_pred,set)],'S4_ORF10000_predictions.tsv')
hh<-hist(orfs$esm2_lgbm_prob,breaks=seq(0,1,.025),plot=FALSE);hd<-data.table(xmin=head(hh$breaks,-1),xmax=tail(hh$breaks,-1),n=hh$counts);write(hd,'S4_score_histogram.tsv')
s4a<-ggplot(hd[n>0])+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=.8,ymax=n),fill=blue,colour='white',linewidth=.15)+scale_y_log10(breaks=c(1,10,100,1000,10000),labels=comma)+geom_vline(xintercept=.5,colour=gold,linetype=2,linewidth=.45)+labs(title='a  Prediction scores in 10,000 sampled ORFs',subtitle='635 candidates at a probability threshold of 0.5',x='ESM2-LGBM probability',y='ORF count (log10)')
ld<-rbindlist(lapply(c('ESM2-LGBM only','Macrel only','Shared'),function(g){v<-orfs[set==g,Length];stopifnot(all(is.finite(v)&v>0));a<-density(log10(v),from=log10(8),to=log10(2500),n=512);data.table(set=paste0(g,' (n = ',length(v),')'),x=a$x,y=a$y)}));write(ld,'S4_length_density.tsv')
s4b<-ggplot(ld,aes(x,y,colour=set))+geom_line(linewidth=.6)+scale_colour_manual(values=c(blue,purple,teal))+scale_x_continuous(breaks=log10(c(10,20,50,200,500,2000)),labels=c('10','20','50','200','500','2,000'))+labs(title='b  Complementary peptide-length coverage',x='Predicted peptide length (aa; log10)',y='Density of log10 length')+theme(legend.position='bottom',legend.text=element_text(size=5.8))+guides(colour=guide_legend(nrow=2))
capture.output(sessionInfo(),file=file.path(O,'work/R_sessionInfo.txt'))
cat('ALL_FIGURES_COMPLETE\n')