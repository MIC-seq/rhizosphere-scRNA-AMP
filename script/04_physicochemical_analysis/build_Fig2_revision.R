suppressPackageStartupMessages({library(ggplot2);library(grid);library(patchwork);library(data.table)})
O <- '/public3/labmember/qianqh/soil_amp_analysis/13_model_test'
source(file.path(O,'scripts/prepare_panels.R'))
clean <- function(p) p+labs(title=NULL,subtitle=NULL,caption=NULL,tag=NULL)+theme(plot.margin=margin(4,4,4,4))
pb<-clean(pB);pb$layers<-pb$layers[1:3];pc<-clean(pC);pc$layers<-pc$layers[-length(pc$layers)]
zin<-function(which){ggplot(curves,aes(if(which=='roc')fpr else tpr,if(which=='roc')tpr else precision,colour=model))+geom_path(linewidth=.35)+scale_colour_manual(values=model_cols)+coord_cartesian(xlim=if(which=='roc')c(.02,.15)else c(.75,.92),ylim=if(which=='roc')c(.91,1)else c(.94,1))+labs(x=NULL,y=NULL)+theme_classic(base_size=4.8,base_family='Arial')+theme(legend.position='none',plot.background=element_rect(fill='white',colour='#A6B3BB',linewidth=.25),axis.line=element_line(linewidth=.2),axis.ticks=element_line(linewidth=.2),axis.ticks.length=unit(.5,'mm'),plot.margin=margin(2,2,2,2))}
pb<-pb+inset_element(zin('roc'),left=.58,bottom=.10,right=.98,top=.40)
pc<-pc+inset_element(zin('pr'),left=.10,bottom=.08,right=.52,top=.40)
pl<-clean(s4b)+theme(legend.position='bottom',legend.text=element_text(size=5),legend.key.width=unit(3,'mm'))+guides(colour=guide_legend(nrow=3))
pu<-ggplot(up,aes(metric,fraction,fill=candidate_set))+geom_col(position=position_dodge(.7),width=.65)+geom_errorbar(aes(ymin=ci_low,ymax=ci_high),position=position_dodge(.7),width=.15,linewidth=.4)+geom_text(aes(y=ci_high+.008,label=observed_n),position=position_dodge(.7),size=2)+scale_fill_manual(values=c('Conservative'=blue,'All expanded'=teal),breaks=c('Conservative','All expanded'),labels=c('Conservative (567)','Complete (2,028)'))+scale_y_continuous(labels=percent,limits=c(0,.215),breaks=c(0,.05,.1,.15,.2),expand=expansion(mult=c(0,0)))+scale_x_discrete(labels=c('≥1 upstream codon'='Upstream\ncodon','5′ complete'="5'-complete\nORF",'Upstream stop observed'='In-frame\nstop'))+labs(x=NULL,y='Fraction of candidates')+theme(legend.position='bottom',legend.text=element_text(size=5))+guides(fill=guide_legend(nrow=2))
pn<-clean(pF %+% nov[set!='All four'])+scale_x_discrete(drop=TRUE,labels=c('Conservative\nn = 567','Expanded-only\nn = 1,461','Complete\nn = 2,028'))+theme(axis.text.x=element_text(size=5.5))+guides(fill=guide_legend(nrow=2))
pg<-clean(pG);pg$data$db_label<-sub('\n.*','',pg$data$db_label)
pg<-pg+labs(x='Bidirectional coverage (%)')
pn<-pn+scale_fill_manual(values=nov_cols,labels=c('Potentially novel','Distant match','Known-like'))
# Preserve original workflow structure, with explicit allocations and pooling dimensions.
roundbox<-function(x,y,w,h,fill,border=NA){grid.roundrect(x,y,w,h,r=unit(1.2,'mm'),gp=gpar(fill=fill,col=border,lwd=.6))}
draw_workflow<-function(){
 bands<-list(c(.23,.43,'Sequence encoding','#D9EAF0'),c(.60,.26,'Sequence pooling','#E7E0F1'),c(.82,.15,'Classification','#E2ECD7'),c(.945,.09,'Output','#F1DFD5'))
 for(v in bands){roundbox(as.numeric(v[1]),.945,as.numeric(v[2]),.10,v[4]);tx(v[3],as.numeric(v[1]),.945,6.2,'bold',just='centre')}
 roundbox(.063,.54,.118,.38,'#EDF4F6',blue);tx('Training\npeptides',.063,.59,6.1,'bold',just='centre');tx('n = 45,731',.063,.42,5.4,just='centre')
 grid.points(seq(.025,.098,length.out=6),rep(.70,6),pch=16,size=unit(1,'mm'),gp=gpar(col=blue))
 roundbox(.179,.54,.086,.23,'#F1F6F8',blue);tx('Tokenizer',.179,.54,5.7,'bold',just='centre');arr(.123,.132,.54)
 roundbox(.319,.54,.154,.58,'#E5EFF3',blue);tx('ESM2 (650M)',.319,.75,6.4,'bold',just='centre');tx('33 transformer layers',.319,.66,5.1,just='centre')
 roundbox(.283,.45,.064,.25,'#CFE0E8');roundbox(.354,.45,.064,.25,'#B5D2DF');tx('1-29',.283,.49,5.5,'bold',just='centre');tx('frozen',.283,.40,5,just='centre');tx('30-33',.354,.49,5.5,'bold',just='centre');tx('fine-tuned',.354,.40,4.8,just='centre');tx('Fine-tuning n = 4,574',.319,.284,4.8,just='centre');arr(.225,.24,.54)
 roundbox(.45,.54,.090,.27,'#F0ECF7',purple);tx('Embeddings',.45,.59,5.7,'bold',just='centre');tx('L x 1,280',.45,.48,5.2,just='centre');arr(.398,.403,.54)
 for(i in 1:3){y<-c(.75,.48,.20)[i];roundbox(.59,y,.137,.22,if(i==3)'#DED5EB'else'#F5F1FA',purple);tx(c('CLS','Mean','Mean + max')[i],.59,y+.037,5.8,'bold',just='centre');tx(c('1,280-D','1,280-D','2,560-D')[i],.59,y-.044,5.1,just='centre');line(.501,y,.518,y,purple,.6)}
 line(.501,.20,.501,.75,purple,.6);arr(.495,.501,.54)
 roundbox(.799,.49,.190,.66,'#EFF4E9','#8AAB74');tx('10-fold LightGBM',.799,.75,6.2,'bold',just='centre');tx('n = 41,157',.799,.67,5.1,just='centre')
 for(i in 1:4){x<-c(.76,.84,.76,.84)[i];y<-c(.53,.53,.34,.34)[i];roundbox(x,y,.062,.13,'#DDE8CB','#9EB783');tx(c('Fold 1','Fold 2','Fold 9','Fold 10')[i],x,y,4.8,just='centre')};tx('...',.8,.43,5.5,just='centre');tx('Mean probability',.799,.20,5.2,just='centre');line(.66,.20,.69,.20,purple,.6);line(.69,.20,.69,.49,purple,.6);arr(.69,.702,.49)
 roundbox(.949,.62,.09,.25,'#F7EEE9','#C28D79');tx('P >= 0.5',.949,.66,5.7,'bold',just='centre');tx('threshold',.949,.57,4.9,just='centre');arr(.897,.902,.62)
 roundbox(.949,.24,.09,.23,'#F1DFD5','#C28D79');tx('Candidate',.949,.28,5.3,'bold',just='centre');tx('AMP',.949,.20,5.5,'bold',just='centre');grid.lines(c(.949,.949),c(.49,.375),arrow=arrow(length=unit(1.2,'mm'),type='closed'),gp=gpar(col=ink,lwd=.6))
}
draw_upset_min<-function(){
 # Quantitative geometry and all 15 intersections are retained; prose moves to caption.
 for(k in c(1,10,100,1000)){yy<-.40+log10(k)/3.6*.45;line(.285,yy,.97,yy,'#E8EDF0',.4);tx(as.character(k),.268,yy,5.2,just='right')}
 tx('Intersection size (log10)',.02,.91,5.7)
 for(i in seq_len(nrow(ints))){x<-ints$x[i];height<-log10(ints$n_sequences[i])/3.6*.45;grid.rect(x,.40+height/2,.031,max(.001,height),gp=gpar(fill=tier_cols[ints$tier[i]],col=NA));tx(as.character(ints$n_sequences[i]),x,.416+height,5.1,just='centre');act<-vapply(c('ESM2-LGBM','Macrel','AMPlify','AMP Scanner v2'),function(m)grepl(m,ints$exact_pattern[i],fixed=TRUE),logical(1));ys<-c(.29,.235,.18,.125);if(sum(act)>1)line(x,min(ys[act]),x,max(ys[act]),tier_cols[ints$tier[i]],1);grid.points(rep(x,4),ys,pch=16,size=unit(1.35,'mm'),gp=gpar(col=ifelse(act,tier_cols[ints$tier[i]],'#DFE5E9')))}
 for(k in 1:4)tx(c('ESM2-LGBM','Macrel','AMPlify','AMP Scanner v2')[k],.26,c(.29,.235,.18,.125)[k],5.5,just='right')
 for(k in 1:3){x<-c(.03,.35,.68)[k];grid.points(x,.03,pch=15,size=unit(1.8,'mm'),gp=gpar(col=tier_cols[k]));tx(names(tier_cols)[k],x+.022,.03,5.1)}
}
place<-function(p,x,y,w,h,tag,H=224){pushViewport(viewport(x=unit(x,'mm'),y=unit(H-y-h,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','bottom')));if(is.function(p))p()else print(p,newpage=FALSE);upViewport();grid.text(tag,x=unit(x,'mm'),y=unit(H-y+1,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=9,fontface='bold',col=ink))}
main<-function(){grid.newpage();place(draw_workflow,4,5,175,47,'a');place(pb,4,58,58,47,'b');place(pc,63,58,57,47,'c');place(clean(s4a),123,58,56,47,'d');
 for(i in seq_along(model_cols)){xx<-9+((i-1)%%4)*29;yy<-108+floor((i-1)/4)*3;grid.lines(unit(c(xx,xx+3),'mm'),unit(rep(224-yy,2),'mm'),gp=gpar(col=model_cols[i],lwd=1.3));grid.text(names(model_cols)[i],unit(xx+4,'mm'),unit(224-yy,'mm'),just='left',gp=gpar(fontfamily='Arial',fontsize=5,col=ink))}
 place(pl,4,118,58,48,'e');place(draw_upset_min,65,118,114,48,'f');place(pu,4,174,48,47,'g');place(pg,55,174,65,47,'h');place(pn+theme(axis.text.x=element_text(size=5))+guides(fill=guide_legend(nrow=3)),123,174,56,47,'i')}
savefig(main,'Fig2_main',183,224)
panels<-list(a=draw_workflow,b=pb,c=pc,d=clean(s4a),e=pl,f=draw_upset_min,g=pu,h=pg,i=pn)
for(n in names(panels)){p<-panels[[n]];if(is.function(p)){fun<-local({pp<-p;function(){grid.newpage();pp()}});savefig(fun,paste0('Fig2',n),if(n=='a')183 else 120,if(n=='a')50 else 60,'panels',FALSE)}else savefig(p,paste0('Fig2',n),if(n=='h')110 else 89,65,'panels',FALSE)}
# Keep original Fig. 3a as eight property panels without headline text.
savefig(clean(pS1),'FigS3_physicochemical',183,113)
pool<-read('fig2_test_metrics.tsv');pol<-melt(pool,id.vars='model',measure.vars=c('AUROC','AUPRC','F1','MCC'),variable.name='metric');pp<-ggplot(pol,aes(metric,value,colour=model,group=model))+geom_line(linewidth=.45)+geom_point(size=1.6)+scale_colour_manual(values=c('CLS pooling'=purple,'Mean pooling'=gold,'Mean-max pooling'=blue))+coord_cartesian(ylim=c(.89,1))+labs(x=NULL,y='Metric value')+theme(legend.position='bottom',legend.text=element_text(size=5.5))
td<-melt(tr,id.vars='threshold',measure.vars=c('Precision','Recall','Specificity','F1'));pt<-ggplot(td,aes(threshold,value,colour=variable))+geom_line(linewidth=.5)+geom_vline(xintercept=.5,linetype=2,linewidth=.3)+scale_colour_manual(values=c(blue,purple,teal,gold))+labs(x='Decision threshold',y='Metric value')+theme(legend.position='bottom')
sup<-data.table(support=c('4/4','3/4','2/4','1/4'),n=c(29,211,1788,5180));sup[,support:=factor(support,levels=support)];ps<-ggplot(sup,aes(support,n))+geom_col(fill=blue,width=.65)+geom_text(aes(label=comma(n)),vjust=-.4,size=2.3)+scale_y_continuous(expand=expansion(mult=c(0,.12)),labels=comma)+labs(x='Number of supporting models',y='Number of sequences')
# Prior S2e/f are retained as complete tables; the score histogram returns to main Fig2d.
# Three radar views: two comparator selections present in source results, plus a pooling view.
radar<-function(dat,ids,palette){
 mm<-c('AUROC','AUPRC','Accuracy','Precision','Recall','Specificity','F1');z<-melt(dat,id.vars='name',measure.vars=mm,variable.name='metric');z[,a:=2*pi*(match(metric,mm)-1)/7];z<-z[name%in%ids];z[,name:=factor(name,levels=ids)];setorder(z,name,a);z[,`:=`(x=sin(a)*value,y=cos(a)*value)];z<-rbind(z,z[a==0][,a:=2*pi]);setorder(z,name,a)
 ang<-2*pi*(0:7)/7;gd<-rbindlist(lapply(c(.25,.5,.75,1),function(v)data.table(x=sin(ang)*v,y=cos(ang)*v,v=v)));axis<-data.table(x=sin(ang[1:7]),y=cos(ang[1:7]),label=mm)
 ggplot(z,aes(x,y,colour=name,group=name))+geom_path(data=gd,aes(x,y,group=v),inherit.aes=FALSE,colour='#DCE3E8',linewidth=.3)+geom_segment(data=axis,aes(x=0,y=0,xend=x,yend=y),inherit.aes=FALSE,colour='#E5E9ED',linewidth=.25)+geom_path(linewidth=.6)+geom_text(data=axis,aes(x=x*1.25,y=y*1.18,label=label),inherit.aes=FALSE,size=1.95,colour=ink)+annotate('text',x=.08,y=c(.25,.5,.75,1),label=c('.25','.50','.75','1.00'),size=1.5,colour='#778895',hjust=0)+scale_colour_manual(values=palette)+coord_fixed(xlim=c(-1.65,1.65),ylim=c(-1.35,1.35),clip='off')+theme_void(base_family='Arial',base_size=6)+theme(legend.title=element_blank(),legend.position='bottom',legend.text=element_text(size=5.8),legend.key.width=unit(3,'mm'),plot.margin=margin(4,3,4,3))+guides(colour=guide_legend(ncol=1))
}
mr<-copy(metrics);mr[,name:=short];pr<-copy(pool);pr[,name:=model]
ra<-radar(mr,names(model_cols),model_cols);rb<-radar(mr,c('ESM2-LGBM','LightGBM','LSTM','Macrel'),model_cols);rc<-radar(pr,c('CLS pooling','Mean pooling','Mean-max pooling'),c('CLS pooling'=purple,'Mean pooling'=gold,'Mean-max pooling'=blue))
radar_page<-function(){
 grid.newpage();plots<-list(ra,rb,rc);ids<-list(names(model_cols),c('ESM2-LGBM','LightGBM','LSTM','Macrel'),c('CLS pooling','Mean pooling','Mean-max pooling'));pal<-list(model_cols,model_cols,c('CLS pooling'=purple,'Mean pooling'=gold,'Mean-max pooling'=blue))
 for(j in 1:3){x<-2+(j-1)*60;place(plots[[j]]+theme(legend.position='none'),x,4,59,49,letters[j],H=80)
  for(k in seq_along(ids[[j]])){yy<-56+(k-1)*3;grid.lines(unit(c(x+13,x+16),'mm'),unit(rep(80-yy,2),'mm'),gp=gpar(col=pal[[j]][ids[[j]][k]],lwd=1.5));grid.text(ids[[j]][k],unit(x+18,'mm'),unit(80-yy,'mm'),just='left',gp=gpar(fontfamily='Arial',fontsize=5.7,col=ink))}
 }
}
supp2<-function(){
 grid.newpage();place(pp,4,5,86,65,'a',H=155);place(pt,94,5,85,65,'b',H=155);place(ra+theme(legend.position='none'),4,78,86,58,'c',H=155);place(ps,94,78,85,72,'d',H=155)
 for(i in seq_along(model_cols)){xx<-12+floor((i-1)/4)*39;yy<-138+((i-1)%%4)*3.2;grid.lines(unit(c(xx,xx+3),'mm'),unit(rep(155-yy,2),'mm'),gp=gpar(col=model_cols[i],lwd=1.4));grid.text(names(model_cols)[i],unit(xx+4,'mm'),unit(155-yy,'mm'),just='left',gp=gpar(fontfamily='Arial',fontsize=5.7,col=ink))}
}
savefig(supp2,'FigS2_model_support',183,155)
savefig(ra+theme(legend.key.height=unit(2.4,'mm'),legend.spacing.y=unit(0,'mm')),'FigS2c_seven_model_radar',89,89,'panels',FALSE)
length_summary<-orfs[set!='Neither',.(n=.N,mean_length=mean(Length),median_length=median(Length)),by=set];write(length_summary,'Supplement_length_summary_recomputed.tsv');write(sup,'Supplement_support_counts.tsv')
capture.output(sessionInfo(),file=file.path(O,'logs/sessionInfo.txt'))
cat('REVISION_COMPLETE\n')
