.libPaths(c('/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/software/Rlib',.libPaths()))
suppressPackageStartupMessages({library(grid);library(data.table);library(svglite);library(jsonlite)})
b<-'/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/fg_MES_20261027';o<-file.path(b,'main_supp_v8');dir.create(o,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(FONTCONFIG_FILE='/etc/fonts/fonts.conf',FONTCONFIG_PATH='/etc/fonts')
font<-'Liberation Sans';ink<-'#303C46'
D<-fread(file.path(b,'tables/f_all45_author_MES.tsv'));rows<-fread(file.path(b,'tables/f_fixed_24_pathways.tsv'))
F<-fread(file.path(b,'source_data/replacement_g_soil_pairs.tsv'));Z<-fread(file.path(b,'source_data/replacement_g_cells.tsv'))
sp<-c('Chloracidobacterium validum','Halomonas sp. JS92-SW72','Pseudomonas aeruginosa');sn<-c('C. validum','Halomonas','P. aeruginosa');cl<-c('Nucleotide metabolism','Amino acid metabolism','Energy metabolism');crops<-c('soybean','rice','wheat')
sc<-c(S1='#567F9B',S2='#D27D63');spcol<-c('#6696B1','#75AA97','#BE829B');gcol<-c('Amino acid metabolism'='#91D6C3','Carbohydrate metabolism'='#EFE59B','Energy metabolism'='#B5ADD7','Lipid metabolism'='#EE978D','Nucleotide metabolism'='#86B4D0');bg<-c('#F1F6F2','#F1F5F9','#FAF6EC');hd<-c('#DFEADF','#DFE8F0','#F0E5CA')
tx<-function(s,x,y,size=5.2,face='plain',just='centre',col=ink,rot=0)grid.text(s,unit(x,'mm'),unit(y,'mm'),just=just,rot=rot,gp=gpar(fontfamily=font,fontsize=size,fontface=face,col=col,lineheight=.88))
rect<-function(x,y,w,h,fill,col=NA)grid.rect(unit(x,'mm'),unit(y,'mm'),unit(w,'mm'),unit(h,'mm'),gp=gpar(fill=fill,col=col,lwd=.4))
seg<-function(x,y,x2,y2,col='#7A858F',lwd=.5)grid.segments(unit(x,'mm'),unit(y,'mm'),unit(x2,'mm'),unit(y2,'mm'),gp=gpar(col=col,lwd=lwd))
pt<-function(x,y,col,size=.7,border=NA)grid.points(unit(x,'mm'),unit(y,'mm'),pch=21,size=unit(size,'mm'),gp=gpar(fill=col,col=border,lwd=.4))
pal<-colorRampPalette(c('#42779A','#F8F7F3','#CF6866'))(401);hot<-colorRampPalette(c('#FFF0AC','#EE9B5F','#B3283B'))(256)
heat<-function(){
 tx('f',0,107,8,'bold','left');tx('Pathway enrichment',5,107,6.2,'bold','left')
 xx<-c(39,45,52.5,58.5,64.5,72,78);ss<-rep(sp,c(2,3,2));uu<-c('S1','S2','S1','S2','S3','S1','S2');tw<-5.7;rh<-3.55;top<-97
 for(i in 1:3){jj<-which(ss==sp[i]);tx(sn[i],mean(xx[jj]),102,5.2,'italic');rect(mean(xx[jj]),99.5,diff(range(xx[jj]))+tw,.6,spcol[i])}
 lim<-ceiling(max(abs(D$MES))*10)/10
 for(i in 1:nrow(rows)){y<-top-(i-.5)*rh;lab<-rows$pathway[i];short<-c('Alanine, aspartate and glutamate metabolism'='Ala/Asp/Glu metabolism','Glycine, serine and threonine metabolism'='Gly/Ser/Thr metabolism','Cysteine and methionine metabolism'='Cys/Met metabolism','Valine, leucine and isoleucine degradation'='Val/Leu/Ile degradation','Amino sugar and nucleotide sugar metabolism'='Amino sugar / nucleotide sugar metab.','Glyoxylate and dicarboxylate metabolism'='Glyoxylate / dicarboxylate metab.');if(lab%in%names(short))lab<-short[[lab]];tx(lab,35,y,5,just='right');rect(82.3,y,1.5,rh-.10,gcol[[rows$group[i]]]);for(j in 1:7){q<-D[species==ss[j]&unit==uu[j]&pathway_id==rows$pathway_id[i]];col<-if(q$pathway_UMI>0)pal[max(1,min(401,round((q$MES+lim)/(2*lim)*400)+1))]else'#E2E5E7';rect(xx[j],y,tw,rh-.13,col)}}
 for(j in 1:7)tx(uu[j],xx[j],10,5)
 # Compact, centered MES bar under the heatmap rather than extending left.
 tx('MES',38.5,6.7,5.2,'bold','right');for(k in 1:100)rect(40+(k-.5)*.20,6.7,.21,1.7,pal[1+floor(400*(k-.5)/100)])
 for(v in c(-lim,0,lim))tx(format(v),40+20*(v+lim)/(2*lim),3.5,5)
 tx('Grey: 0 detected UMI',82.8,6.7,5,just='right')
 labs<-c('Amino acid','Carbohydrate','Energy','Lipid','Nucleotide');lx<-c(0,20,44,59,70)
 for(i in 1:5){rect(lx[i]+.8,.4,1.4,1.4,gcol[i]);tx(labs[i],lx[i]+2,.4,5,just='left')}
}
path<-function(){
 tx('g',0,107,8,'bold','left');tx('Single-cell ssGSEA',5,107,5.6,'bold','left');tx('Paired state MES',79,107,5.8,'bold')
 H<-25;XL<-4;XR<-c(48,69,90);W<-19;Y<-c(70,39,8)
 for(j in 1:3){rect(XR[j]+W/2,54,W,96,bg[j]);rect(XR[j]+W/2,100,W,5,hd[j]);tx(c('Soybean','Rice','Wheat')[j],XR[j]+W/2,100,5.5,'bold')}
 lims<-list(c(3.4,3.55),c(.35,.55),c(.5,1.75));ticks<-list(seq(3.4,3.55,.05),seq(.35,.55,.05),seq(.5,1.75,.25))
 for(i in 1:3){z<-copy(Z[species==sp[i]]);q<-copy(F[species==sp[i]]);yb<-Y[i]
 tx(paste0(sn[i],' | ',cl[i]),XL,yb+29,5.3,'italic','left')
 rx<-range(z$UMAP_1);ry<-range(z$UMAP_2);span<-max(diff(rx),diff(ry))*1.12;xx<-function(v)XL+H*((v-mean(rx))/span+.5);yy<-function(v)yb+H*((v-mean(ry))/span+.5)
 pt(xx(z$UMAP_1),yy(z$UMAP_2),'#D3D8DC',.68);obs<-z[pathway_UMI>0][order(raw_ssgsea)];rr<-range(obs$raw_ssgsea);cols<-hot[pmax(1,pmin(256,1+floor(255*(obs$raw_ssgsea-rr[1])/diff(rr))))];pt(xx(obs$UMAP_1),yy(obs$UMAP_2),cols,.68)
 cen<-z[,.(u=median(UMAP_1),v=median(UMAP_2)),by=state_code];for(k in 1:nrow(cen)){x<-xx(cen$u[k]);y<-yy(cen$v[k]);rect(x,y,3.7,2.7,'#FFFFFFE8','#87929A');tx(cen$state_code[k],x,y,5)}
 seg(XL,yb,XL+H,yb);seg(XL,yb,XL,yb+H)
 for(k in 1:100)rect(XL+7+(k-.5)*.13,yb+26.5,.14,1.4,hot[1+floor(255*(k-.5)/100)])
 tx(sprintf('%.2f',rr[1]),XL+5.8,yb+26.5,5,just='right');tx(sprintf('%.2f',rr[2]),XL+21.2,yb+26.5,5,just='left')
 lim<-lims[[i]];lim<-lim+c(-1,1)*diff(lim)*.035;yt<-function(v)yb+H*(v-lim[1])/diff(lim)
 for(j in 1:3){x0<-XR[j];a<-q[crop==crops[j]];setorder(a,sample);n<-nrow(a);rect(x0+W/2,yb+26.5,W,3,hd[j]);tx(paste0('n = ',n),x0+W/2,yb+26.5,5)
 for(v in ticks[[i]])seg(x0,yt(v),x0+W,yt(v),'#DCE3E6',.4);seg(x0,yb,x0+W,yb)
 if(j==1){seg(x0,yb,x0,yb+H);for(v in ticks[[i]]){seg(x0-.7,yt(v),x0,yt(v));tx(sprintf('%.2f',v),x0-1.3,yt(v),5,just='right')}}
 sx<-c(x0+4,x0+15);for(k in 1:2){seg(sx[k],yb,sx[k],yb+H,'#E2E7E8',.35);tx(paste0('S',k),sx[k],yb-2,5)}
 if(n){off<-if(n==1)0 else seq(-.8,.8,length.out=n);for(k in 1:n){xc<-sx+off[k];yc<-yt(c(a$MES_S1[k],a$MES_S2[k]));seg(xc[1],yc[1],xc[2],yc[2],'#9AA9B5',.75);pt(xc,yc,unname(sc),1.3,'white')}}else tx('No paired data',x0+W/2,yb+H/2,5,col='#89949C')
 }
 }
 tx('UMAP 2',.7,51,5,rot=90);tx('UMAP 1',16.5,3.8,5);tx('MES',36,51,5.2,rot=90);tx('COG state',79,3.8,5)
 pt(48,.5,'#D3D8DC',1);tx('0 UMI',50,.5,5,just='left');pt(71,.5,sc[1],1);tx('S1',73,.5,5,just='left');pt(91,.5,sc[2],1);tx('S2',93,.5,5,just='left')
}
draw<-function(){grid.newpage();pushViewport(viewport(x=unit(0,'mm'),y=unit(2,'mm'),width=unit(84,'mm'),height=unit(110,'mm'),just=c('left','bottom')));heat();popViewport();pushViewport(viewport(x=unit(88,'mm'),y=unit(2,'mm'),width=unit(110,'mm'),height=unit(110,'mm'),just=c('left','bottom')));path();popViewport()}
cairo_pdf(file.path(o,'Fig5fg_compact.pdf'),width=198/25.4,height=112/25.4,family=font);draw();dev.off()
svglite(file.path(o,'Fig5fg_compact.svg'),width=198/25.4,height=112/25.4);draw();dev.off()
png(file.path(o,'Fig5fg_compact_300dpi.png'),width=198,height=112,units='mm',res=300,type='cairo');draw();dev.off()
# Vector-only PDF placement. All original panels are embedded without rerasterization.
pre<-function(h)c('\\documentclass{article}',sprintf('\\usepackage[paperwidth=595.276bp,paperheight=%sbp,margin=0bp]{geometry}',h),'\\usepackage{graphicx}','\\pagestyle{empty}','\\setlength{\\parindent}{0pt}','\\setlength{\\unitlength}{1bp}','\\begin{document}',sprintf('\\noindent\\begin{picture}(595.276,%s)',h))
place<-function(file,src,x,y,pageh,scale=1,sw=595.276,sh=841.89){w<-src[3]-src[1];h<-src[4]-src[2];sprintf('\\put(%.5f,%.5f){\\includegraphics[trim=%.5fbp %.5fbp %.5fbp %.5fbp,clip,width=%.5fbp]{%s}}',x,pageh-y-h*scale,src[1],sh-src[4],sw-src[3],src[2],w*scale,file)}
mainh<-910
main<-c(pre(mainh),place('../input_main.pdf',c(0,0,595.276,558.5),0,0,mainh),place('../input_main.pdf',c(348,558.5,595.276,576.5),348,558.5,mainh),sprintf('\\put(15,%.5f){\\includegraphics[width=%.5fbp]{Fig5fg_compact.pdf}}',mainh-582-112*72/25.4,198*72/25.4),'\\end{picture}','\\end{document}')
supph<-841.89
supp<-c(pre(supph),place('../input_supp.pdf',c(0,0,595.276,480),0,0,supph),place('../input_main.pdf',c(12,569.5,232,831),20,515,supph),place('../input_main.pdf',c(250,580,575,831),250,525.5,supph),sprintf('\\put(20,%.5f){\\sffamily\\fontsize{8}{9}\\selectfont\\textbf{c}}',supph-508),sprintf('\\put(250,%.5f){\\sffamily\\fontsize{8}{9}\\selectfont\\textbf{d}}',supph-508),'\\end{picture}','\\end{document}')
writeLines(main,file.path(o,'Fig5_main_MES.tex'));writeLines(supp,file.path(o,'FigS5_with_original_fg.tex'));old<-getwd();setwd(o)
for(n in c('Fig5_main_MES','FigS5_with_original_fg')){rc<-system2('/usr/bin/pdflatex',c('-interaction=nonstopmode','-halt-on-error',paste0(n,'.tex')),stdout=paste0(n,'_build.log'),stderr=paste0(n,'_error.log'));stopifnot(rc==0);system2('/usr/bin/pdftocairo',c('-png','-singlefile','-r','150',paste0(n,'.pdf'),n))}
system2('/usr/bin/pdftocairo',c('-png','-singlefile','-r','180','Fig5fg_compact.pdf','Fig5fg_compact'));setwd(old)
fwrite(D[pathway_id%in%rows$pathway_id],file.path(o,'f_source.tsv'),sep='\t');fwrite(F,file.path(o,'g_source.tsv'),sep='\t')
# Render geometry: same H and y origin for each UMAP and its three crop panels.
geom<-rbindlist(lapply(1:3,function(i)data.table(id=c(paste0('u',i),paste0('p',i,1:3)),row=i,x=c(92,136,157,178),y=c(70,39,8)[i]+2,w=c(25,19,19,19),h=25)))
stopifnot(all(geom[,diff(range(y)),by=row]$V1==0),all(geom$h==25));fwrite(geom,file.path(o,'panel_geometry_mm.tsv'),sep='\t')
writeLines(c('Style and placement changes only; 168 heatmap values, 14 soil pairs and 1708 cells preserved.','New main f/g: 198 x 112 mm, minimum new text size 5 pt. Main canvas height extended to 910 bp to retain a-e at original scale.','Original main f -> Supplement c; original main g -> Supplement d; Supplement a/b retained.','Original embedded panels retain their source fonts including pre-existing small text; no claim of submission compliance for inherited typography.','Original UMI-share panels are complementary descriptive analyses, not identical estimands or validation of MES.','Main g selected broad classes remain CV nucleotide / Halomonas amino acid / PA energy; displaced old g remains CV energy / Halomonas amino acid / PA nucleotide.','Vector placement through LaTeX graphicx driven by R; no scientific panel rasterization.'),file.path(o,'layout_notes.txt'))
cat('MAIN_SUPP_V8_COMPLETE\n')

