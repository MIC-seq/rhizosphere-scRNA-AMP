#!/usr/bin/env Rscript
# R-only figure; calibrated display change, not a change in enrichment estimates.
# Core claim: species-specific AMP programs differ in representation across fixed COG states.
# All 88 cells of the matrix retained: 72 finite values plus 16 zero-observation cells.
# Final size 75 x 80 mm; serif-free vector PDF/SVG and 600-dpi PNG.
args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=', '', grep('^--file=',commandArgs(),value=TRUE)[1])
base <- dirname(normalizePath(script))
out <- if(length(args)) args[1] else file.path(base,'rendered')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table);library(jsonlite);library(scales)})
source(file.path(base,'panel_alignment.R'))
d <- fread(file.path(base,'old_vs_corrected_all_88_tiles.tsv'))
stopifnot(nrow(d)==88,sum(d$total_cells==0)==16,all(d[total_cells==0,enrichment_ratio]==0),all(is.finite(d[total_cells>0,log2_enrichment])))
font <- 'Liberation Sans'
theme_set(theme_void(base_family=font))
sp <- c('Chloracidobacterium validum','Halomonas sp. JS92-SW72','Pseudomonas aeruginosa')
states <- list(c('Envelope/signaling','Motility/trafficking'),c('Transcription/signaling','Replication/metabolism','Ion transport'),c('Turnover/signaling','Cell cycle/translation'))
labels <- list(c('S1 [M]','S2 [N]'),c('S1 [K]','S2 [L]','S3 [P]'),c('S1 [O]','S2 [D]'))
titles <- c("italic('C. validum')","italic('Halomonas')~'sp. JS92-SW72'","italic('P. aeruginosa')")
zero_grey <- '#D4D4D4'
low <- '#357493'; mid <- '#FFFFFF'; high <- '#B85C4A'
colour_scale_s <- 0.2
nonlinear <- scales::trans_new(name='asinh_log2_s0.2',transform=function(x)asinh(x/colour_scale_s),inverse=function(y)colour_scale_s*sinh(y),domain=c(-Inf,Inf))
make_panel <- function(i,mode='asinh'){
 z <- copy(d[species==sp[i]])
 programs <- unique(z$AMP_program)
 programs <- programs[order(programs=='Other',suppressWarnings(as.integer(sub('Program ','',programs))))]
 z[,`:=`(x=match(AMP_program,programs),y=match(functional_state,rev(states[[i]])))]
 # A category-specific grey is used only for zero observed cells.
 # Finite values retain their original log2 units and are never replaced by zero.
 z[,display:=ifelse(total_cells==0,NA_real_,log2_enrichment)]
 selected_trans <- if(mode=='asinh') nonlinear else 'identity'
 selected_breaks <- if(mode=='asinh') c(-2.5,-.2,0,.2,2.5) else c(-2.5,0,2.5)
 p <- ggplot(z,aes(x,y))+geom_tile(aes(fill=display),width=.96,height=.94,colour='white',size=.18)+
  geom_tile(data=z[total_cells==0],aes(alpha='No observed cells'),fill=zero_grey,width=.96,height=.94,colour='white',size=.18,show.legend=TRUE)+
  scale_fill_gradient2(low=low,mid=mid,high=high,midpoint=0,limits=c(-2.5,2.5),na.value=zero_grey,trans=selected_trans,breaks=selected_breaks,labels=as.character(selected_breaks),name=if(mode=='asinh')'Log2 enrichment (asinh colour scale)' else 'Log2 enrichment (linear colour scale)')+
  scale_alpha_manual(values=c('No observed cells'=1),name=NULL,guide=guide_legend(order=2,override.aes=list(fill=zero_grey,colour='white',alpha=1)))+
  scale_x_continuous(breaks=seq_along(programs),labels=sub('Program ','',programs),expand=expansion(add=.51))+
  scale_y_continuous(breaks=seq_along(states[[i]]),labels=rev(labels[[i]]),expand=expansion(add=.51))+
  labs(x=NULL,y=NULL,title=parse(text=titles[i]))+
  theme_void(base_size=6,base_family=font)+
  theme(axis.text.x=element_text(size=5.2,colour='#333333',margin=margin(t=2)),axis.text.y=element_text(size=5.6,colour='#333333',margin=margin(r=3)),plot.title=element_text(size=6.6,hjust=0,margin=margin(b=3)),legend.title=element_text(size=5.8),legend.text=element_text(size=5.2),legend.key.width=unit(2.4,'mm'),legend.key.height=unit(2.4,'mm'),legend.spacing.y=unit(.5,'mm'),legend.box='vertical',plot.margin=margin(2,2,3,2),legend.position='bottom')+
  guides(fill=guide_colorbar(order=1,title.position='top',title.hjust=.5,barwidth=unit(42,'mm'),barheight=unit(2.1,'mm')))
 p
}
make_stack <- function(mode='asinh') wrap_plots(lapply(1:3,make_panel,mode=mode),ncol=1,heights=c(2,3,2),guides='collect')+
 plot_annotation(title='e',caption='AMP-expression program ID',theme=theme(plot.title=element_text(family=font,size=8,face='bold'),plot.caption=element_text(family=font,size=5.3,hjust=.5,margin=margin(t=1)),plot.margin=margin(2,2,2,2))) & theme(legend.position='bottom',legend.box='vertical')
export <- function(p,name,w=75,h=80){
 prefix <- file.path(out,name)
 require_patchwork_panel_alignment(p,paste0(prefix,'.alignment-layout.json'),paste0(prefix,'.alignment.json'),w/25.4,h/25.4,panel_ids=c('C','H','P'),column_groups=list(c('C','H','P')),audit_script=file.path(base,'audit_panel_alignment.py'),python=Sys.getenv('FIGURE_QA_PYTHON',Sys.which('python')),strict=TRUE)
 svglite::svglite(paste0(prefix,'.svg'),width=w/25.4,height=h/25.4,system_fonts=list(sans=font));print(p);dev.off()
 grDevices::cairo_pdf(paste0(prefix,'.pdf'),width=w/25.4,height=h/25.4,family=font);print(p);dev.off()
 grDevices::png(paste0(prefix,'_600dpi.png'),width=w,height=h,units='mm',res=600,type='cairo');print(p);dev.off()
}
export(make_stack('asinh'),'Fig5e_corrected_asinh_main')
export(make_stack('linear'),'Fig5e_corrected_linear_reference')
fwrite(d,file.path(out,'Fig5e_Source_Data_UNCHANGED.tsv'),sep='\t',na='NA')
settings <- list(statistic='Original corrected log2 enrichment unchanged',display_only='asinh(log2_enrichment / 0.2)',shared_across_all_species=TRUE,finite_colour_limits=c(-2.5,2.5),legend_tick_labels_are_original_log2_units=TRUE,zero_observations='Distinct grey, no in-tile text',n_zero_observations=16,n_finite_values=72)
write_json(settings,file.path(out,'colour_mapping.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
message('Server R render complete. Inspect the asinh main panel and linear reference at actual placement size, then run PDF font/collision checks.')
