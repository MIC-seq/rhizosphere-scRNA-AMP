#!/usr/bin/env Rscript
# Reference phylogeny with family sectors, phylum band and pooled transcript-read abundance.
# All taxa passing documented filters are shown; only tip labels are thinned for readability.
suppressPackageStartupMessages(library(ape))
.libPaths(c('software/Rlib',.libPaths()))
options(stringsAsFactors=FALSE)
tab <- read.delim('source_data/tree_species.tsv',check.names=FALSE)
tr <- read.tree('tree/species_phylogeny.nwk')
tr <- ladderize(tr)
tr$tip.label <- sub('\\.(fna|faa)$','',tr$tip.label)
stopifnot(!anyDuplicated(tr$tip.label),all(tr$tip.label %in% tab$tip))
excluded <- tab[!tab$tip %in% tr$tip.label,,drop=FALSE]
excluded$tree_exclusion <- rep('not_retained_after_phylogenetic_marker_QC',nrow(excluded))
write.table(excluded,'source_data/marker_qc_exclusions.tsv',sep='\t',row.names=FALSE,quote=FALSE)
tab <- tab[match(tr$tip.label,tab$tip),]
stopifnot(!anyNA(tab$tip),all(is.finite(tab$pooled_percent)),all(tab$pooled_percent>0))
n <- Ntip(tr); root <- setdiff(tr$edge[,1],tr$edge[,2]);stopifnot(length(root)==1)
children <- split(tr$edge[,2],tr$edge[,1])
desc <- function(node) if(node<=n) node else unlist(lapply(children[[as.character(node)]],desc))
ord <- desc(root);pos <- integer(n);pos[ord] <- seq_len(n)
# Eight-degree opening; branch lengths remain proportional to substitutions per site.
tipangle <- (94+(pos-0.5)*352/n)*pi/180
angle <- numeric(n+tr$Nnode);angle[seq_len(n)] <- tipangle
for(node in (n+1):(n+tr$Nnode)) angle[node] <- mean(range(tipangle[desc(node)]))
depth <- node.depth.edgelength(tr);stopifnot(max(depth)>0)
radius <- 0.10+0.81*depth/max(depth)
palette <- c('#4C9FBA','#E2AC47','#67AE9C','#B49BC7','#DC8585','#7F9CB7','#B5B96E','#CC9B77','#A0B8B5','#BA99AD','#98A4C0','#86A788')
ph <- names(sort(tapply(tab$reads,tab$phylum,sum),decreasing=TRUE))
phcols <- setNames(rep(palette,length.out=length(ph)),ph)
fam <- names(sort(tapply(tab$reads,tab$family,sum),decreasing=TRUE))
fam <- head(setdiff(fam,'Unclassified'),10)
famcols <- setNames(adjustcolor(rep(palette,length.out=length(fam)),alpha.f=0.24),fam)
famcols <- c(famcols,Other='#F2F2F2')
group <- ifelse(tab$family %in% fam,tab$family,'Other')
xy <- function(r,a) cbind(r*cos(a),r*sin(a))
arc <- function(r,a,b,col,lwd=0.4,lty=1) {
  aa <- seq(a,b,length.out=max(3,ceiling(abs(b-a)*100)))
  lines(xy(r,aa),col=col,lwd=lwd,lty=lty)
}
sector <- function(r0,r1,a,b,col) {
  aa <- seq(a,b,length.out=max(3,ceiling(abs(b-a)*70)))
  polygon(rbind(xy(r0,aa),xy(r1,rev(aa))),col=col,border=NA)
}
bar0 <- 1.035; barspan <- .22
# Log10 axis starts at the a priori abundance inclusion threshold (0.01%).
lo <- -2; hi <- max(-1,ceiling(log10(max(tab$pooled_percent))))
barend <- bar0+barspan*(log10(tab$pooled_percent)-lo)/(hi-lo)
da <- 352*pi/180/n
selected <- head(order(tab$reads,decreasing=TRUE),12)
draw <- function(){
 par(mar=c(.5,.5,1,.3),family='Helvetica',ps=7,xpd=NA)
 plot.new();plot.window(xlim=c(-1.66,2.5),ylim=c(-1.55,1.6),asp=1)
 # Contiguous family sectors; disjoint occurrences remain separate sectors.
 g <- group[ord]; rr <- rle(g); end <- cumsum(rr$lengths);start <- c(1,head(end,-1)+1)
 for(k in seq_along(start)){
   ids <- ord[start[k]:end[k]];r0 <- .12
   if(length(ids)>1){m <- getMRCA(tr,ids);if(!is.null(m))r0 <- max(.12,radius[m])}
   sector(r0,.965,min(tipangle[ids])-da/2,max(tipangle[ids])+da/2,famcols[rr$values[k]])
 }
 for(node in (n+1):(n+tr$Nnode)) {
   kids <- children[[as.character(node)]]
   arc(radius[node],min(angle[kids]),max(angle[kids]),'#4E5357',.35)
 }
 for(k in seq_len(nrow(tr$edge))){
   p <- tr$edge[k,1];ch <- tr$edge[k,2]
   z <- xy(c(radius[p],radius[ch]),angle[ch]);lines(z,col='#4E5357',lwd=.35)
 }
 for(i in seq_len(n)){
   z <- xy(c(radius[i],.967),tipangle[i]);lines(z,col='#D6D9DB',lwd=.22,lty=3)
   sector(.975,1.017,tipangle[i]-da/2,tipangle[i]+da/2,phcols[tab$phylum[i]])
 }
 for(tick in seq(lo,hi))arc(bar0+barspan*(tick-lo)/(hi-lo),94*pi/180,446*pi/180,'#E3E7EA',.3)
 for(i in seq_len(n)) sector(bar0,barend[i],tipangle[i]-.35*da,tipangle[i]+.35*da,'#5A93B4')
 # Direct labels use separate left/right vertical placement to avoid collisions.
 for(side in c(-1,1)){
   ids <- selected[sign(cos(tipangle[selected]))==side]
   ids <- ids[order(sin(tipangle[ids]))]
   if(!length(ids))next
   yy <- 1.35*sin(tipangle[ids]);gap <- .095
   for(j in seq_along(yy))if(j>1)yy[j] <- max(yy[j],yy[j-1]+gap)
   if(max(yy)>1.38)yy <- yy-(max(yy)-1.38)
   for(j in seq_along(ids)){
     i <- ids[j];z <- xy(max(1.27,barend[i]+.02),tipangle[i]);xx <- side*1.38
     lines(rbind(z,c(side*1.32,yy[j]),c(xx,yy[j])),col=phcols[tab$phylum[i]],lwd=.45)
     text(xx,yy[j],tab$species[i],adj=c(if(side==1)0 else 1,.5),cex=.74,font=3,col='#34393E')
   }
 }
 # Compact scale in the opening at 90 degrees.
 for(tick in seq(lo,hi)){
   rt <- bar0+barspan*(tick-lo)/(hi-lo)
   text(.018,rt,format(10^tick,scientific=FALSE,trim=TRUE),adj=c(0,.5),cex=.72)
 }
 text(-.02,1.44,'Relative transcript-read abundance (%)',adj=c(.5,.5),cex=.72)
 text(-1.64,1.54,'b',adj=0,font=2,cex=1.25)
 text(-1.45,1.54,sprintf('Reference phylogeny of %s detected species',n),adj=0,font=2,cex=.9)
 # Legends use ranked categories; remaining families retain all their tips in neutral shading.
 lx <- 1.84;ly <- 1.18
 text(lx,ly,'Phylum',adj=0,font=2,cex=.85);ly <- ly-.10
 for(p in ph){rect(lx,ly-.026,lx+.06,ly+.026,col=phcols[p],border=NA);text(lx+.09,ly,p,adj=0,cex=.74);ly <- ly-.075}
 ly <- ly-.045;text(lx,ly,'Family',adj=0,font=2,cex=.85);ly <- ly-.10
 for(f in c(fam,'Other')){rect(lx,ly-.026,lx+.06,ly+.026,col=famcols[f],border=NA);text(lx+.09,ly,f,adj=0,cex=.74);ly <- ly-.075}
 text(lx,ly-.08,'49 datasets pooled',adj=0,cex=.72)
 # Branch-length scale is exactly proportional to radial branch lengths.
 scale <- signif(max(depth)/5,1);x0 <- -.27;y0 <- -1.45
 segments(x0,y0,x0+.81*scale/max(depth),y0,lwd=.8)
 text(x0,y0-.07,paste(scale,'substitutions/site'),adj=0,cex=.74)
}
dir.create('figures',showWarnings=FALSE)
# 183 mm main-figure width; all drawing/exporting performed by R.
cairo_pdf('figures/Figure1b_phylogeny.pdf',width=183/25.4,height=145/25.4,family='Helvetica');draw();dev.off()
svglite::svglite('figures/Figure1b_phylogeny.svg',width=183/25.4,height=145/25.4);draw();dev.off()
ragg::agg_png('figures/Figure1b_phylogeny.png',width=183,height=145,units='mm',res=600);draw();dev.off()
ragg::agg_tiff('figures/Figure1b_phylogeny.tiff',width=183,height=145,units='mm',res=600,compression='lzw');draw();dev.off()
write.table(tab,'source_data/plotted_species.tsv',sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c(sprintf('Fig. 1b | Reference-genome phylogeny of %d prokaryotic species detected across 49 agricultural-soil microbial scRNA-seq datasets.',n),
 'All 49 datasets were treated as a pooled metatranscriptomic collection. Only non-rRNA R2 reads were classified; barcode/UMI reads were not used as transcript sequences.',
 'Family background sectors and the phylum ring annotate reference taxonomy. The outer bars show pooled relative transcript-read abundance on a log10 scale, not cellular abundance. The denominator includes all reads classified to Bacteria or Archaea before species filtering; retained tips are not renormalized.',
 'Species required at least 100 reads and 0.01% pooled prokaryotic transcript-read abundance, a suitable local reference genome and at least 100 phylogenetic markers. All exclusions and retained abundance coverage are recorded in source data.',
 'Proteins were predicted consistently from all reference genomes using Prodigal2.6.3 (single-genome mode, genetic code11). PhyloPhlAn 400-marker analysis used DIAMOND blastp (highest-scoring HSP per protein-target pair), MAFFT, trimAl and IQ-TREE (LG+G4; 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates; seed 42). Branch lengths are substitutions per site. The circular representation does not assert an ancestral root. Only the 12 most abundant retained species are labelled; all retained taxa are plotted.'), 'figures/Figure1b_legend.txt')
coverage <- sum(tab$pooled_percent)
writeLines(c(paste('Tree tips:',n),paste('Retained abundance coverage (% of prokaryotic classified reads):',coverage),paste('Excluded after marker QC:',nrow(excluded)),'Status: generated; visual review required before publication.'),'figures/QA.txt')
