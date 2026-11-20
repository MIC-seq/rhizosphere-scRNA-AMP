.libPaths(c('/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/imeta_mes_20261025/Rlib','/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions/software/Rlib',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(GSVA);library(GSEABase)})
b<-'/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions';o<-file.path(b,'six_class_joint_mes_20261026')
for(d in c('tables','objects','figures','scripts','reports','source_data','provenance','logs'))dir.create(file.path(o,d),recursive=TRUE,showWarnings=FALSE)
source(file.path(b,'imeta_mes_20261025/author_code_only/scripts/author_functions_unmodified.R'))
Oall<-readRDS(file.path(b,'inputs/analysis_objects.rds'));catg<-fread(file.path(b,'broad_metabolism_20261018/tables/KEGG_class_pathway_catalogue.tsv'));U<-fread(file.path(b,'joint_genome_20261026/tables/joint_UMAP_all_selected_cells.tsv'))
cl<-c('Amino acid metabolism','Carbohydrate metabolism','Energy metabolism','Lipid metabolism','Nucleotide metabolism','Metabolism of cofactors and vitamins');ids<-c('amino_acid','carbohydrate','energy','lipid','nucleotide','cofactor');names(cl)<-ids
stopifnot(all(cl%in%catg$group),nrow(U)==1708,uniqueN(U$cell_id)==1708,all(is.finite(U$joint_UMAP_1)))
full<-cells<-members<-audit<-list()
for(si in seq_along(Oall)){
 sp<-names(Oall)[si];O<-Oall[[si]];M<-as.matrix(O$counts);G<-as.data.table(O$genes);C<-as.data.table(O$cells);rownames(M)<-G$gene_id;colnames(M)<-C$cell_id
 stopifnot(!anyDuplicated(rownames(M)),!anyDuplicated(colnames(M)),setequal(C$cell_id,U[species==sp,cell_id]))
 gp<-lapply(G$KEGG_Pathway,function(x)unique(sub('^ko','map',trimws(strsplit(x,',',fixed=TRUE)[[1]]))))
 ix<-setNames(lapply(cl,function(g){pp<-catg[group==g,pathway_id];which(vapply(gp,function(p)any(p%in%pp),TRUE))}),ids)
 stopifnot(all(lengths(ix)>0));sets<-GeneSetCollection(lapply(ids,function(k)GeneSet(G$gene_id[ix[[k]]],setName=k,geneIdType=NullIdentifier())))
 for(k in ids){j<-ix[[k]];members[[paste(si,k)]]<-cbind(data.table(species=sp,class_id=k,class=cl[[k]]),G[j,setdiff(names(G),'species'),with=FALSE],data.table(class_gene_UMI=rowSums(M[j,,drop=FALSE]),class_gene_detected_cells=rowSums(M[j,,drop=FALSE]>0)))}
 cache<-file.path(o,'objects',paste0('species_',si,'_six_classes.rds'));seed<-2026102600L+si
 if(file.exists(cache)){z<-readRDS(cache);raw<-z$raw;null<-z$random_scores;norm<-z$norm;random_sets<-z$random_sets}else{
 set.seed(seed);cat(format(Sys.time()),'START',sp,'\n');flush.console()
 raw<-raw_mes(geneset_input=sets,count_matrix=M,method='ssGSEA');ng<-get_gene_count(geneset_input=sets);random_sets<-get_random_geneset(pathway_gsva=raw,gene_count=ng,count_matrix=M)
 jobs<-parallel::mclapply(seq_len(nrow(raw)),function(ii){rr<-get_permutation(pathway_gsva=raw[ii,,drop=FALSE],count_matrix=M,total_pathway_random=random_sets);cat('RANDOM_CLASS_COMPLETE',sp,rownames(raw)[ii],'\n');flush.console();rr[[1]]},mc.cores=6,mc.set.seed=FALSE)
 stopifnot(!any(vapply(jobs,inherits,FALSE,'try-error')));null<-setNames(jobs,rownames(raw));norm<-author_normalize(raw,null)
 saveRDS(list(species=sp,seed=seed,raw=raw,random_scores=null,random_sets=random_sets,norm=norm,gene_sets=sets,cell_ids=C$cell_id,gene_ids=G$gene_id),cache)
 }
 for(k in ids){j<-ix[[k]];v<-colSums(M[j,,drop=FALSE]);c<-merge(C,U[,.(cell_id,joint_UMAP_1,joint_UMAP_2,joint_cluster)],by='cell_id',all.x=TRUE,sort=FALSE);c<-c[match(colnames(M),cell_id)]
  c[,`:=`(class_id=k,class=cl[[k]],pathway_UMI=as.numeric(v),raw_ssgsea=as.numeric(raw[k,]),observed=as.numeric(v)>0)]
  cells[[paste(si,k)]]<-c
  full[[paste(si,k)]]<-data.table(species=sp,class_id=k,class=cl[[k]],raw_mean=mean(as.numeric(raw[k,])),MES=norm$MES[k,1],P_background=norm$P[k,1],cells=ncol(M),member_genes=length(j),detected_genes=sum(rowSums(M[j,,drop=FALSE])>0),detected_cells=sum(v>0),detection_fraction=mean(v>0),pathway_UMI=sum(v),detected_soils=uniqueN(c$sample[v>0]),soils=uniqueN(c$sample),seed=seed,n_random=1000)
 }
 cat(format(Sys.time()),'SPECIES_COMPLETE',sp,'\n');flush.console()
}
D<-rbindlist(full);D[,q_background_BH:=p.adjust(P_background,'BH')];D[,q_background_BH6:=p.adjust(P_background,'BH'),by=species]
S<-rbindlist(cells,fill=TRUE);J<-rbindlist(members,fill=TRUE)
fwrite(D,file.path(o,'tables/six_class_species_MES.tsv'),sep='\t');fwrite(S,file.path(o,'tables/six_class_cell_raw_scores.tsv'),sep='\t');fwrite(J,file.path(o,'tables/six_class_gene_membership.tsv'),sep='\t');fwrite(catg[group%in%cl],file.path(o,'tables/six_class_pathway_catalogue.tsv'),sep='\t')
fwrite(S[,.(cells=.N,detected_cells=sum(observed),pathway_UMI=sum(pathway_UMI),mean_raw=mean(raw_ssgsea)),by=.(species,sample,crop,class_id,class)],file.path(o,'tables/six_class_soil_raw_profiles.tsv'),sep='\t')
fwrite(D[,.(species,class_id,detected_cells,member_genes,detected_genes,detected_soils,MES,P_background)],file.path(o,'reports/preselection_coverage.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(o,'reports/MES_sessionInfo.txt'))
cat('SIX_CLASS_AUTHOR_MES_COMPLETE\n');print(D)
