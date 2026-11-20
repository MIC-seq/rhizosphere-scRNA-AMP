base<-'/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions';prev<-file.path(base,'imeta_mes_20261025');out<-file.path(prev,'author_code_only');for(d in c('tables','objects','figures','scripts','reports','provenance'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(prev,'Rlib'),file.path(base,'software/Rlib'),.libPaths()));suppressPackageStartupMessages({library(data.table);library(GSVA);library(GSEABase)})
source(file.path(out,'scripts/author_functions_unmodified.R'))
Oall<-readRDS(file.path(base,'inputs/analysis_objects.rds'));catg<-fread(file.path(base,'inputs/pathway_catalogue.tsv'))
args<-commandArgs(TRUE);mode<-if(length(args))args[1]else'overall';allresults<-list();audit<-list()
selection<-if(length(args)>=2)as.integer(args[2])else seq_along(Oall);workers<-if(length(args)>=3)as.integer(args[3])else 6L
for(si in selection){
 s<-names(Oall)[si];O<-Oall[[si]];M<-as.matrix(O$counts);G<-as.data.table(O$genes);C<-as.data.table(O$cells);rownames(M)<-G$gene_id;colnames(M)<-C$cell_id
 sets<-GeneSetCollection(lapply(names(O$pathway_indices),function(f)GeneSet(G$gene_id[O$pathway_indices[[f]]],setName=f,geneIdType=NullIdentifier())))
 K<-fread(file.path(base,'inputs/state_key.tsv'));C<-merge(C,K[,.(species,functional_state,state_code)],by=c('species','functional_state'),all.x=TRUE,sort=FALSE);C<-C[match(O$cells$cell_id,cell_id)]
 units<-if(mode=='overall')'all'else if(mode=='state')sort(unique(C$state_code))else c('soybean','rice','wheat')
 for(u in units){
  j<-if(u=='all')seq_len(ncol(M))else if(mode=='state')which(C$state_code==u)else which(C$crop==u);tag<-paste0('species_',si,'_',u);cat(format(Sys.time()),'START',tag,'cells',length(j),'\n');flush.console()
  if(length(j)<2){audit[[tag]]<-data.table(species=s,unit=u,cells=length(j),status='No expression matrix with >=2 cells');next}
  cache<-file.path(out,'objects',paste0(tag,'.rds'))
  if(file.exists(cache)){r<-readRDS(cache);allresults[[tag]]<-r$result;audit[[tag]]<-r$audit;next}
  X<-M[,j,drop=FALSE];seed<-20261025L+si*100L+match(u,c('all','soybean','rice','wheat','S1','S2','S3'));set.seed(seed);t0<-Sys.time()
  gsva_es<-raw_mes(geneset_input=sets,count_matrix=X,method='ssGSEA')
  gmt_gene_count<-get_gene_count(geneset_input=sets)
  random_set<-get_random_geneset(pathway_gsva=gsva_es,gene_count=gmt_gene_count,count_matrix=X)
  # Dispatch independent pathways in parallel; original function and 1000-set matrices are unchanged.
  jobs<-parallel::mclapply(seq_len(nrow(gsva_es)),function(ii){zz<-get_permutation(pathway_gsva=gsva_es[ii,,drop=FALSE],count_matrix=X,total_pathway_random=random_set);cat('PATH_DONE',tag,rownames(gsva_es)[ii],'\\n');flush.console();zz[[1]]},mc.cores=workers,mc.set.seed=FALSE)
  stopifnot(!any(vapply(jobs,inherits,FALSE,'try-error')));first_permutation<-setNames(jobs,rownames(gsva_es));rm(jobs)
  norm<-author_normalize(gsva_es,first_permutation)
  res<-data.table(species=s,unit=u,pathway_id=rownames(gsva_es),raw_mean=rowMeans(gsva_es),MES=norm$MES[,1],P=norm$P[,1],seed=seed,n_random=1000,cells=length(j),genes=nrow(X))
  res[,`:=`(author_activated=MES>0,author_significant=P<.01)]
  res[,pathway_UMI:=vapply(pathway_id,function(f)sum(X[O$pathway_indices[[f]],,drop=FALSE]),0.0)]
  res[,detected_cells:=vapply(pathway_id,function(f)sum(colSums(X[O$pathway_indices[[f]],,drop=FALSE])>0),0.0)]
  res<-merge(res,catg,by='pathway_id',all.x=TRUE,sort=FALSE)
  aa<-data.table(species=s,unit=u,cells=length(j),genes=nrow(X),pathways=nrow(gsva_es),status='complete',seconds=as.numeric(difftime(Sys.time(),t0,units='secs')))
  saveRDS(list(result=res,audit=aa,raw_cell_scores=gsva_es,random_sets=random_set,random_scores=first_permutation),cache,compress=TRUE)
  fwrite(res,file.path(out,'tables',paste0(tag,'.tsv')),sep='\t');allresults[[tag]]<-res;audit[[tag]]<-aa
  cat(format(Sys.time()),'DONE',tag,'seconds',aa$seconds,'\n');flush.console();rm(first_permutation,random_set);gc()
 }
}
fwrite(rbindlist(allresults,fill=TRUE),file.path(out,'tables',paste0(mode,if(length(selection)==1)paste0('_species_',selection)else'', '_MES.tsv')),sep='\t');fwrite(rbindlist(audit,fill=TRUE),file.path(out,'reports',paste0(mode,if(length(selection)==1)paste0('_species_',selection)else'', '_run_audit.tsv')),sep='\t');capture.output(sessionInfo(),file=file.path(out,'reports/sessionInfo.txt'))
