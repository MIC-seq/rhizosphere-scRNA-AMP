#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Matrix);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=3)stop("root config outdir")
root<-normalizePath(args[1]);cfg<-fromJSON(args[2]);out<-args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(as.integer(cfg$embedding$random_seed))
species_order<-c("Chloracidobacterium validum","Halomonas sp. JS92-SW72","Pseudomonas aeruginosa")
species_short<-setNames(c("C. validum","Halomonas sp. JS92-SW72","P. aeruginosa"),species_order)
slug<-function(x)gsub("[^A-Za-z0-9_.-]+","_",x)
state_map<-list(
 "Chloracidobacterium validum"=c(`0`="Envelope/signaling",`1`="Motility/trafficking"),
 "Halomonas sp. JS92-SW72"=c(`0`="Transcription/signaling",`1`="Replication/metabolism",`2`="Ion transport"),
 "Pseudomonas aeruginosa"=c(`0`="Turnover/signaling",`1`="Cell cycle/translation"))
state_order<-c("Envelope/signaling","Motility/trafficking","Transcription/signaling","Replication/metabolism","Ion transport","Turnover/signaling","Cell cycle/translation")
state_cols<-setNames(c("#4477AA","#EE6677","#228833","#CCBB44","#66CCEE","#AA3377","#BBBBBB"),state_order)

cells<-rbindlist(lapply(species_order,function(sp){
 d<-fread(file.path(root,"results/12_integrated_amp_function/species_embeddings",slug(sp),"cog_umap_cells.tsv"))
 d[,`:=`(species=sp,functional_state=unname(state_map[[sp]][as.character(seurat_clusters)]),AMP_cluster=as.character(AMP_cluster))]
 d
}),fill=TRUE)
min_cells<-as.integer(cfg$association$minimum_major_amp_program_cells)
program_counts<-cells[!is.na(AMP_cluster)&AMP_cluster!="",.N,by=AMP_cluster]
major<-program_counts[N>=min_cells,AMP_cluster]
cells[,AMP_program:=fifelse(is.na(AMP_cluster)|AMP_cluster=="","Unclassified",fifelse(AMP_cluster%in%major,paste("Program",AMP_cluster),"Other"))]
program_levels<-c(paste("Program",sort(as.integer(major))),"Other","Unclassified")
program_levels<-program_levels[program_levels%in%cells$AMP_program]
base_program_cols<-setNames(hcl.colors(max(3,length(major)),"Dark 3"),paste("Program",sort(as.integer(major))))
program_cols<-c(base_program_cols[intersect(program_levels,names(base_program_cols))],Other="#C9C9C9",Unclassified="#E6E6E6");program_cols<-program_cols[names(program_cols)%in%program_levels]

# Leave-one-sample-out projection: infer the functional-state preference of each
# existing AMP program from all other samples, without changing either clustering.
cells[,AMP_inferred_state:="Unassigned"];projection_maps<-list()
for(sp in species_order) for(heldout in unique(cells[species==sp,sample])){
 train<-cells[species==sp&sample!=heldout&!is.na(functional_state)&!is.na(AMP_program)]
 states<-unique(train$functional_state);programs<-unique(train$AMP_program)
 grid<-CJ(AMP_program=programs,functional_state=states)
 obs<-train[,.N,by=.(AMP_program,functional_state)]
 grid<-merge(grid,obs,by=c("AMP_program","functional_state"),all.x=TRUE);grid[is.na(N),N:=0L]
 grid[,program_total:=sum(N),by=AMP_program];grid[,state_total:=sum(N),by=functional_state]
 total<-nrow(train);ns<-length(states)
 grid[,preference_score:=log2(((N+.5)/(program_total+.5*ns))/((state_total+.5)/(total+.5*ns)))]
 best<-grid[order(AMP_program,-preference_score),.SD[1],by=AMP_program]
 best[,`:=`(species=sp,heldout_sample=heldout)]
 projection_maps[[paste(sp,heldout)]]<-best
 idx<-which(cells$species==sp&cells$sample==heldout)
 pred<-setNames(best$functional_state,best$AMP_program)[cells$AMP_program[idx]]
 cells[idx,AMP_inferred_state:=fifelse(is.na(pred),"Unassigned",unname(pred))]
}
fwrite(rbindlist(projection_maps),file.path(out,"AMP_program_leave_one_sample_out_state_map.tsv"),sep="\t")
fwrite(cells,file.path(out,"figure5_cells_with_state_program.tsv.gz"),sep="\t",compress="gzip")
fwrite(data.table(functional_state=names(state_cols),color=unname(state_cols)),file.path(out,"functional_state_palette.tsv"),sep="\t")
fwrite(data.table(AMP_program=names(program_cols),color=unname(program_cols)),file.path(out,"AMP_program_palette.tsv"),sep="\t")
fwrite(program_counts[order(-N)],file.path(out,"AMP_program_original_counts.tsv"),sep="\t")

cramers_v<-function(a,b){tab<-table(a,b);n<-sum(tab);if(n==0||min(dim(tab))<2)return(NA_real_);sqrt(unname(suppressWarnings(chisq.test(tab,correct=FALSE)$statistic))/(n*min(nrow(tab)-1,ncol(tab)-1)))}
nmi<-function(a,b){tab<-table(a,b);p<-tab/sum(tab);pa<-rowSums(p);pb<-colSums(p);nz<-which(p>0,arr.ind=TRUE);mi<-sum(p[nz]*log(p[nz]/(pa[nz[,1]]*pb[nz[,2]])));ha<--sum(pa[pa>0]*log(pa[pa>0]));hb<--sum(pb[pb>0]*log(pb[pb>0]));if(ha+hb==0)NA_real_ else 2*mi/(ha+hb)}
ari<-function(a,b){tab<-table(a,b);choose2<-function(x)x*(x-1)/2;n<-sum(tab);idx<-sum(choose2(tab));ra<-sum(choose2(rowSums(tab)));cb<-sum(choose2(colSums(tab)));expected<-ra*cb/choose2(n);den<-(ra+cb)/2-expected;if(den==0)NA_real_ else (idx-expected)/den}
B<-as.integer(cfg$association$stratified_permutations);metrics<-list();cont_long<-list()
for(sp in species_order){
 d<-cells[species==sp&!is.na(functional_state)&!is.na(AMP_program)]
 tab<-as.data.table(table(d$functional_state,d$AMP_program));setnames(tab,c("functional_state","AMP_program","N"));tab[,species:=sp]
 fwrite(tab,file.path(out,paste0("contingency_",slug(sp),".tsv")),sep="\t");cont_long[[sp]]<-tab
 obs<-cramers_v(d$functional_state,d$AMP_program)
 perm<-rep(NA_real_,B)
 for(i in seq_len(B)){pp<-d$AMP_program;for(ix in split(seq_len(nrow(d)),d$sample))pp[ix]<-sample(pp[ix]);perm[i]<-cramers_v(d$functional_state,pp)}
 metrics[[sp]]<-data.table(species=sp,n_cells=nrow(d),n_states=uniqueN(d$functional_state),n_programs=uniqueN(d$AMP_program),cramers_v=obs,NMI=nmi(d$functional_state,d$AMP_program),ARI=ari(d$functional_state,d$AMP_program),permutations=B,permutation_p=(1+sum(perm>=obs,na.rm=TRUE))/(1+sum(is.finite(perm))),permuted_v_mean=mean(perm,na.rm=TRUE),permuted_v_sd=sd(perm,na.rm=TRUE))
}
metrics<-rbindlist(metrics);fwrite(metrics,file.path(out,"state_program_association_metrics.tsv"),sep="\t");fwrite(rbindlist(cont_long),file.path(out,"state_program_contingency_long.tsv"),sep="\t")

# E: calculate enrichment within each sample, then average samples equally.
e_sample<-cells[!is.na(functional_state)&!is.na(AMP_program),{
 fs<-functional_state;ap<-AMP_program;sample_n<-.N;states<-unique(fs);programs<-unique(ap);z<-CJ(functional_state=states,AMP_program=programs)
 z[,N:=mapply(function(st,pr)sum(fs==st&ap==pr),functional_state,AMP_program)]
 z[,state_total:=vapply(functional_state,function(st)sum(fs==st),numeric(1))]
 z[,program_total:=vapply(AMP_program,function(pr)sum(ap==pr),numeric(1))]
 z[,`:=`(state_fraction=N/state_total,log2_enrichment=log2(((N+.5)/(state_total+1))/((program_total+.5)/(sample_n+1))))]
 z
},by=.(species,sample)]
e_summary<-e_sample[,.(mean_state_fraction=mean(state_fraction),mean_log2_enrichment=mean(log2_enrichment),sd_log2_enrichment=sd(log2_enrichment),samples=.N,total_cells=sum(N)),by=.(species,functional_state,AMP_program)]
e_summary<-merge(e_summary,metrics[,.(species,association_p=permutation_p)],by="species",all.x=TRUE)
fwrite(e_sample,file.path(out,"figure5E_within_sample_enrichment.tsv"),sep="\t");fwrite(e_summary,file.path(out,"figure5E_enrichment_summary.tsv"),sep="\t")

# F: log-normalized COG scores; z-score within each sample and COG, then average samples equally.
read_sparse<-function(path){con<-gzfile(path,"rt");on.exit(close(con));as(readMM(con),"dgCMatrix")}
f_all<-list();f_sample_all<-list()
for(sp in species_order){
 sd<-file.path(root,"results/12_integrated_amp_function/species_matrices",slug(sp));meta<-fread(file.path(sd,"cells.tsv"));cogs<-fread(file.path(sd,"cogs.tsv"))$COG_category;mat<-read_sparse(file.path(sd,"cog_by_cell.mtx.gz"));rownames(mat)<-cogs;colnames(mat)<-meta$cell_id
 lab<-cells[species==sp,.(cell_id,sample,AMP_program)];keep<-intersect(lab$cell_id,colnames(mat));lab<-lab[match(keep,cell_id)];mat<-mat[,keep,drop=FALSE]
 totals<-Matrix::colSums(mat);norm<-log1p(sweep(as.matrix(mat),2,pmax(totals,1),"/")*1e4);det<-as.matrix(mat)>0
 long<-rbindlist(lapply(seq_len(nrow(norm)),function(i)data.table(cell_id=colnames(norm),COG_category=rownames(norm)[i],score=norm[i,],detected=as.integer(det[i,]))))
 long<-merge(long,lab,by="cell_id")
 long[,z_score:=if(.N>1&&sd(score)>0)as.numeric(scale(score)) else 0,by=.(sample,COG_category)]
 fs<-long[,.(mean_z=mean(z_score),detection_fraction=mean(detected),cells=.N),by=.(sample,AMP_program,COG_category)]
 fg<-fs[,.(standardized_COG_score=mean(mean_z),signal_fraction=mean(detection_fraction),sd_across_samples=sd(mean_z),samples=.N,total_cell_records=sum(cells),pvalue=if(.N>=3)tryCatch(t.test(mean_z,mu=0)$p.value,error=function(e)NA_real_)else NA_real_),by=.(AMP_program,COG_category)]
 fg[,species:=sp];fs[,species:=sp];f_all[[sp]]<-fg;f_sample_all[[sp]]<-fs
}
f_summary<-rbindlist(f_all);f_summary[,FDR:=p.adjust(pvalue,"BH"),by=species]
fwrite(rbindlist(f_sample_all),file.path(out,"figure5F_within_sample_COG_scores.tsv"),sep="\t");fwrite(f_summary,file.path(out,"figure5F_COG_program_profiles.tsv"),sep="\t")
fwrite(data.table(parameter=c("random_seed","minimum_major_amp_program_cells","stratified_permutations"),value=c(cfg$embedding$random_seed,min_cells,B)),file.path(out,"analysis_parameters.tsv"),sep="\t")
cat("cells",nrow(cells),"major programs",length(major),"permutations",B,"\n");print(metrics)
