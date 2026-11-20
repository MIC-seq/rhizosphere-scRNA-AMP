#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Matrix); library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else '/public3/labmember/qianqh/soil_amp_analysis/14_metabolic_functions'
cfg <- fromJSON(file.path(root,'config/analysis.json'))
set.seed(cfg$seed)
tab <- function(name, d) fwrite(as.data.table(d),file.path(root,'tables',paste0(name,'.tsv')),sep='\t',na='NA')
rd <- function(path) fread(path,na.strings=c('NA',''),data.table=FALSE)
cells.all <- rd(file.path(root,'inputs/cells.tsv'))
members <- rd(file.path(root,'inputs/pathway_gene_membership.tsv'))
catalog <- rd(file.path(root,'inputs/pathway_catalogue.tsv'))
amp <- rd(file.path(root,'inputs/candidate_counts.tsv.gz'))
candidates <- rd(file.path(root,'inputs/existing_candidate_panel.tsv'))
candidates$scope <- ifelse(candidates$species %in% cfg$candidate_analysis_species,'exploratory_candidate_panel','aggregate_pathway_context_only')
tab('candidate_analysis_scope',candidates)
capture.output(sessionInfo(),file=file.path(root,'audit/R_analysis_sessionInfo.txt'))
all.elig <- all.cover <- all.scores <- all.states <- all.effects <- all.meta <- all.cand.effects <- all.cand.meta <- list()
all.member <- all.qc <- all.loso <- list()
objects <- list()

partial_r <- function(x,y,depth,ngene) {
  ok <- is.finite(x)&is.finite(y)&is.finite(depth)&is.finite(ngene)
  x<-x[ok];y<-y[ok];depth<-depth[ok];ngene<-ngene[ok]
  if(length(x)<10L || length(unique(x))<2L || length(unique(y))<2L) return(NA_real_)
  Z<-cbind(1,rank(log1p(depth),ties.method='average'),rank(log1p(ngene),ties.method='average'))
  rx<-lm.fit(Z,rank(x,ties.method='average'))$residuals
  ry<-lm.fit(Z,rank(y,ties.method='average'))$residuals
  if(sd(rx)<1e-8 || sd(ry)<1e-8) return(NA_real_)
  max(-0.999999,min(0.999999,cor(rx,ry)))
}

summarise_effect <- function(d,sp,pid,label,predictor,variant) {
  d<-d[is.finite(d$partial_r),,drop=FALSE]
  k<-nrow(d)
  if(!k) return(NULL)
  z<-atanh(d$partial_r);mz<-mean(z)
  se<-if(k>1) sd(z)/sqrt(k) else NA_real_
  qualified <- k>=cfg$min_samples_inference && is.finite(se) && se>1e-12
  p<-if(qualified) 2*pt(-abs(mz/se),df=k-1) else NA_real_
  ci<-if(qualified) tanh(mz+c(-1,1)*qt(.975,k-1)*se) else c(NA_real_,NA_real_)
  flips<-as.matrix(expand.grid(rep(list(c(-1,1)),k)))
  pflip<-mean(abs(as.vector(flips %*% z)/k)>=abs(mz)-1e-12)
  data.frame(species=sp,pathway_id=pid,pathway=label,predictor=predictor,variant=variant,
             n_samples=k,n_cells=sum(d$n_cells),mean_partial_r=tanh(mz),ci_low=ci[1],ci_high=ci[2],
             p_value=p,signflip_p=pflip,same_direction_samples=sum(sign(d$partial_r)==sign(mz)),
             min_sample_r=min(d$partial_r),max_sample_r=max(d$partial_r),
             eligible_inference=qualified,samples=paste(d$sample,collapse=';'),stringsAsFactors=FALSE)
}

for(sp in cfg$species) {
  message('Analysing ',sp)
  slug<-gsub('[^A-Za-z0-9_.-]+','_',sp)
  dir<-file.path(root,'inputs/matrices',slug)
  G<-rd(file.path(dir,'genes.tsv'));C<-rd(file.path(dir,'cells.tsv'))
  con<-gzfile(file.path(dir,'gene_by_cell.mtx.gz'),'rt');M<-as.matrix(readMM(con));close(con)
  stopifnot(nrow(M)==nrow(G),ncol(M)==nrow(C),all(colSums(M)==C$retained_gene_umi))
  stopifnot(all(colSums(M>0)==C$retained_detected_genes),!anyDuplicated(C$cell_id))
  depth<-colSums(M);ngene<-colSums(M>0)
  legacy<-as.logical(G$legacy_retained)
  legacy.depth<-colSums(M[legacy,,drop=FALSE]);legacy.ngene<-colSums(M[legacy,,drop=FALSE]>0)
  xburden<-log1p(as.numeric(C$AMP_UMI_Fig4))
  observed<-rowSums(M>0);totals<-rowSums(M)
  all.qc[[sp]]<-data.frame(species=sp,cells=nrow(C),samples=length(unique(C$sample)),
                          median_gene_UMI=median(depth),median_detected_genes=median(ngene),
                          q25_detected_genes=unname(quantile(ngene,.25)),q75_detected_genes=unname(quantile(ngene,.75)),
                          cells_at_least_20_genes=sum(ngene>=20),cells_at_least_50_genes=sum(ngene>=50),
                          zero_retained_UMI=sum(depth==0),restored_UMI=sum(depth-legacy.depth))
  P<-matrix(0,nrow(C),nrow(catalog),dimnames=list(C$cell_id,catalog$pathway_id))
  PS<-P;path.indices<-list();path.dominant<-list()
  for(j in seq_len(nrow(catalog))) {
    pid<-catalog$pathway_id[j];label<-catalog$pathway[j]
    idx<-which(G$gene_id %in% members$gene_id[members$species==sp & members$pathway_id==pid])
    path.indices[[pid]]<-idx
    n.supported<-sum(observed[idx]>=cfg$min_gene_detected_cells)
    eligible<-n.supported>=cfg$min_pathway_genes
    pc<-if(length(idx)) colSums(M[idx,,drop=FALSE]) else rep(0,nrow(C))
    score<-log1p(10000*pc/pmax(depth,1));P[,j]<-pc;PS[,j]<-score
    dominant<-if(length(idx)) idx[which.max(totals[idx])] else integer()
    path.dominant[[pid]]<-dominant
    dfrac<-if(sum(pc)>0) max(totals[idx])/sum(pc) else NA_real_
    all.elig[[length(all.elig)+1]]<-data.frame(species=sp,pathway_id=pid,pathway=label,group=catalog$group[j],
       annotated_observed_genes=length(idx),genes_detected_in_at_least_3_cells=n.supported,
       detected_cells=sum(pc>0),cell_prevalence=mean(pc>0),pathway_UMI=sum(pc),
       dominant_gene_fraction=dfrac,eligible_gene_set=eligible,
       exclusion_reason=if(eligible)'PASS' else 'fewer_than_5_members_detected_in_3_cells')
    all.scores[[length(all.scores)+1]]<-data.frame(cell_id=C$cell_id,species=sp,sample=C$sample,crop=C$crop,
       functional_state=C$functional_state,pathway_id=pid,pathway_UMI=pc,retained_gene_UMI=depth,
       detected_pathway_genes=if(length(idx)) colSums(M[idx,,drop=FALSE]>0) else 0,score=score,
       eligible_gene_set=eligible)
    for(sm in unique(C$sample)) {
      take<-which(C$sample==sm)
      n.dgene<-if(length(idx)) sum(rowSums(M[idx,take,drop=FALSE])>0) else 0
      all.cover[[length(all.cover)+1]]<-data.frame(species=sp,sample=sm,crop=C$crop[take[1]],pathway_id=pid,
         pathway=label,n_cells=length(take),positive_cells=sum(pc[take]>0),detected_pathway_genes=n.dgene,
         pathway_UMI=sum(pc[take]),total_gene_UMI=sum(depth[take]),
         normalized_pathway_UMI=10000*sum(pc[take])/max(sum(depth[take]),1))
      for(st in unique(C$functional_state)) {
        si<-take[C$functional_state[take]==st]
        if(!length(si)) next
        all.states[[length(all.states)+1]]<-data.frame(species=sp,sample=sm,crop=C$crop[si[1]],functional_state=st,
           pathway_id=pid,pathway=label,group=catalog$group[j],n_cells=length(si),pathway_UMI=sum(pc[si]),
           total_gene_UMI=sum(depth[si]),detected_cells=sum(pc[si]>0),cell_prevalence=mean(pc[si]>0),
           detected_pathway_genes=if(length(idx)) sum(rowSums(M[idx,si,drop=FALSE])>0) else 0,
           pseudobulk_score=log1p(10000*sum(pc[si])/max(sum(depth[si]),1)),
           eligible_display=eligible && length(si)>=cfg$min_state_sample_cells)
      }
    }
    if(!eligible) next
    variants<-list(primary=list(idx=idx,depth=depth,ngene=ngene),
                   legacy_filter=list(idx=idx[legacy[idx]],depth=legacy.depth,ngene=legacy.ngene),
                   remove_dominant=list(idx=setdiff(idx,dominant),depth=depth,ngene=ngene))
    for(v in names(variants)) {
      vi<-variants[[v]]$idx;vd<-variants[[v]]$depth;vn<-variants[[v]]$ngene
      vp<-if(length(vi)) colSums(M[vi,,drop=FALSE]) else rep(0,nrow(C))
      vs<-log1p(10000*vp/pmax(vd,1));effects<-list()
      for(sm in unique(C$sample)) {
        ii<-which(C$sample==sm & vd>0)
        dg<-if(length(vi)) sum(rowSums(M[vi,ii,drop=FALSE])>0) else 0
        ok<-length(ii)>=cfg$min_sample_cells && sum(vp[ii]>0)>=cfg$min_positive_pathway_cells_per_sample && dg>=cfg$min_detected_pathway_genes_per_sample
        r<-if(ok) partial_r(xburden[ii],vs[ii],vd[ii],vn[ii]) else NA_real_
        effects[[length(effects)+1]]<-data.frame(species=sp,pathway_id=pid,pathway=label,predictor='total_cAMP_UMI',variant=v,
          sample=sm,n_cells=length(ii),positive_pathway_cells=sum(vp[ii]>0),detected_pathway_genes=dg,
          partial_r=r,eligible=ok,reason=if(ok)'PASS' else 'insufficient_sample_cells_or_pathway_detection')
      }
      ed<-rbindlist(effects);all.effects[[length(all.effects)+1]]<-ed
      meta<-summarise_effect(ed,sp,pid,label,'total_cAMP_UMI',v)
      if(!is.null(meta)) all.meta[[length(all.meta)+1]]<-meta
      if(v=='primary') {
        ef<-as.data.frame(ed[is.finite(partial_r)])
        if(nrow(ef)>=3) for(sm in ef$sample) {
          oo<-ef[ef$sample!=sm,,drop=FALSE]
          all.loso[[length(all.loso)+1]]<-data.frame(species=sp,pathway_id=pid,omitted_sample=sm,
             n_remaining_samples=nrow(oo),mean_partial_r=tanh(mean(atanh(oo$partial_r))))
        }
      }
    }
    # Candidate detail is limited to the fixed environmental-species panel.
    if(sp %in% cfg$candidate_analysis_species) {
      for(aid in candidates$amp_id[candidates$species==sp]) {
        ar<-amp[amp$amp_id==aid,,drop=FALSE]
        ax<-ar$umi_count[match(C$cell_id,ar$cell_id)];ax[is.na(ax)]<-0
        effects<-list()
        for(sm in unique(C$sample)) {
          ii<-which(C$sample==sm & depth>0)
          dg<-sum(rowSums(M[idx,ii,drop=FALSE])>0)
          ok<-length(ii)>=cfg$min_sample_cells && sum(pc[ii]>0)>=cfg$min_positive_pathway_cells_per_sample &&
             dg>=cfg$min_detected_pathway_genes_per_sample && sum(ax[ii]>0)>=cfg$min_candidate_positive_cells_per_sample && sum(ax[ii]==0)>=5
          rr<-if(ok) partial_r(log1p(ax[ii]),score[ii],depth[ii],ngene[ii]) else NA_real_
          effects[[length(effects)+1]]<-data.frame(species=sp,pathway_id=pid,pathway=label,predictor=aid,variant='primary',
              sample=sm,n_cells=length(ii),positive_candidate_cells=sum(ax[ii]>0),partial_r=rr,eligible=ok)
        }
        ed<-rbindlist(effects);all.cand.effects[[length(all.cand.effects)+1]]<-ed
        meta<-summarise_effect(ed,sp,pid,label,aid,'primary')
        if(!is.null(meta)) all.cand.meta[[length(all.cand.meta)+1]]<-meta
      }
      # Member-level output is descriptive; no gene-wise significance or target prioritization.
      for(g in idx) for(sm in unique(C$sample)) {
        ii<-which(C$sample==sm)
        all.member[[length(all.member)+1]]<-data.frame(species=sp,pathway_id=pid,pathway=label,sample=sm,
          gene_id=G$gene_id[g],gene_name=G$Preferred_name[g],product=G$product[g],
          n_cells=length(ii),gene_UMI=sum(M[g,ii]),detected_cells=sum(M[g,ii]>0),
          pathway_UMI=sum(pc[ii]),pathway_UMI_fraction=if(sum(pc[ii])>0)sum(M[g,ii])/sum(pc[ii]) else NA_real_)
      }
    }
  }
  objects[[sp]]<-list(cells=C,genes=G,counts=M,pathway_counts=P,pathway_scores=PS,
                      pathway_indices=path.indices,dominant=path.dominant)
}

tab('cell_depth_summary',rbindlist(all.qc))
tab('pathway_eligibility',rbindlist(all.elig))
tab('sample_pathway_coverage',rbindlist(all.cover))
fwrite(rbindlist(all.scores),file.path(root,'tables/cell_pathway_scores.tsv.gz'),sep='\t')
state<-rbindlist(all.states);tab('sample_state_pathway_scores',state)
state.mean<-state[eligible_display==TRUE,.(mean_pseudobulk_score=mean(pseudobulk_score),
  sd_pseudobulk_score=sd(pseudobulk_score),mean_cell_prevalence=mean(cell_prevalence),
  n_samples=.N,n_cells=sum(n_cells),min_detected_genes=min(detected_pathway_genes)),
  by=.(species,functional_state,pathway_id,pathway,group)]
state.mean[,within_species_row_z:=if(.N>1 && sd(mean_pseudobulk_score)>0) as.numeric(scale(mean_pseudobulk_score)) else 0,
           by=.(species,pathway_id)]
tab('state_pathway_sample_balanced',state.mean)
tab('sample_burden_pathway_effects',rbindlist(all.effects))
met<-rbindlist(all.meta,fill=TRUE)
met[,fdr:=p.adjust(p_value,'BH'),by=variant]
met[,signflip_fdr:=p.adjust(ifelse(eligible_inference,signflip_p,NA_real_),'BH'),by=variant]
meta.primary<-met[variant=='primary']
loso<-rbindlist(all.loso);tab('burden_pathway_leave_one_sample_out',loso)
loso.summary<-loso[,.(loso_min_r=min(mean_partial_r),loso_max_r=max(mean_partial_r),
                     loso_direction_consistent=all(mean_partial_r>0)|all(mean_partial_r<0)),by=.(species,pathway_id)]
meta.primary<-merge(meta.primary,loso.summary,by=c('species','pathway_id'),all.x=TRUE)
for(v in c('legacy_filter','remove_dominant')) {
  sub<-met[variant==v,.(species,pathway_id,mean_partial_r,n_samples)]
  setnames(sub,c('mean_partial_r','n_samples'),paste0(v,c('_r','_n_samples')))
  meta.primary<-merge(meta.primary,sub,by=c('species','pathway_id'),all.x=TRUE)
}
tab('burden_pathway_meta_all_variants',met)
tab('burden_pathway_primary',meta.primary)
tab('sample_candidate_pathway_effects',rbindlist(all.cand.effects))
cmet<-rbindlist(all.cand.meta,fill=TRUE)
if(nrow(cmet)) {
  cmet[,fdr:=p.adjust(p_value,'BH')]
  cmet[,signflip_fdr:=p.adjust(ifelse(eligible_inference,signflip_p,NA_real_),'BH')]
}
tab('candidate_pathway_exploratory',cmet)
tab('environmental_pathway_gene_contributions',rbindlist(all.member))
saveRDS(objects,file.path(root,'inputs/analysis_objects.rds'))
stats<-data.frame(metric=c('cells_retained','species','independent_samples','curated_pathways',
 'eligible_species_pathway_sets','primary_tests','primary_FDR_below_0.05','primary_signflip_FDR_below_0.05',
 'candidate_tests','candidate_FDR_below_0.05'),
 value=c(nrow(cells.all),length(cfg$species),length(unique(cells.all$sample)),nrow(catalog),
 sum(rbindlist(all.elig)$eligible_gene_set),sum(is.finite(meta.primary$p_value)),
 sum(meta.primary$fdr<.05,na.rm=TRUE),sum(meta.primary$signflip_fdr<.05,na.rm=TRUE),
 sum(is.finite(cmet$p_value)),sum(cmet$fdr<.05,na.rm=TRUE)))
tab('analysis_summary',stats)
print(stats)
print(meta.primary[order(fdr)][1:min(12,.N),.(species,pathway,mean_partial_r,n_samples,p_value,fdr,signflip_p,legacy_filter_r,remove_dominant_r)])
