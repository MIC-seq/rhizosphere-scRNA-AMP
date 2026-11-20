raw_mes <- function(geneset_input = NA,count_matrix = NA,method = "ssGSEA"){
       if (method == "ssGSEA") {
           library(GSVA)
           library(GSEABase)
           gsva_es <- gsva(as.matrix(count_matrix), geneset_input,method=c("ssgsea"), kcdf=c("Poisson")) 
           gsva_es <- as.data.frame(gsva_es)
         }
       if (method == "GSVA") {
           library(GSVA)
           library(GSEABase)
           gsva_es <- gsva(as.matrix(count_matrix), geneset_input, method=c("gsva"), kcdf=c("Poisson"))
           gsva_es <- as.data.frame(gsva_es)
         }
       if (method == "AUCell") {
           library(AUCell)
           library(GSEABase)
           cells_rankings <- AUCell_buildRankings(as.matrix(count_matrix), plotStats=F)
           cells_AUC <- AUCell_calcAUC(geneset_input, cells_rankings)
           gsva_es <- data.frame(getAUC(cells_AUC))
         }
       if (method == "VISION") {
           library(VISION)
           n.umi <- colSums(count_matrix)
           scaled_counts <- t(t(count_matrix) / n.umi) * median(n.umi)
           vis <- Vision(scaled_counts, signatures = geneset_input)
           vis <- analyze(vis)
           gsva_es<-data.frame(t(vis@SigScores))
         }
       return(gsva_es)
     }

get_gene_count <- function(geneset_input = NA){
  gene_count <- matrix(data = NA,nrow = length(geneset_input),ncol = 2)
  for (i in 1:length(geneset_input)){
    length(c(geneset_input[[i]]@geneIds)) -> gene_count[i,1]
    gene_count[i,2] <- geneset_input[[i]]@setName
  }
  
  gene_count <- as.data.frame(gene_count)
  colnames(gene_count) <- c('gene_number','pathway_name')
  return(gene_count)
}

get_random_geneset <- function(pathway_gsva = NA,gene_count = NA,count_matrix = NA){
  total_pathway_random <- c()
  for (q in 1:nrow(pathway_gsva)){
    if (q %%10==0){print(q)}
    gene_select_number <- gene_count$gene_number[which(gene_count$pathway_name==rownames(pathway_gsva)[q])]
    gene_select_number <- as.numeric(gene_select_number)
    pathway_random <- c()
    for (i in 1:1000){
      sample(rownames(count_matrix),gene_select_number) -> gene_select
      pathway_random[[i]] <- c(gene_select)
      names(pathway_random)[[i]] <- paste('random',i,sep = "_")
    }
    pathway_random -> total_pathway_random[[rownames(pathway_gsva)[q]]]
  }
  return(total_pathway_random)
}


get_permutation <- function(pathway_gsva = NA,count_matrix = NA,total_pathway_random = NA){
  first_permutation <- c()
  for (q in 1:nrow(pathway_gsva)){
    if (q %%10==0){print(q)}
    gene_gsva <- gsva(as.matrix(count_matrix), total_pathway_random[[rownames(pathway_gsva)[q]]],method=c("ssgsea"), kcdf=c("Poisson"))
    first_permutation[[rownames(pathway_gsva)[q]]] <- gene_gsva
  }
  return(first_permutation)
}

author_normalize <- function(gsva_es, first_permutation) {
gsva_es_mean <- as.data.frame(apply(gsva_es,1,mean))
permutation_first_result <- matrix(data = NA,nrow = length(first_permutation),ncol = 1)
get_zcore <- matrix(data = NA,nrow = length(first_permutation),ncol = 1)
for (i in 1:length(first_permutation)){
  rank_data <- as.data.frame(first_permutation[[i]])
  apply(rank_data,1,mean) -> rank_data
  rank_data <- as.data.frame(rank_data)
  rank_data <- rbind(gsva_es_mean[i,1],rank_data)
  get_zcore_n <- scale(rank_data$rank_data)
  get_zcore[i,1] <- get_zcore_n[1,1]
  rank_test <- rank(-(rank_data[c(1:1001),]))
  for (q in 1:length(rank_test)){
    if(rank_test[q]>1){rank_test[q] <- rank_test[q]-1}
  }
  rank_test <- rank_test/1000
  rank_test[1] -> permutation_first_result[i,1]
}
permutation_first_result <- as.data.frame(permutation_first_result)
rownames(permutation_first_result) <- rownames(gsva_es_mean)
get_zcore <- as.data.frame(get_zcore)
rownames(get_zcore) <- rownames(gsva_es_mean)



return(list(MES=get_zcore, P=permutation_first_result))
}
