#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(jsonlite)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("config matrix_root output_root")
cfg <- fromJSON(args[1], simplifyVector=TRUE)
matrix_root <- args[2]
output_root <- args[3]
dir.create(output_root, recursive=TRUE, showWarnings=FALSE)

cog_names <- c(
  A="RNA processing", B="Chromatin structure", C="Energy production",
  D="Cell cycle", E="Amino acid metabolism", F="Nucleotide metabolism",
  G="Carbohydrate metabolism", H="Coenzyme metabolism", I="Lipid metabolism",
  J="Translation", K="Transcription", L="Replication and repair",
  M="Cell envelope biogenesis", N="Cell motility", O="Protein turnover",
  P="Inorganic ion transport", Q="Secondary metabolism", R="General prediction",
  S="Unknown function", T="Signal transduction", U="Intracellular trafficking",
  V="Defense mechanisms", W="Extracellular structures", X="Mobilome",
  Y="Nuclear structure", Z="Cytoskeleton"
)

slug <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
read_sparse <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  as(readMM(con), "dgCMatrix")
}

run_embedding <- function(species, representation, input_dir, output_dir, pars) {
  feature_file <- if (representation == "gene") "genes.tsv" else "cogs.tsv"
  matrix_file <- if (representation == "gene") "gene_by_cell.mtx.gz" else "cog_by_cell.mtx.gz"
  features <- fread(file.path(input_dir, feature_file))
  metadata <- fread(file.path(input_dir, "cells.tsv"))
  mat <- read_sparse(file.path(input_dir, matrix_file))
  feature_ids <- if (representation == "gene") features$gene_id else features$COG_category
  rownames(mat) <- make.unique(feature_ids)
  colnames(mat) <- metadata$cell_id
  metadata <- as.data.frame(metadata); rownames(metadata) <- metadata$cell_id
  minimum <- if (representation == "gene") pars$min_gene_detection_cells else pars$min_cog_detection_cells
  keep_features <- Matrix::rowSums(mat > 0) >= minimum
  mat <- mat[keep_features, , drop=FALSE]
  keep_cells <- Matrix::colSums(mat) > 0
  mat <- mat[, keep_cells, drop=FALSE]
  metadata <- metadata[colnames(mat), , drop=FALSE]
  if (nrow(mat) < 4 || ncol(mat) < 50) stop(species, " ", representation, " matrix is too sparse for UMAP")

  object <- CreateSeuratObject(counts=mat, meta.data=metadata, project=slug(species), min.cells=0, min.features=0)
  object <- NormalizeData(object, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
  if (representation == "gene") {
    nvar <- min(as.integer(pars$gene_variable_features), nrow(object))
    object <- FindVariableFeatures(object, selection.method="vst", nfeatures=nvar, verbose=FALSE)
    use_features <- VariableFeatures(object)
  } else {
    use_features <- rownames(object)
    VariableFeatures(object) <- use_features
  }
  object <- ScaleData(object, features=use_features, verbose=FALSE)
  npcs <- min(as.integer(pars$max_pca_dimensions), length(use_features)-1, ncol(object)-1)
  object <- RunPCA(object, features=use_features, npcs=npcs, seed.use=as.integer(pars$random_seed), verbose=FALSE)
  dims <- seq_len(min(npcs, as.integer(pars$max_pca_dimensions)))
  object <- RunUMAP(object, reduction="pca", dims=dims,
                    n.neighbors=min(as.integer(pars$umap_neighbors), ncol(object)-1),
                    min.dist=as.numeric(pars$umap_min_dist),
                    seed.use=as.integer(pars$random_seed), verbose=FALSE)
  object <- FindNeighbors(object, reduction="pca", dims=dims,
                          k.param=min(20, ncol(object)-1), verbose=FALSE)
  object <- FindClusters(object, resolution=as.numeric(pars$cluster_resolution),
                         random.seed=as.integer(pars$random_seed), verbose=FALSE)

  embedding <- as.data.table(Embeddings(object, "umap"), keep.rownames="cell_id")
  setnames(embedding, c("cell_id", "UMAP_1", "UMAP_2"))
  meta_out <- as.data.table(object@meta.data, keep.rownames="cell_id")
  keep_meta <- c("cell_id", "sample", "crop", "AMP_cluster", "AMP_UMI_Fig4",
                 "stable_AMP_count_Fig4", "gene_amp_load", "gene_total_non_rrna_umi",
                 "gene_detected_non_rrna_genes", "seurat_clusters")
  meta_out <- meta_out[, ..keep_meta]
  result <- merge(embedding, meta_out, by="cell_id", all.x=TRUE, sort=FALSE)
  result[, `:=`(species=species, representation=representation)]
  fwrite(result, file.path(output_dir, paste0(representation, "_umap_cells.tsv")), sep="\t")

  markers <- tryCatch(
    as.data.table(FindAllMarkers(object, only.pos=TRUE, min.pct=0.05,
                                 logfc.threshold=0.1, verbose=FALSE), keep.rownames="marker_row"),
    error=function(e) data.table(error=conditionMessage(e))
  )
  fwrite(markers, file.path(output_dir, paste0(representation, "_cluster_markers.tsv")), sep="\t")

  diagnostics <- data.table(
    species=species, representation=representation, input_cells=nrow(metadata),
    embedded_cells=ncol(object), retained_features=nrow(mat), variable_features=length(use_features),
    clusters=length(unique(Idents(object))), pca_dimensions=length(dims),
    cor_umap1_log_umi=cor(result$UMAP_1, log1p(as.numeric(result$gene_total_non_rrna_umi)), use="complete.obs"),
    cor_umap2_log_umi=cor(result$UMAP_2, log1p(as.numeric(result$gene_total_non_rrna_umi)), use="complete.obs")
  )

  if (representation == "cog") {
    avg <- AverageExpression(object, assays="RNA", group.by="seurat_clusters", slot="data", verbose=FALSE)$RNA
    state_rows <- rbindlist(lapply(colnames(avg), function(cluster) {
      values <- sort(avg[, cluster], decreasing=TRUE)
      top <- head(names(values), 3)
      data.table(cluster=cluster, rank=seq_along(top), COG_category=top,
                 COG_name=unname(cog_names[top]), mean_normalized_expression=unname(values[top]))
    }))
    fwrite(state_rows, file.path(output_dir, "cog_state_top_categories.tsv"), sep="\t")
  }
  diagnostics
}

all_diagnostics <- list()
for (species in cfg$species) {
  input_dir <- file.path(matrix_root, slug(species))
  output_dir <- file.path(output_root, slug(species)); dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  for (representation in c("gene", "cog")) {
    message("Running ", species, " ", representation)
    all_diagnostics[[length(all_diagnostics)+1]] <- run_embedding(
      species, representation, input_dir, output_dir, cfg$embedding)
  }
}
fwrite(rbindlist(all_diagnostics), file.path(output_root, "embedding_diagnostics.tsv"), sep="\t")
sessionInfo()
