#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 31_candidate_amp_function_association.R CONFIG.json")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("R package jsonlite is required")
cfg <- jsonlite::fromJSON(args[[1]])
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE)
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}
clean_cog <- function(x) gsub("[^A-Z]", "", ifelse(is.na(x), "", x))
has_cog <- function(x, code) grepl(code, clean_cog(x), fixed = TRUE)

message("Reading cell and AMP data")
cells <- read_tsv(cfg$cell_metadata)
amp_long <- read_tsv(cfg$amp_counts)
cells <- cells[!duplicated(cells$cell_id), ]
amp_long <- amp_long[amp_long$cell_id %in% cells$cell_id & amp_long$umi_count > 0, ]

species_dirs <- list.dirs(cfg$species_matrix_root, recursive = FALSE, full.names = TRUE)
if (!length(species_dirs)) stop("No species matrix directories found")

module_labels <- c(
  defense = "Defense mechanisms",
  envelope = "Cell wall/membrane/envelope biogenesis",
  secretion = "Secretion/trafficking",
  signaling = "Signal transduction",
  secondary = "Secondary metabolite biosynthesis",
  bacteriocin = "Bacteriocin transport/immunity",
  t6ss = "Type VI secretion system",
  quorum_biofilm = "Quorum sensing/biofilm",
  stress = "Stress response"
)

module_gene_rows <- list()
module_score_rows <- list()
filter_rows <- list()
association_rows <- list()

for (sp_dir in species_dirs) {
  sp_cells <- read_tsv(file.path(sp_dir, "cells.tsv"))
  sp_cells <- merge(sp_cells, cells[, c("cell_id", "species")], by = "cell_id", all.x = TRUE)
  species_name <- unique(na.omit(sp_cells$species))
  if (length(species_name) != 1L) {
    warning("Skipping ambiguous species directory: ", sp_dir)
    next
  }
  species_name <- species_name[[1]]
  genes <- read_tsv(file.path(sp_dir, "genes.tsv"))
  mtx_path <- file.path(sp_dir, "gene_by_cell.mtx.gz")
  con <- gzfile(mtx_path, open = "rt")
  counts <- Matrix::readMM(con)
  close(con)
  counts <- as(counts, "dgCMatrix")
  if (nrow(counts) != nrow(genes) || ncol(counts) != nrow(sp_cells)) {
    stop("Matrix dimension mismatch in ", sp_dir)
  }

  text <- tolower(paste(
    ifelse(is.na(genes$product), "", genes$product),
    ifelse(is.na(genes$COG_description), "", genes$COG_description)
  ))
  module_index <- list(
    defense = which(has_cog(genes$COG_category, "V")),
    envelope = which(has_cog(genes$COG_category, "M")),
    secretion = which(has_cog(genes$COG_category, "U")),
    signaling = which(has_cog(genes$COG_category, "T")),
    secondary = which(has_cog(genes$COG_category, "Q")),
    bacteriocin = which(grepl(
      "bacteriocin|lantibiotic|microcin|colicin|competence.*bacteriocin|bacteriocin.*(transport|immunity)|toxin[- ]antitoxin|toxin.*immunity",
      text, perl = TRUE
    )),
    t6ss = which(grepl(
      "type vi secretion|type 6 secretion|\\bt6ss\\b|\\bhcp\\b|\\bvgrg\\b|\\bpaar\\b|\\btss[a-m][0-9]*\\b",
      text, perl = TRUE
    )),
    quorum_biofilm = which(grepl(
      "quorum|autoinducer|biofilm|diguanylate cyclase|c-di-gmp|cyclic-di-gmp|\\bggdef\\b|\\beal domain\\b|\\bpilz\\b|\\blux[irqsp]\\b",
      text, perl = TRUE
    )),
    stress = which(grepl(
      "stress|heat shock|cold shock|oxidative|peroxide|superoxide|catalase|thioredoxin|glutaredoxin|chaperon|osmotic|universal stress|dna damage|sos response",
      text, perl = TRUE
    ))
  )
  module_index <- lapply(module_index, unique)
  for (mn in names(module_index)) {
    idx <- module_index[[mn]]
    if (length(idx)) {
      module_gene_rows[[length(module_gene_rows) + 1L]] <- data.frame(
        species = species_name, module_id = mn, module = unname(module_labels[mn]),
        gene_id = genes$gene_id[idx], product = genes$product[idx],
        COG_category = genes$COG_category[idx], COG_description = genes$COG_description[idx],
        stringsAsFactors = FALSE
      )
    }
  }

  total_gene_umi <- Matrix::colSums(counts)
  score_df <- data.frame(
    cell_id = sp_cells$cell_id, sample = sp_cells$sample, species = species_name,
    total_gene_umi = total_gene_umi, stringsAsFactors = FALSE
  )
  for (mn in names(module_labels)) {
    idx <- module_index[[mn]]
    if (length(idx) >= cfg$min_module_genes) {
      module_umi <- Matrix::colSums(counts[idx, , drop = FALSE])
      score_df[[mn]] <- log1p(1e4 * module_umi / pmax(total_gene_umi, 1))
    } else {
      score_df[[mn]] <- NA_real_
    }
  }
  module_score_rows[[length(module_score_rows) + 1L]] <- score_df

  sp_amp <- amp_long[amp_long$cell_id %in% score_df$cell_id, ]
  if (!nrow(sp_amp)) next
  amp_ids <- unique(sp_amp$amp_id)
  for (amp_id in amp_ids) {
    ad <- sp_amp[sp_amp$amp_id == amp_id, c("cell_id", "umi_count")]
    dat <- merge(score_df, ad, by = "cell_id", all.x = TRUE)
    dat$umi_count[is.na(dat$umi_count)] <- 0
    dat$amp_expr <- log1p(1e4 * dat$umi_count / pmax(dat$total_gene_umi, 1))
    n_pos <- sum(dat$umi_count > 0)
    prevalence <- n_pos / nrow(dat)
    pos_samples <- length(unique(dat$sample[dat$umi_count > 0]))
    amp_sd <- sd(dat$amp_expr)
    reasons <- character()
    if (n_pos < cfg$min_positive_cells) reasons <- c(reasons, "too_few_positive_cells")
    if (pos_samples < cfg$min_positive_samples) reasons <- c(reasons, "too_few_positive_samples")
    if (prevalence < cfg$min_prevalence) reasons <- c(reasons, "prevalence_below_minimum")
    if (prevalence > cfg$max_prevalence) reasons <- c(reasons, "near_ubiquitous")
    if (!is.finite(amp_sd) || amp_sd < cfg$min_amp_expression_sd) reasons <- c(reasons, "insufficient_expression_variation")
    passed <- !length(reasons)
    filter_rows[[length(filter_rows) + 1L]] <- data.frame(
      species = species_name, amp_id = amp_id, n_cells = nrow(dat),
      positive_cells = n_pos, positive_samples = pos_samples,
      prevalence = prevalence, amp_expression_sd = amp_sd, passed = passed,
      exclusion_reason = if (passed) "" else paste(reasons, collapse = ";"),
      stringsAsFactors = FALSE
    )
    if (!passed) next

    dat$amp_z <- safe_z(dat$amp_expr)
    dat$depth_z <- safe_z(log1p(dat$total_gene_umi))
    for (mn in names(module_labels)) {
      y <- dat[[mn]]
      usable <- is.finite(y) & is.finite(dat$amp_z) & is.finite(dat$depth_z)
      d <- dat[usable, , drop = FALSE]
      y <- y[usable]
      if (nrow(d) < 30L || sd(y) == 0 || sd(d$amp_z) == 0) next
      d$module_z <- safe_z(y)
      fit <- try(lm(module_z ~ amp_z + depth_z + factor(sample), data = d), silent = TRUE)
      if (inherits(fit, "try-error")) next
      co <- summary(fit)$coefficients
      if (!"amp_z" %in% rownames(co)) next
      association_rows[[length(association_rows) + 1L]] <- data.frame(
        species = species_name, amp_id = amp_id, module_id = mn,
        module = unname(module_labels[mn]), n_cells = nrow(d),
        positive_cells = n_pos, positive_samples = pos_samples,
        prevalence = prevalence, beta = co["amp_z", "Estimate"],
        standard_error = co["amp_z", "Std. Error"],
        t_value = co["amp_z", "t value"], p_value = co["amp_z", "Pr(>|t|)"],
        model = "scale(module_score) ~ scale(AMP_expression) + scale(log1p(total_gene_UMI)) + sample",
        stringsAsFactors = FALSE
      )
    }
  }
}

module_genes <- if (length(module_gene_rows)) do.call(rbind, module_gene_rows) else data.frame()
module_scores <- if (length(module_score_rows)) do.call(rbind, module_score_rows) else data.frame()
filters <- if (length(filter_rows)) do.call(rbind, filter_rows) else data.frame()
associations <- if (length(association_rows)) do.call(rbind, association_rows) else data.frame()
if (!nrow(associations)) stop("No AMP-module models were estimable")
associations$fdr <- ave(associations$p_value, associations$species, FUN = function(x) p.adjust(x, "BH"))
associations$reliable <- associations$fdr < cfg$association_fdr

message("Mapping AMP reference and prediction annotations")
ref_meta <- read_tsv(cfg$amp_reference_metadata)
membership <- read_tsv(cfg$amp_reference_membership)
novel_members <- membership[membership$source == "Novel", c("amp_id", "source_id")]
pred <- read_tsv(cfg$four_model_predictions)
pred_keep <- intersect(c(
  "sequence_id", "sequence", "length", "macrel_hemolytic",
  "macrel_hemolytic_probability", "esm2_probability", "support_count_4",
  "support_tier", "exact_pattern"
), names(pred))
pred <- pred[, pred_keep, drop = FALSE]
novel_pred <- merge(novel_members, pred, by.x = "source_id", by.y = "sequence_id", all.x = TRUE)
if (nrow(novel_pred)) {
  ord <- order(novel_pred$amp_id, -ifelse(is.na(novel_pred$esm2_probability), -Inf, novel_pred$esm2_probability))
  novel_pred <- novel_pred[ord, ]
  novel_pred <- novel_pred[!duplicated(novel_pred$amp_id), ]
}

annotation <- merge(ref_meta, novel_pred, by = "amp_id", all.x = TRUE)
annotation$novelty_class <- ifelse(
  annotation$n_APD6 > 0, "known-like",
  ifelse(annotation$n_AMPSphere > 0, "remote homolog", "potentially novel")
)
annotation$prediction_grade <- ifelse(
  !is.na(annotation$support_tier), annotation$support_tier,
  ifelse(annotation$n_APD6 > 0, "curated AMP reference",
         ifelse(annotation$n_AMPSphere > 0, "AMPSphere reference", "not available"))
)
annotation$toxicity_prediction <- "not available"
annotation$hemolysis_prediction <- ifelse(
  is.na(annotation$macrel_hemolytic), "not available", annotation$macrel_hemolytic
)
annotation$esm2_lightgbm_score <- annotation$esm2_probability

amp_summary <- aggregate(cbind(total_umi = amp_long$umi_count) ~ amp_id, amp_long, sum)
amp_pos <- aggregate(cell_id ~ amp_id, amp_long, function(x) length(unique(x)))
names(amp_pos)[2] <- "positive_cells_all_targets"
annotation <- merge(annotation, amp_summary, by = "amp_id", all.x = TRUE)
annotation <- merge(annotation, amp_pos, by = "amp_id", all.x = TRUE)

assoc_anno <- merge(associations, annotation, by = "amp_id", all.x = TRUE)
sig_summary <- aggregate(
  cbind(max_abs_beta = abs(assoc_anno$beta), min_fdr = assoc_anno$fdr) ~ species + amp_id,
  assoc_anno,
  FUN = function(x) c(max = max(x, na.rm = TRUE), min = min(x, na.rm = TRUE))
)
# Recompute explicitly to avoid matrix-valued aggregate columns across R versions.
candidate_stats <- do.call(rbind, lapply(split(assoc_anno, interaction(assoc_anno$species, assoc_anno$amp_id, drop = TRUE)), function(d) {
  data.frame(
    species = d$species[1], amp_id = d$amp_id[1],
    max_abs_beta = max(abs(d$beta), na.rm = TRUE),
    max_abs_beta_reliable = if (any(d$reliable)) max(abs(d$beta[d$reliable]), na.rm = TRUE) else NA_real_,
    min_fdr = min(d$fdr, na.rm = TRUE), n_reliable_modules = sum(d$reliable),
    positive_cells = d$positive_cells[1], positive_samples = d$positive_samples[1],
    prevalence = d$prevalence[1], stringsAsFactors = FALSE
  )
}))
candidate_table <- merge(candidate_stats, annotation, by = "amp_id", all.x = TRUE)
candidate_table$novelty_priority <- match(candidate_table$novelty_class,
                                           c("known-like", "remote homolog", "potentially novel")) - 1
candidate_table$safety_priority <- ifelse(
  grepl("^NonHemo$", candidate_table$hemolysis_prediction, ignore.case = TRUE), 1, 0
)
candidate_table$model_priority <- ifelse(is.na(candidate_table$esm2_lightgbm_score), 0,
                                          candidate_table$esm2_lightgbm_score)
candidate_table$ranking_score <-
  3 * as.numeric(candidate_table$n_reliable_modules > 0) +
  pmin(candidate_table$max_abs_beta, 1) +
  0.5 * candidate_table$model_priority +
  0.25 * candidate_table$novelty_priority +
  0.25 * candidate_table$safety_priority +
  0.1 * log1p(candidate_table$positive_cells)
candidate_table <- candidate_table[order(
  -as.numeric(candidate_table$n_reliable_modules > 0),
  -candidate_table$ranking_score, candidate_table$min_fdr
), ]

# Balanced first pass across species, followed by global ranking.
first <- unlist(lapply(split(seq_len(nrow(candidate_table)), candidate_table$species), function(ii) {
  head(ii, cfg$max_candidates_per_species_first_pass)
}), use.names = FALSE)
selected_idx <- unique(c(first, seq_len(nrow(candidate_table))))
selected_idx <- head(selected_idx, min(cfg$top_candidates, nrow(candidate_table)))
candidate_table$selected_for_plot <- FALSE
candidate_table$selected_for_plot[selected_idx] <- TRUE
candidate_table$recommend_for_validation <- FALSE
recommendable <- which(
  candidate_table$n_reliable_modules > 0 &
    candidate_table$hemolysis_prediction != "Hemo"
)
recommendable <- head(recommendable, 5L)
candidate_table$recommend_for_validation[recommendable] <- TRUE

write_tsv(module_genes, file.path(cfg$output_dir, "tables", "module_gene_membership.tsv"))
write_tsv(module_scores, file.path(cfg$output_dir, "tables", "cell_module_scores.tsv"))
write_tsv(filters, file.path(cfg$output_dir, "tables", "amp_expression_filter_audit.tsv"))
write_tsv(assoc_anno, file.path(cfg$output_dir, "tables", "amp_function_module_associations.tsv"))
write_tsv(candidate_table, file.path(cfg$output_dir, "tables", "candidate_amp_annotations_and_ranking.tsv"))
write_tsv(candidate_table[candidate_table$recommend_for_validation, ],
          file.path(cfg$output_dir, "tables", "recommended_candidates.tsv"))

summary_lines <- c(
  "# Candidate AMP–function module analysis",
  "",
  paste0("- Target species: ", length(unique(associations$species))),
  paste0("- AMP–species pairs tested: ", length(unique(paste(associations$species, associations$amp_id)))),
  paste0("- Association models: ", nrow(associations)),
  paste0("- Reliable associations (within-species BH FDR < ", cfg$association_fdr, "): ",
         sum(associations$reliable)),
  paste0("- Candidates selected for the figure: ", sum(candidate_table$selected_for_plot)),
  paste0("- Candidates recommended for experimental prioritization: ",
         sum(candidate_table$recommend_for_validation)),
  "",
  "Association coefficients are adjusted associations, not evidence of antimicrobial activity or causality.",
  "Toxicity is reported as unavailable because no traceable toxicity result was found in the supplied upstream table."
)
writeLines(summary_lines, file.path(cfg$output_dir, "analysis_summary.md"))
message("Association analysis complete: ", cfg$output_dir)
