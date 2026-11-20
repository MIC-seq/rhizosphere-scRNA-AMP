#!/usr/bin/env Rscript
# Run from the extracted package root: Rscript scripts/01_render_figure1b.R
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script)==1L)
root <- dirname(dirname(normalizePath(script)))
for (pkg in c("ape", "svglite", "ragg")) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing R package: ", pkg)
}
if (!capabilities("cairo")) stop("R must support Cairo PDF output")
setwd(file.path(root, "supplemented207"))
source(file.path(root, "scripts", "plot_tree_reviewed.R"), chdir=FALSE)
