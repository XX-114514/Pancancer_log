#!/usr/bin/env Rscript
# Per-unit malignancy evidence, two modes called in turn by malignancy_source_v1.py: "scatomic"
# (pan-cancer classifier, no CNV option) and then "infercnv"
# (reference = immune/stromal cells of the same unit; a held-out part of the reference is scored as
# "control" observations to calibrate the CNV threshold). Inputs written by malignancy_source_v1.py:
#   counts.mtx (genes x cells, raw), genes.txt, cells.tsv (cell, lineage, infercnv_group)
# infercnv_group: reference / control / observation / none. Outputs: scatomic.tsv, cnv_scores.tsv, status.json
args <- commandArgs(trailingOnly = TRUE)
d <- args[[1]]; gene_pos <- args[[2]]; ncores <- as.integer(args[[3]]); seed <- as.integer(args[[4]]); mode <- args[[5]]  # scatomic | infercnv
suppressPackageStartupMessages({library(Matrix); library(jsonlite)})
m <- as(readMM(file.path(d, "counts.mtx")), "CsparseMatrix")
rownames(m) <- readLines(file.path(d, "genes.txt"))
cells <- read.delim(file.path(d, "cells.tsv"), stringsAsFactors = FALSE)
colnames(m) <- cells$cell
stf <- file.path(d, "status.json")
st <- if (file.exists(stf)) fromJSON(stf) else list()

if (mode == "scatomic") {
# scATOMIC
t0 <- Sys.time()
r <- try({
  suppressPackageStartupMessages(library(scATOMIC))
  pred <- run_scATOMIC(rna_counts = m, mc.cores = ncores)
  res <- create_summary_matrix(prediction_list = pred, use_CNVs = FALSE, modify_results = TRUE,
                               mc.cores = ncores, raw_counts = m, min_prop = 0.5)
  res$cell <- rownames(res)
  write.table(res, file.path(d, "scatomic.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  TRUE
}, silent = TRUE)
st$scatomic <- if (isTRUE(r)) "done" else paste("failed:", substr(as.character(r), 1, 300))
st$scatomic_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

if (mode == "infercnv") {
# inferCNV
t0 <- Sys.time()
use <- cells$infercnv_group %in% c("reference", "control", "observation")
if (sum(cells$infercnv_group == "reference") < 50 || sum(cells$infercnv_group == "observation") < 1) {
  st$infercnv <- "not_evaluable_too_few_reference_or_observation"
} else {
  r <- try({
    suppressPackageStartupMessages(library(infercnv))
    ann <- data.frame(group = cells$infercnv_group[use], row.names = cells$cell[use])
    obj <- CreateInfercnvObject(raw_counts_matrix = m[, use], annotations_file = ann, delim = "\t",
                                gene_order_file = gene_pos, ref_group_names = "reference")
    obj <- infercnv::run(obj, cutoff = 0.1, out_dir = file.path(d, "infercnv"), cluster_by_groups = FALSE,
                         denoise = TRUE, HMM = FALSE, analysis_mode = "samples", no_plot = TRUE,
                         num_threads = ncores, write_expr_matrix = FALSE, save_rds = FALSE, save_final_rds = FALSE)
    x <- obj@expr.data
    score <- colMeans((x - 1)^2)
    # CNV correlation (Tirosh 2016 / Neftel 2019): Pearson correlation of each cell's CNV profile with the
    # mean profile of the top 5% scoring observation cells of the unit (clonal signal, sensitive where few
    # large events give a modest mean-squared score)
    obs_cells <- intersect(colnames(x), cells$cell[cells$infercnv_group == "observation"])
    top <- obs_cells[score[obs_cells] >= quantile(score[obs_cells], 0.95)]
    prof <- rowMeans(x[, top, drop = FALSE] - 1)
    corr <- as.numeric(cor(x - 1, prof))
    write.table(data.frame(cell = names(score), cnv_score = score, cnv_corr = corr), file.path(d, "cnv_scores.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    TRUE
  }, silent = TRUE)
  st$infercnv <- if (isTRUE(r)) "done" else paste("failed:", substr(as.character(r), 1, 300))
}
st$infercnv_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
unlink(file.path(d, "infercnv"), recursive = TRUE)
}
st$scATOMIC_version <- tryCatch(as.character(packageVersion("scATOMIC")), error = function(e) NA)
st$infercnv_version <- tryCatch(as.character(packageVersion("infercnv")), error = function(e) NA)
write_json(st, stf, auto_unbox = TRUE, pretty = TRUE)
