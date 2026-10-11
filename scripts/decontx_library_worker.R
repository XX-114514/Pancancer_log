#!/usr/bin/env Rscript
# Actual library counts only. No global pools, fixed contaminant-gene lists,
# or forced contamination priors. Original inputs are never modified.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
directory <- args[[1]]; seed <- as.integer(args[[2]])
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(decontX))
suppressPackageStartupMessages(library(jsonlite))
mat <- as(readMM(file.path(directory, "input_counts_private.mtx")), "CsparseMatrix")
cells <- read.delim(file.path(directory, "cells_private.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
genes <- readLines(file.path(directory, "genes_private.tsv"))
stopifnot(ncol(mat)==nrow(cells), nrow(mat)==length(genes), !anyDuplicated(cells$stable_cell_id))
colnames(mat) <- cells$stable_cell_id; rownames(mat) <- genes
z <- as.integer(factor(cells$broad_group))
stopifnot(length(unique(z))>=2)
started <- format(Sys.time(), tz="UTC", usetz=TRUE)
result <- decontX(mat, z=z, maxIter=500, delta=c(10,10), estimateDelta=TRUE,
                  convergence=0.001, iterLogLik=10, varGenes=5000, seed=seed,
                  logfile=file.path(directory,"decontX_method.log"),verbose=TRUE)
stopifnot(identical(dim(result$decontXcounts), dim(mat)), all(is.finite(result$contamination)))
writeMM(result$decontXcounts,file.path(directory,"corrected_counts_private.mtx"))
write.table(data.frame(stable_cell_id=cells$stable_cell_id, contamination=result$contamination,
                      cluster=result$z),file.path(directory,"decontX_cells_private.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
saveRDS(result$estimates,file.path(directory,"decontX_estimates_private.rds"))
method_log <- readLines(file.path(directory,"decontX_method.log"),warn=FALSE)
iteration_lines <- grep("Completed iteration:.*converge:",method_log,value=TRUE)
final_line <- tail(iteration_lines,1)
convergence_metric <- if(length(final_line)) as.numeric(sub(".*converge: ([0-9.eE+-]+).*","\\1",final_line)) else NA_real_
iterations <- if(length(final_line)) as.integer(sub(".*Completed iteration: ([0-9]+).*","\\1",final_line)) else NA_integer_
write_json(list(status="method_completed",start_time=started,end_time=format(Sys.time(),tz="UTC",usetz=TRUE),
    decontX_version=as.character(packageVersion("decontX")),R_version=R.version.string,
    n_cells=ncol(mat),n_genes=nrow(mat),n_clusters=length(unique(result$z)),
    background="no_empty_droplets_filtered_only_estimate",
    convergence_metric=convergence_metric,iterations=iterations,
    converged=is.finite(convergence_metric) && convergence_metric<0.001,
    median_contamination=median(result$contamination), max_contamination=max(result$contamination),
    params=list(seed=seed,maxIter=500,estimateDelta=TRUE,delta=c(10,10),convergence=0.001)),
    file.path(directory,"decontX_summary.json"),pretty=TRUE,auto_unbox=TRUE)
