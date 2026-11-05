# One iNMF random start in its own process: RcppPlanc aborts the whole R process
# ("solve(): solution not found") on some starts, so run_r_method.R calls this per start.
# Usage: Rscript liger_start.R <run_dir> <K> <start_seed> <min_samples_per_dataset> <out.rds> <overlay_lib>
a <- commandArgs(TRUE); .libPaths(c(a[6], .libPaths())); suppressPackageStartupMessages(library(RcppPlanc))
X <- as.matrix(read.delim(file.path(a[1], "inputs/X_minmax.tsv"), row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(a[1], "inputs/metadata.tsv"), check.names = FALSE)
n_src <- table(meta$source); ds <- ifelse(n_src[meta$source] >= as.integer(a[4]), meta$source, "pooled_small_sources")
groups <- split(seq_len(ncol(X)), ds)
set.seed(as.integer(a[3]))
m <- inmf(lapply(groups, function(i) X[, i, drop = FALSE]), k = as.integer(a[2]), lambda = 5, niter = 30L, nCores = 1L)
saveRDS(list(m = m, groups = groups), a[5])
