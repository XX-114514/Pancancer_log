args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3L)
query_dir <- args[[1]]; cache <- args[[2]]; ncores <- as.integer(args[[3]])
suppressPackageStartupMessages({library(Matrix); library(SingleR)})
# ExperimentHub cache metadata identifies this pair as EH3492 / EH3493,
# Human Primary Cell Atlas logcounts / colData. Read directly, without downloads.
ref <- readRDS(file.path(cache, '1d748702c9eed_3508'))
labels <- readRDS(file.path(cache, '1d7483d3c0df_3509'))
stopifnot(ncol(ref) == nrow(labels), ncol(ref) == 713L,
          !anyDuplicated(rownames(ref)), !anyDuplicated(colnames(ref)))
# Cached colData has no rownames; its positional ordering is defined by the
# paired ExperimentHub resource, not independently recoverable sample IDs.
if (!is.null(rownames(labels))) stopifnot(identical(colnames(ref), rownames(labels)))
q <- as(readMM(file.path(query_dir, 'query.mtx')), 'CsparseMatrix')
rownames(q) <- readLines(file.path(query_dir, 'genes.txt'))
colnames(q) <- readLines(file.path(query_dir, 'cells.txt'))
stopifnot(!anyDuplicated(rownames(q)), !anyDuplicated(colnames(q)))
common <- intersect(rownames(q), rownames(ref)); stopifnot(length(common) >= 1000L)
out <- data.frame(cell_id=colnames(q))
for (level in c('main', 'fine')) {
    pred <- SingleR(test=q[common,,drop=FALSE], ref=ref[common,,drop=FALSE],
                    labels=labels[[paste0('label.', level)]],
                    de.method='classic', fine.tune=TRUE, prune=TRUE,
                    num.threads=ncores, BPPARAM=BiocParallel::SerialParam())
    stopifnot(identical(rownames(pred), colnames(q)))
    out[[paste0('hpca_', level, '_raw')]] <- pred$labels
    out[[paste0('hpca_', level, '_pruned')]] <- pred$pruned.labels
    out[[paste0('hpca_', level, '_delta_next')]] <- pred$delta.next
}
write.table(out, file.path(query_dir,'predictions.tsv'), sep='\t',
            row.names=FALSE, quote=FALSE, na='')
cat('SingleR',as.character(packageVersion('SingleR')),'cells',ncol(q),
    'common_genes',length(common),'\n')
