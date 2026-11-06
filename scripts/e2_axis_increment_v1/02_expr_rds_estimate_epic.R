# 生成 expr.rds, ESTIMATE 和 EPIC CAF。不拟合生存模型。
a <- commandArgs(TRUE); RUN <- a[1]; lib <- file.path(RUN, "private_R_lib"); .libPaths(c(lib, .libPaths())); IN <- file.path(RUN, "inputs")
E <- readRDS(file.path(IN, "expr.rds")); s <- colnames(E); stopifnot(!anyNA(E)); cat("expr.rds", dim(E), "\n")
library(parallel); library(tidyestimate)
f <- as.data.frame(E[rownames(E) %in% unlist(tidyestimate::common_genes$hgnc_symbol), ])
fc <- filter_common_genes(f, id = "hgnc_symbol", tidy = FALSE, tell_missing = TRUE, find_alias = FALSE)  # 第 1 列为基因符号
idx <- split(seq_len(ncol(E)), cut(seq_len(ncol(E)), 8, labels = FALSE))
es <- do.call(rbind, mclapply(idx, function(i) as.data.frame(estimate_score(fc[, c(1, i + 1), drop = FALSE], is_affymetrix = FALSE)), mc.cores = 8))
print(head(es, 2)); write.table(es, file.path(IN, "estimate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
library(EPIC)
idx2 <- split(seq_len(ncol(E)), cut(seq_len(ncol(E)), 8, labels = FALSE))
ep <- do.call(rbind, mclapply(idx2, function(i) {
  tp <- pmax(2^E[, i, drop = FALSE] - 0.001, 0)
  EPIC(bulk = tp, reference = "TRef", withOtherCells = TRUE)$cellFractions
}, mc.cores = 8))
ep <- ep[match(s, rownames(ep)), ]; print(colnames(ep))
write.table(data.frame(sample = s, caf_fraction = ep[, "CAFs"]), file.path(IN, "epic_caf.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
