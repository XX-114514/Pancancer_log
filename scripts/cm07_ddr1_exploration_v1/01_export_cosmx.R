# 导出 CosMx 全部基因(只读原对象)
suppressMessages({library(Seurat); library(Matrix)})
a <- commandArgs(TRUE); RUN <- a[1]; SRC <- a[2]
o <- readRDS(SRC); cnt <- GetAssayData(o, assay = "Nanostring", layer = "counts")
writeMM(cnt, file.path(RUN, "inputs/cosmx_counts.mtx"))
writeLines(rownames(cnt), file.path(RUN, "inputs/cosmx_genes.txt")); writeLines(colnames(cnt), file.path(RUN, "inputs/cosmx_cells.txt"))
cat(dim(cnt), "\n"); print(colnames(o@meta.data)); print(names(o@images)[1:3])
print(head(rownames(cnt)[grepl("^NegPrb|^Neg", rownames(cnt))]))
print(table(o$nb_clus, o$annotation)[1:6, ]); for (cn in c("niche","Niche","patient","Patient","sample","Sample")) if (cn %in% colnames(o@meta.data)) { cat(cn, "\n"); print(table(o@meta.data[[cn]])) }
