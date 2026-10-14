#!/usr/bin/env Rscript
# Official pinned source, private library. Does not change shared R packages.
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
lib <- file.path(root,'library')
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
packages <- c('PMA_1.2-4.tar.gz','unikn_1.0.0.tar.gz','ppcor_1.1.tar.gz','UpSetR_1.4.1.tar.gz','lmerTest_3.1-3.tar.gz','DIALOGUE-9c146cc.tar.gz')
for (filename in packages) {
  install.packages(file.path(root,'sources',filename),repos=NULL,type='source',lib=lib)
}
required <- c('DIALOGUE','lme4','lmerTest','PMA','plyr','matrixStats','psych','stringi','RColorBrewer','unikn','reshape2','ggplot2','ppcor','Hmisc','grid','beanplot','UpSetR','assertthat')
available <- vapply(required,requireNamespace,logical(1),quietly=TRUE)
write.table(data.frame(package=required,available=available),file.path(root,'runtime_dependencies.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
if (!all(available))stop('Isolated dependency gate failed; previous installation evidence preserved')
library(DIALOGUE)
set.seed(20261014)
# Synthetic 24 samples, two cell types; no real expression data or clinical labels.
samples <- rep(sprintf('S%02d',1:24),each=30)
latent <- rep(scale(seq_len(24))[,1],each=30)
rA <- lapply(c('T','Myeloid'),function(name) {
  x <- matrix(rnorm(720*8),nrow=720,ncol=8)
  x[,1] <- latent+rnorm(720,sd=.2)
  cells <- paste0(name,'_',seq_len(720))
  rownames(x) <- cells
  gene <- matrix(rnorm(100*720),nrow=100,ncol=720)
  gene[1:12,] <- gene[1:12,]+rep(latent,each=12)
  rownames(gene) <- paste0('GENE',seq_len(100)); colnames(gene) <- cells
  make.cell.type(name=name,tpm=gene,samples=samples,X=x,metadata=data.frame(row.names=cells),cellQ=rnorm(720))
})
out <- file.path(root,'synthetic_smoke_v1')
dir.create(out,showWarnings=FALSE)
param <- DLG.get.param(k=2,results.dir=out,seed1=20261014,abn.c=20,n.genes=40,conf='cellQ',pheno=NULL,plot.flag=FALSE)
result <- DIALOGUE.run(rA,main='SyntheticAPICheck',param=param,plot.flag=FALSE)
saveRDS(result,file.path(out,'result.rds'))
writeLines(c(paste('R',getRversion()),paste('DIALOGUE',packageVersion('DIALOGUE')),paste('message',result$message),paste('real_cells',0)),file.path(root,'synthetic_smoke_summary.txt'))
cat('SYNTHETIC_API_FINISHED_REAL_CELLS_ZERO\n')
