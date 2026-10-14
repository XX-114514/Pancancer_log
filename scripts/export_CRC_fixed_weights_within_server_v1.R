# Serialize already frozen coefficients within the authorized server only.
args<-commandArgs(trailingOnly=TRUE)
root<-args[1];dest<-args[2]
stopifnot(!dir.exists(dest));dir.create(dest)
for(type in c('Tnode9','Macnode18')){
 for(program in c('MCP1','MCP2')){
  path<-file.path(root,'private_frozen_projection_v1',paste0(type,'_',program,'_fixed_weights.rds'))
  if(!file.exists(path))next
  x<-readRDS(path)
  stopifnot(length(x$genes)==length(x$weight),all(is.finite(x$weight)),all(x$gene_sd>0))
  write.table(data.frame(gene=x$genes,weight=x$weight,mean=x$gene_mean,sd=x$gene_sd),
    file.path(dest,paste0(type,'_',program,'_genes.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(term=names(x$quality_coefficients),value=unname(x$quality_coefficients)),
    file.path(dest,paste0(type,'_',program,'_quality.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 }
}
cat('FIXED_WEIGHTS_SERIALIZED_WITHIN_SERVER_ONLY\n')
