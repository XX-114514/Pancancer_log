args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L)
suppressPackageStartupMessages(library(SingleR))
refs<-read.delim(args[[1]],stringsAsFactors=FALSE,check.names=FALSE)
probes<-read.delim(args[[2]],stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(as.character(packageVersion('SingleR'))==args[[5]])
audit<-list()
for(i in seq_len(nrow(refs))) {
  r<-refs[i,];matrix<-readRDS(r$matrix_path);labels<-as.data.frame(readRDS(r$labels_path))
  stopifnot(identical(colnames(matrix),rownames(labels)))
  for(j in seq_len(nrow(probes))) {
    p<-probes[j,];query<-as.matrix(read.delim(gzfile(p$profile_path),row.names=1,check.names=FALSE))
    query<-query[,rev(seq_len(ncol(query))),drop=FALSE]
    original<-read.delim(file.path(args[[3]],paste0(p$source,'_',r$reference,'_private_SingleR_proposals_v1.tsv')),stringsAsFactors=FALSE,check.names=FALSE)
    for(level in c('main','fine')) {
      predicted<-SingleR(test=query,ref=matrix,labels=labels[[paste0('label.',level)]],de.method='classic',fine.tune=TRUE,prune=TRUE,num.threads=2L,BPPARAM=BiocParallel::SerialParam())
      old<-original[original$label_level==level,,drop=FALSE]
      old<-old[match(sub('^C','',rownames(predicted)),as.character(old$cluster)),,drop=FALSE]
      stopifnot(identical(as.character(old$raw_label),as.character(predicted$labels)),identical(as.character(old$pruned_label),as.character(predicted$pruned.labels)))
      stopifnot(max(abs(old$delta_next-predicted$delta.next))<1e-10)
      audit[[length(audit)+1L]]<-data.frame(source=p$source,reference=r$reference,level=level,n_profiles=nrow(old),labels_equal=TRUE,pruned_equal=TRUE,delta_equal=TRUE,query_axis_reversed=TRUE)
    }
  }
}
write.table(do.call(rbind,audit),args[[4]],sep='\t',row.names=FALSE,quote=FALSE)
cat('INDEPENDENT_SINGLER_FULL_SOURCE_REVERSED_QUERY_REPLAY_PASS\n')
