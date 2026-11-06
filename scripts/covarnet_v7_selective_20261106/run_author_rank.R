# Independent rank jobs: author-normalized within-parent frequencies, 30 nsNMF starts.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L);root<-normalizePath(args[1]);k<-as.integer(args[2])
suppressPackageStartupMessages({library(NMF);library(dplyr);library(reshape2);library(igraph)})
cfg<-jsonlite::fromJSON(file.path(root,'config','paper_config.json'));stopifnot(k%in%cfg$ranks)
out<-file.path(root,'rank_survey',sprintf('K%02d',k));stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
definitions<-parse(file.path(root,'author','CoVarNet.R'));e<-new.env(parent=globalenv())
for(name in c('freq_normalize')){
 i<-which(vapply(definitions,function(x)is.call(x)&&length(x)>=3&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name)),logical(1)));stopifnot(length(i)==1L);eval(definitions[[i]],e)
}
fq<-read.delim(file.path(root,'inputs','frequency_PRIVATE.tsv'),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
X<-e$freq_normalize(fq,normalize='minmax');stopifnot(all(is.finite(X)),min(X)>=0,max(X)<=1,nrow(X)>=max(cfg$ranks))
start<-format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC');seed<-as.integer(cfg$seed_base+k)
fit<-nmf(X,rank=k,method='nsNMF',nrun=as.integer(cfg$nrun),seed=seed,.options='-p',maxIter=as.integer(cfg$maxIter),theta=cfg$theta)
saveRDS(fit,file.path(out,'model_PRIVATE.rds'))
W<-basis(fit);H<-coef(fit);S<-smoothing(NMF::fit(fit));stopifnot(max(abs(W%*%S%*%H-as.matrix(fitted(fit))))<1e-7)
metric<-data.frame(rank=k,cophenetic=as.numeric(cophcor(fit)),RSS=sum((X-as.matrix(fitted(fit)))^2),nrun=as.integer(NMF::nrun(fit)),seed=seed)
stopifnot(is.finite(metric$cophenetic),metric$nrun==cfg$nrun);write.table(metric,file.path(out,'rank_metrics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='PASS_rank_fit_written',UTC_start=start,UTC_end=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),rank=k,nrun=metric$nrun,n_nodes=nrow(X),n_samples=ncol(X),theta=cfg$theta,maxIter=cfg$maxIter,NMF_version=as.character(packageVersion('NMF')),cophenetic=metric$cophenetic,effective_WS_reconstruction_verified=TRUE),file.path(out,'complete.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
