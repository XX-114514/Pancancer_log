# Controlled K14 diagnostics: each resample/seed still uses 30 nsNMF starts.
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);base<-normalizePath(args[2]);jobid<-args[3]
suppressPackageStartupMessages({library(NMF);library(dplyr);library(reshape2);library(igraph)})
jobs<-jsonlite::fromJSON(file.path(root,'jobs.json'));j<-jobs[jobs$job==jobid,,drop=FALSE];stopifnot(nrow(j)==1L);out<-file.path(root,'fits',jobid);dir.create(out);start<-format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC')
f<-read.delim(file.path(base,'inputs/frequency_PRIVATE.tsv'),row.names=1,check.names=FALSE);src<-read.delim(file.path(root,'sample_source_PRIVATE.tsv'),row.names=1,check.names=FALSE);keep<-colnames(f)[-1];if(j$kind=='LOSO')keep<-keep[src[keep,'gse_id']!=j$held_source];f<-f[,c('majorCluster',keep),drop=FALSE];raw<-as.matrix(f[,-1]);span<-apply(raw,1,max)-apply(raw,1,min);removed<-rownames(raw)[span==0];f<-f[span>0,,drop=FALSE];stopifnot(nrow(f)>14)
def<-parse(file.path(base,'author/CoVarNet.R'));e<-new.env(parent=globalenv())
for(name in c('freq_normalize','pair_correlation','cm_network')){
 i<-which(vapply(def,function(x)is.call(x)&&length(x)>=3&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name)),logical(1)));x<-def[[i]]
 if(name=='cm_network'){s<-paste(deparse(x,width.cutoff=500),collapse='\n');s<-sub('rownames(weight) = as.list(1:nrow(weight))','rownames(weight) = as.character(seq_len(nrow(weight)))',s,fixed=TRUE);x<-parse(text=s)[[1]]};eval(x,e)
};e$subgraph<-igraph::induced_subgraph;e$delete.vertices<-igraph::delete_vertices
X<-e$freq_normalize(f,normalize='minmax');fit<-nmf(X,rank=14,method='nsNMF',nrun=30,seed=as.integer(j$seed),.options='-p',maxIter=2000,theta=.5);colnames(basis(fit))<-sprintf('F%02d',1:14);rownames(coef(fit))<-sprintf('F%02d',1:14);saveRDS(fit,file.path(out,'model_PRIVATE.rds'));W<-basis(fit);A<-W%*%smoothing(NMF::fit(fit));H<-coef(fit);stopifnot(max(abs(A%*%H-as.matrix(fitted(fit))))<1e-7)
for(name in c('W','A','H'))write.table(get(name),file.path(out,paste0(name,'_PRIVATE.tsv')),sep='\t',quote=FALSE,col.names=NA)
pairs<-e$pair_correlation(f,method='pearson');net<-e$cm_network(fit,pairs,top_n=15,corr=.2,pval_fdr=.05);write.table(net$each$node,file.path(out,'each_node.tsv'),sep='\t',quote=FALSE,row.names=FALSE);write.table(net$each$edge,file.path(out,'each_edge.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='PASS_controlled_K14_refit',UTC_start=start,UTC_end=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),job=jobid,kind=j$kind,held_source=j$held_source,seed=j$seed,nrun=NMF::nrun(fit),n_samples=ncol(X),n_nodes=nrow(X),removed_constant_nodes=removed,CCC=as.numeric(cophcor(fit)),rank_reoptimized=FALSE),file.path(out,'complete.json'),auto_unbox=TRUE,pretty=TRUE)
