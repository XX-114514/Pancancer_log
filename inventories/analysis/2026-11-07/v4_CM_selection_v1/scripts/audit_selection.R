args<-commandArgs(trailingOnly=TRUE);r<-normalizePath(args[1]);suppressPackageStartupMessages(library(NMF));jobs<-jsonlite::fromJSON(file.path(r,'jobs.json'))
for(j in jobs$job){
 p<-file.path(r,'fits',j);f<-readRDS(file.path(p,'model_PRIVATE.rds'));z<-jsonlite::fromJSON(file.path(p,'complete.json'));C<-consensus(f);D<-as.dist(1-C);CCC<-cor(as.vector(D),as.vector(cophenetic(hclust(D,method='average'))));stopifnot(NMF::nrun(f)==30,abs(CCC-z$CCC)<1e-10,ncol(basis(f))==14,max(abs(basis(f)%*%smoothing(NMF::fit(f))%*%coef(f)-as.matrix(fitted(f))))<1e-7)
 for(name in c('W','A','H')){expected<-switch(name,W=basis(f),A=basis(f)%*%smoothing(NMF::fit(f)),H=coef(f));saved<-as.matrix(read.delim(file.path(p,paste0(name,'_PRIVATE.tsv')),row.names=1,check.names=FALSE));stopifnot(identical(dim(saved),dim(expected)),max(abs(saved-expected))<1e-8)}
}
jsonlite::write_json(list(status='PASS_independent_R_all_18_models_540_starts_CCC_WSH_and_exports',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),n_jobs=nrow(jobs),n_starts_each=30,rank=14),file.path(r,'diagnostics/independent_R_audit.json'),auto_unbox=TRUE,pretty=TRUE)
