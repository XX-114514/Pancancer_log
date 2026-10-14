# Fixed two original identities, patient discovery/holdout. No clinical label use.
args<-commandArgs(trailingOnly=TRUE)
root<-args[1]
source(args[2],local=new.env(parent=globalenv()))
library(Matrix)
types<-c('Tnode9','Macnode18')
rA<-lapply(types,function(type){
 prefix<-file.path(root,'private_inputs',paste0(type,'_discovery'))
 md<-read.delim(paste0(prefix,'_metadata.tsv'),check.names=FALSE)
 genes<-readLines(paste0(prefix,'_genes.tsv'))
 X<-as.matrix(read.delim(paste0(prefix,'_PC.tsv'),row.names=1,check.names=FALSE))
 tpm<-as.matrix(readMM(gzfile(paste0(prefix,'_tpm.mtx.gz'))))
 stopifnot(nrow(tpm)==length(genes),ncol(tpm)==nrow(md),identical(rownames(X),md$cell),all(is.finite(tpm)))
 rownames(tpm)<-genes;colnames(tpm)<-md$cell
 make.cell.type(name=type,tpm=tpm,samples=md$sample,X=X,metadata=data.frame(geneQ=md$geneQ,row.names=md$cell),cellQ=md$cellQ)
})
out<-file.path(root,'private_discovery_model_v1')
dir.create(out,showWarnings=FALSE)
param<-DLG.get.param(k=2,results.dir=paste0(out,'/'),seed1=20261014,abn.c=20,n.genes=200,conf=c('cellQ','geneQ'),covar=c('cellQ','geneQ','tme.qc'),pheno=NULL,plot.flag=FALSE)
result<-DIALOGUE.run(rA,main='CRC_FrozenPair015',param=param,plot.flag=FALSE)
saveRDS(result,file.path(out,'result.rds'))
summary<-data.frame(status=if(is.null(result$scores))'NO_PROGRAMS_OR_REQUIRES_REVIEW' else 'DISCOVERY_FIT_DONE_HOLDOUT_PENDING',n_discovery_patients=26,n_holdout_patients=12,nonempty_MCPs=if(is.null(result$MCPs))0 else sum(vapply(result$MCPs,function(x)!is.null(x)&&length(x)>0,logical(1))),clinical_outcomes_used=FALSE)
write.table(summary,file.path(root,'discovery_fit_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('REAL_CRC_DISCOVERY_TERMINAL_HOLDOUT_NOT_YET_VALIDATED\n')
