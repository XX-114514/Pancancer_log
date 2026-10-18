# Official fixed DIALOGUE implementation; no clinical outcomes or heldout fitting.
args<-commandArgs(trailingOnly=TRUE)
input<-args[1];out<-args[2];runtime<-args[3];k<-as.integer(args[4]);seed<-as.integer(args[5]);adjust.source<-args[6]=='source_adjusted'
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE,mode='0700')
.libPaths(c(file.path(runtime,c('library_v7','library_v5','library_v2','library')),.libPaths()))
stopifnot(as.character(packageVersion('DIALOGUE'))=='1.0',as.character(packageVersion('lme4'))=='1.1.35.5')
stopifnot(startsWith(normalizePath(find.package('lme4')),normalizePath(file.path(runtime,'library_v5'))))
suppressPackageStartupMessages(library(DIALOGUE));suppressPackageStartupMessages(library(Matrix))
views<-c('R_FOXP3_capture','E_effector_capture','Myeloid')
rA<-lapply(views,function(view){
 p<-file.path(input,'private_inputs',paste0(view,'_discovery'))
 md<-read.delim(paste0(p,'_metadata.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
 genes<-readLines(paste0(p,'_genes.tsv'));pc<-as.matrix(read.delim(paste0(p,'_PC.tsv'),row.names=1,check.names=FALSE))
 mat<-as.matrix(readMM(gzfile(paste0(p,'_tpm.mtx.gz'))))
 stopifnot(nrow(mat)==length(genes),ncol(mat)==nrow(md),identical(rownames(pc),md$cell),all(is.finite(mat)),all(is.finite(pc)))
 stopifnot(length(unique(md$sample))==38,all(table(md$sample)>=20),all(md$partition=='reviewed_case_discovery'))
 stopifnot(setequal(md$dataset,c('GSE286827','GSE294109')))
 rownames(mat)<-genes;colnames(mat)<-md$cell
 make.cell.type(name=view,tpm=mat,samples=md$sample,X=pc,
  metadata=data.frame(geneQ=md$geneQ,sourceCC=md$sourceCC,fraction=md$fraction,logN=md$logN,row.names=md$cell),cellQ=md$cellQ)
})
controls<-c('cellQ','geneQ',if(adjust.source)'sourceCC','fraction','logN')
param<-DLG.get.param(k=k,results.dir=paste0(out,'/'),seed1=seed,abn.c=20,n.genes=200,
 conf=controls,covar=controls,pheno=NULL,plot.flag=FALSE)
# Preserve official intermediate/full/reduced outputs; no cleanup branch.
param$full.version<-TRUE
name<-paste0('RNA_MCP_20261018_K',k,'_S',seed,'_',if(adjust.source)'sourceAdjusted' else 'sourcePermissive')
start<-Sys.time();result<-DIALOGUE.run(rA,main=name,param=param,plot.flag=FALSE)
slots<-list()
for(view in views){
 gp<-result$gene.pval[[view]]
 if(!is.null(gp)&&nrow(gp)>0)stopifnot(all(c('program','coef')%in%colnames(gp)))
 for(program in paste0('MCP',seq_len(k))){
  count<-if(is.null(gp)||nrow(gp)==0)0 else sum(gp$program==program&is.finite(gp$coef)&gp$coef>0,na.rm=TRUE)
  slots[[paste(view,program)]]<-data.frame(view=view,program=program,positive_NNLS_genes=count)
 }
}
write.table(do.call(rbind,slots),file.path(out,'private_all_planned_gene_slots.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
status<-data.frame(status=if(is.null(result$scores))'OFFICIAL_RUN_COMPLETE_NO_SCORES' else 'OFFICIAL_DISCOVERY_FIT_COMPLETE_REVIEW_PENDING',
 k=k,seed=seed,source_mode=if(adjust.source)'source_adjusted' else 'source_permissive_sensitivity',
 discovery_cases=38,TNBC_cells_used_in_fitting=0,views=length(views),
 nonempty_gene_slots=sum(vapply(slots,function(x)x$positive_NNLS_genes>0,logical(1))),
 elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')),
 heldout_or_clinical_validation=FALSE)
write.table(status,file.path(out,'private_fit_execution_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c(paste('R',getRversion()),paste('DIALOGUE',packageVersion('DIALOGUE')),paste('lme4',packageVersion('lme4'))),file.path(out,'private_runtime_versions.txt'))
cat('OFFICIAL_SINGLE_FIT_FINISHED_STABILITY_AND_SOURCE_HOLDOUT_PENDING\n')
