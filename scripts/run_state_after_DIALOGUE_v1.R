# Official frozen implementation; all14 fixed coverage contexts, no outcome selection.
args<-commandArgs(trailingOnly=TRUE);input<-args[1];out<-args[2];runtime<-args[3]
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE,mode='0700')
.libPaths(c(file.path(runtime,c('library_v7','library_v5','library_v2','library')),.libPaths()))
stopifnot(as.character(packageVersion('lme4'))=='1.1.35.5')
suppressPackageStartupMessages(library(DIALOGUE));suppressPackageStartupMessages(library(Matrix))
scope<-read.delim(file.path(input,'public_input_scope_summary.tsv'),check.names=FALSE);contexts<-sort(unique(scope$context));stopifnot(length(contexts)==14)
units<-read.delim(gzfile(args[4]),check.names=FALSE)
major_totals<-aggregate(n_cells~sample_uid+major,data=units,FUN=sum)
major_keys<-paste(major_totals$sample_uid,major_totals$major,sep='||')
status<-list();slots<-list()
for(context in contexts){
 name<-sprintf('context_%04d',context);root<-file.path(input,name);od<-file.path(out,name);dir.create(od,mode='0700')
 run_one<-function(){
  rA<-lapply(c('A','B'),function(side){
   p<-file.path(root,paste0('type',side,'_discovery'));md<-read.delim(paste0(p,'_metadata.tsv'),check.names=FALSE);genes<-readLines(paste0(p,'_genes.tsv'));pc<-as.matrix(read.delim(paste0(p,'_PC.tsv'),row.names=1,check.names=FALSE));mat<-as.matrix(readMM(gzfile(paste0(p,'_tpm.mtx.gz'))));stopifnot(ncol(mat)==nrow(md),nrow(mat)==length(genes),identical(rownames(pc),md$cell),all(is.finite(mat)),length(unique(md$sample))>=20,all(table(md$sample)>=20));rownames(mat)<-genes;colnames(mat)<-md$cell
   um<-units[match(md$unit_index,units$unit_index),];stopifnot(!anyNA(um$unit_index),identical(as.character(um$sample_uid),as.character(md$sample_uid)))
   frac<-um$n_cells/major_totals$n_cells[match(paste(um$sample_uid,um$major,sep='||'),major_keys)];logN<-log1p(um$n_cells);stopifnot(all(is.finite(frac)),all(frac>0))
   make.cell.type(name=paste0('type',side),tpm=mat,samples=md$sample,X=pc,metadata=data.frame(geneQ=md$geneQ,fraction=frac,logN=logN,row.names=md$cell),cellQ=md$cellQ)
  })
  param<-DLG.get.param(k=2,results.dir=paste0(od,'/'),seed1=20261016,abn.c=20,n.genes=200,conf=c('cellQ','geneQ','fraction','logN'),covar=c('cellQ','geneQ','fraction','logN'),pheno=NULL,plot.flag=FALSE)
  result<-DIALOGUE.run(rA,main=name,param=param,plot.flag=FALSE);saveRDS(result,file.path(od,'discovery_result.rds'))
  nscore<-0
  for(type in c('typeA','typeB')){
   gp<-result$gene.pval[[type]]
   for(program in c('MCP1','MCP2')){
    n<-if(is.null(gp))0 else sum(gp$program==program&is.finite(gp$coef)&gp$coef>0,na.rm=TRUE)
    slots[[paste(context,type,program)]]<<-data.frame(context=context,type=type,program=program,positive_NNLS_genes=n)
    if(n>0)nscore<-nscore+1
   }
  }
  data.frame(context=context,status=if(is.null(result$scores))'NO_PROGRAMS_OR_UNSUPPORTED' else 'DISCOVERY_FIT_DONE',gene_supported_slots=nscore,message=if(is.null(result$message))'' else paste(result$message,collapse=';'),heldout_projection_pending=any(scope$context==context&scope$partition=='holdout'))
 }
 logcon<-file(file.path(od,'fit.log'),open='wt');sink(logcon);sink(logcon,type='message')
 value<-tryCatch(run_one(),error=function(e)data.frame(context=context,status='FAILED_REVIEW_REQUIRED',gene_supported_slots=0,message=conditionMessage(e),heldout_projection_pending=any(scope$context==context&scope$partition=='holdout')))
 sink(type='message');sink();close(logcon);status[[name]]<-value
 write.table(do.call(rbind,status),file.path(out,'public_fit_status.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
 if(length(slots)>0)write.table(do.call(rbind,slots),file.path(out,'public_gene_slot_scope.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 rm(value);gc()
}
writeLines('ALL14_DISCOVERY_CONTEXTS_ATTEMPTED_HOLDOUT_AND_BIOLOGICAL_REVIEW_PENDING',file.path(out,'execution_complete.txt'))
