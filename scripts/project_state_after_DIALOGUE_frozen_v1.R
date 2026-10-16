# Frozen reconstruction and internal heldout projection; no clinical/LR selection.
args<-commandArgs(trailingOnly=TRUE);input<-args[1];models<-args[2];out<-args[3];runtime<-args[4];units_path<-args[5]
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE,mode='0700')
.libPaths(c(file.path(runtime,c('library_v7','library_v5','library_v2','library')),.libPaths()))
suppressPackageStartupMessages(library(Matrix))
units<-read.delim(gzfile(units_path),check.names=FALSE)
mt<-aggregate(n_cells~sample_uid+major,data=units,FUN=sum);mk<-paste(mt$sample_uid,mt$major,sep='||')
status<-read.delim(file.path(models,'public_fit_status.tsv'),check.names=FALSE)
checks<-list();scope<-list();scores<-list();genes_out<-list();errors<-list()
read_part<-function(root,side,part){
 p<-file.path(root,paste0('type',side,'_',part));if(!file.exists(paste0(p,'_metadata.tsv')))return(NULL)
 md<-read.delim(paste0(p,'_metadata.tsv'),check.names=FALSE);genes<-readLines(paste0(p,'_genes.tsv'));m<-readMM(gzfile(paste0(p,'_tpm.mtx.gz')));stopifnot(nrow(m)==length(genes),ncol(m)==nrow(md),!anyDuplicated(md$cell));rownames(m)<-genes;colnames(m)<-md$cell
 u<-units[match(md$unit_index,units$unit_index),];stopifnot(!anyNA(u$unit_index),identical(as.character(u$sample_uid),as.character(md$sample_uid)))
 md$fraction<-u$n_cells/mt$n_cells[match(paste(u$sample_uid,u$major,sep='||'),mk)];md$logN<-log1p(u$n_cells)
 stopifnot(all(is.finite(as.matrix(md[,c('cellQ','geneQ','fraction','logN')]))));list(md=md,m=m)
}
for(context in sort(status$context)){
 name<-sprintf('context_%04d',context);root<-file.path(input,name);path<-file.path(models,name,'discovery_result.rds');if(!file.exists(path)){errors[[name]]<-data.frame(context=context,status='NO_DISCOVERY_RESULT',message='See original fit status; no reconstruction possible');next};r<-readRDS(path)
 for(side in c('A','B')){
  type<-paste0('type',side);tr<-read_part(root,side,'discovery');te<-read_part(root,side,'holdout');stopifnot(length(unique(tr$md$sample))>=20)
  if(!is.null(te))stopifnot(length(intersect(tr$md$sample,te$md$sample))==0,identical(rownames(tr$m),rownames(te$m)))
  for(program in c('MCP1','MCP2')){
   key<-paste(context,type,program,sep='_');gp<-r$gene.pval[[type]];active<-if(is.null(gp))NULL else gp[gp$program==program&is.finite(gp$coef)&gp$coef>0,,drop=FALSE]
   n<-nrow(active);if(is.null(n))n<-0
   scope[[key]]<-data.frame(context=context,type=type,program=program,active_genes=n,has_internal_holdout=!is.null(te),heldout_blocks=if(is.null(te))0 else length(unique(te$md$sample)),status=if(n==0)'NO_GENE_SUPPORTED_SCORE' else 'PENDING_RECONSTRUCTION')
   if(n==0)next
   one<-tryCatch({
    stopifnot(!anyDuplicated(active$genes),all(active$genes%in%rownames(tr$m)))
    train<-as.matrix(tr$m[active$genes,,drop=FALSE]);mu<-rowMeans(train);sigma<-apply(train,1,sd);stopifnot(all(is.finite(sigma)),all(sigma>0));weight<-active$coef*ifelse(active$up,1,-1)
    raw<-as.numeric(crossprod(weight,sweep(sweep(train,1,mu,'-'),1,sigma,'/')))
    design<-model.matrix(~cellQ+geneQ+fraction+logN,data=tr$md);stopifnot(qr(design)$rank==ncol(design))
    fit<-lm(raw~cellQ+geneQ+fraction+logN,data=tr$md);sc<-residuals(fit);saved<-r$scores[[type]];ix<-match(tr$md$cell,saved$cells);stopifnot(!anyNA(ix),program%in%names(saved));err<-max(abs(sc-saved[ix,program]));stopifnot(is.finite(err),err<1e-6,all(is.finite(sc)))
    local_scores<-list();for(part in c('discovery','holdout')){
     z<-if(part=='discovery')tr else te;if(is.null(z))next
     if(part=='discovery')v<-sc else {test<-as.matrix(z$m[active$genes,,drop=FALSE]);rh<-as.numeric(crossprod(weight,sweep(sweep(test,1,mu,'-'),1,sigma,'/')));v<-rh-predict(fit,newdata=z$md)}
     stopifnot(all(is.finite(v)));agg<-aggregate(score~sample,data=data.frame(sample=z$md$sample,score=v),FUN=mean);agg$context<-context;agg$type<-type;agg$program<-program;agg$partition<-part;local_scores[[part]]<-agg
    }
    saveRDS(list(genes=active$genes,weight=weight,gene_mean=mu,gene_sd=sigma,quality_coefficients=coef(fit),confounders=c('cellQ','geneQ','fraction','logN')),file.path(out,paste0(key,'_weights.rds')))
    list(check=data.frame(context=context,type=type,program=program,discovery_cells=length(sc),holdout_cells=if(is.null(te))0 else nrow(te$md),max_saved_score_error=err,active_genes=n,training_only_transform=TRUE),scores=local_scores,genes=data.frame(context=context,type=type,program=program,gene=active$genes,weight=weight))
   },error=function(e)list(error=conditionMessage(e)))
   if(!is.null(one$error)){scope[[key]]$status<-'FAILED_PARITY_OR_DESIGN_REVIEW';errors[[key]]<-data.frame(context=context,status=scope[[key]]$status,message=one$error);next}
   checks[[key]]<-one$check;scores[[key]]<-do.call(rbind,one$scores);genes_out[[key]]<-one$genes;scope[[key]]$status<-'FROZEN_RECONSTRUCTION_PASS'
  }
 }
}
write_rows<-function(x,file){if(length(x)>0)write.table(do.call(rbind,x),file.path(out,file),sep='\t',row.names=FALSE,quote=TRUE)}
if(length(scope)==0)write.table(data.frame(context=integer(),type=character(),program=character(),active_genes=integer(),has_internal_holdout=logical(),heldout_blocks=integer(),status=character()),file.path(out,'public_slot_scope.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
write_rows(scope,'public_slot_scope.tsv');write_rows(checks,'public_reconstruction_checks.tsv');write_rows(errors,'public_projection_errors.tsv');write_rows(genes_out,'public_fixed_gene_weights.tsv');write_rows(scores,'private_block_mean_scores.tsv')
writeLines('ALL_AVAILABLE_SLOTS_ATTEMPTED_REVIEW_SCOPE_AND_ERRORS_NOT_AUTOMATIC_HELDOUT_VALIDATION',file.path(out,'execution_complete.txt'))
