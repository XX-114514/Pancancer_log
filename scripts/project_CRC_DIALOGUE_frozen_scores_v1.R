# Reconstruct all discovery NNLS scores, then project untouched heldout cells.
args<-commandArgs(trailingOnly=TRUE);root<-args[1]
source(args[2],local=new.env(parent=globalenv()))
library(Matrix)
r<-readRDS(file.path(root,'private_discovery_model_v1/result.rds'))
out<-file.path(root,'private_frozen_projection_v1');dir.create(out,showWarnings=FALSE)
checks<-list();gene_names<-list();patient_scores<-list();slot_scope<-list()
for(type in c('Tnode9','Macnode18')){
 read_partition<-function(part){
  p<-file.path(root,'private_inputs',paste0(type,'_',part));md<-read.delim(paste0(p,'_metadata.tsv'),check.names=FALSE);g<-readLines(paste0(p,'_genes.tsv'));m<-readMM(gzfile(paste0(p,'_tpm.mtx.gz')));rownames(m)<-g;colnames(m)<-md$cell;list(md=md,m=m)
 }
 tr<-read_partition('discovery');te<-read_partition('holdout');stopifnot(length(unique(tr$md$sample))==26,length(unique(te$md$sample))==12,length(intersect(tr$md$sample,te$md$sample))==0)
 for(program in c('MCP1','MCP2')){
  gp<-r$gene.pval[[type]];gp<-gp[gp$program==program,];active<-gp[is.finite(gp$coef)&gp$coef>0,];slot_scope[[paste(type,program)]]<-data.frame(cell_type=type,program=program,n_gene_evidence=nrow(gp),n_active_NNLS_genes=nrow(active),gene_score_supported=nrow(active)>0)
  if(nrow(active)==0)next
  stopifnot(!anyDuplicated(active$genes),all(active$genes%in%rownames(tr$m)),all(active$genes%in%rownames(te$m)))
  train<-as.matrix(tr$m[active$genes,,drop=FALSE]);test<-as.matrix(te$m[active$genes,,drop=FALSE]);mu<-rowMeans(train);sigma<-apply(train,1,sd);stopifnot(all(sigma>0))
  weight<-active$coef*ifelse(active$up,1,-1);zt<-sweep(sweep(train,1,mu,'-'),1,sigma,'/');zh<-sweep(sweep(test,1,mu,'-'),1,sigma,'/');raw<-as.numeric(crossprod(weight,zt));raw_h<-as.numeric(crossprod(weight,zh));fit<-lm(raw~cellQ+geneQ,data=tr$md);score<-residuals(fit);held<-raw_h-predict(fit,newdata=te$md)
  saved<-r$scores[[type]];match_rows<-match(tr$md$cell,saved$cells);stopifnot(!anyNA(match_rows));expected<-saved[match_rows,program];err<-max(abs(score-expected));stopifnot(err<1e-6,all(is.finite(score)),all(is.finite(held)))
  checks[[paste(type,program)]]<-data.frame(cell_type=type,program=program,n_discovery_cells=length(score),n_holdout_cells=length(held),max_saved_score_error=err,n_weighted_genes=nrow(active),training_only_normalization_and_quality_fit=TRUE)
  for(part in c('discovery','holdout')){
   md<-if(part=='discovery')tr$md else te$md;s<-if(part=='discovery')score else held;v<-data.frame(sample=md$sample,score=s);v<-aggregate(score~sample,data=v,FUN=mean);v$cell_type<-type;v$program<-program;v$partition<-part;patient_scores[[paste(type,program,part)]]<-v
  }
  saveRDS(list(genes=active$genes,weight=weight,gene_mean=mu,gene_sd=sigma,quality_coefficients=coef(fit)),file.path(out,paste0(type,'_',program,'_fixed_weights.rds')))
  gene_names[[paste(type,program)]]<-data.frame(cell_type=type,program=program,genes=active$genes,direction=ifelse(active$up,'up','down'))
 }
}
write.table(do.call(rbind,checks),file.path(root,'frozen_projection_numeric_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,slot_scope),file.path(root,'gene_supported_program_slots.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,gene_names),file.path(root,'discovery_active_gene_signatures.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,patient_scores),file.path(out,'all_patient_mean_scores_private.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('DISCOVERY_SCORE_PARITY_AND_FROZEN_HOLDOUT_PROJECTION_PASS\n')
