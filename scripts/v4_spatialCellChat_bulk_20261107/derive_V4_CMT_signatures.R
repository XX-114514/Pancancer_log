args<-commandArgs(trailingOnly=TRUE);r<-normalizePath(args[1]);inp<-file.path(r,'bulk/signature_inputs');out<-file.path(r,'bulk/frozen_signatures');stopifnot(!dir.exists(out));dir.create(out)
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(Matrix));set.seed(2026110703)
x<-readMM(file.path(inp,'mean_CP10K_PRIVATE.mtx'));genes<-readLines(file.path(inp,'genes.tsv'));meta<-read.delim(file.path(inp,'samples_PRIVATE.tsv'),row.names=1,check.names=FALSE);rownames(x)<-genes;colnames(x)<-rownames(meta)
a<-CreateSeuratObject(counts=x,meta.data=meta,min.cells=0,min.features=0);a<-SetAssayData(a,slot='data',new.data=log1p(x));Idents(a)<-'CMT';program<-list();scope<-list()
for(cm in sprintf('CM%02d',1:14)){
 n<-sum(meta$CMT==cm);gg<-data.frame();status<-'TOO_FEW_CMT_SAMPLES'
 if(n>=3 && n<=nrow(meta)-3){
  de<-FindMarkers(a,ident.1=cm,only.pos=TRUE,logfc.threshold=0.25,min.pct=0.1,test.use='wilcox',verbose=FALSE);de$gene<-rownames(de);write.table(de,file.path(out,paste0(cm,'_CMT_DE.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  fc<-if('avg_log2FC'%in%names(de))'avg_log2FC' else 'avg_logFC';de<-de[is.finite(de$p_val_adj)&de$p_val_adj<0.05&de[[fc]]>0.25,,drop=FALSE];de<-de[order(-de[[fc]],de$gene),,drop=FALSE];gg<-head(de,10)
  status<-if(nrow(gg)==10)'FROZEN10GENE_PROGRAM' else 'FEWER_THAN10_QUALIFIED_GENES'
 }
 if(nrow(gg)>0){gg$cm<-cm;program[[cm]]<-gg}
 scope[[cm]]<-data.frame(cm=cm,n_high_samples=n,n_low_samples=nrow(meta)-n,n_signature_genes=nrow(gg),status=status,original_empty_network=cm%in%c('CM01','CM02','CM05'))
}
sig<-if(length(program))do.call(rbind,program) else data.frame(cm=character(),gene=character());write.table(sig,file.path(out,'all14_frozen_CMT_signature_genes.tsv'),sep='\t',quote=FALSE,row.names=FALSE);write.table(do.call(rbind,scope),file.path(out,'all14_signature_eligibility.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='ALL14_CMT_GENE_PROGRAMS_FROZEN_BEFORE_BULK_OUTCOMES',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),Seurat=as.character(packageVersion('Seurat')),n_factors=14,program_not_exact_CM_NMF_activity=TRUE,top10_positive_after_adjP05=TRUE),file.path(out,'complete.json'),pretty=TRUE,auto_unbox=TRUE)
