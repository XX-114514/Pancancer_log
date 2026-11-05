# Independent arithmetic checks on the synthetic fixture only.
suppressPackageStartupMessages({library(Matrix);library(Seurat);library(jsonlite)})
a<-commandArgs(TRUE);root<-a[1];out<-file.path(root,'result')
input<-readRDS(file.path(root,'input.rds'));obj<-readRDS(file.path(out,'PRIVATE_lineage_review_object.rds'))
raw<-input$counts[,colnames(obj),drop=FALSE];src<-as.character(input$metadata[colnames(obj),input$source_column]);av<-input$availability
cl<-as.character(Idents(obj));lib<-Matrix::colSums(raw)
e<-read.delim(file.path(out,'measured_marker_expression.tsv'),colClasses=c(cluster='character'),check.names=FALSE)
max_err<-0
for(i in seq_len(nrow(e))){
 g<-e$gene[i];idx<-which(cl==e$cluster[i]);measured<-idx[av[g,src[idx]]]
 stopifnot(e$n_cells[i]==length(idx),e$n_measured[i]==length(measured))
 if(length(measured)){
  pct<-sum(raw[g,measured]>0)/length(measured)
  avg<-mean(as.numeric(raw[g,measured])/lib[measured]*10000)
  stopifnot(abs(pct-e$pct_expressing_measured[i])<1e-10)
  err<-abs(avg-e$mean_cp10k_measured[i]);max_err<-max(max_err,err);stopifnot(err<1e-8)
 }else stopifnot(is.na(e$pct_expressing_measured[i]),is.na(e$mean_cp10k_measured[i]))
}
joint<-read.delim(file.path(out,'same_cell_core_positive_evidence.tsv'),colClasses=c(cluster='character'),check.names=FALSE)
joint<-joint[joint$panel=='T_Treg',]
for(i in seq_len(nrow(joint))){
 idx<-which(cl==joint$cluster[i]&src==joint$source[i])
 for(th in c(1,2)){
  truth<-sum((raw['FOXP3',idx]>=th)&(Matrix::colSums(raw[c('IL2RA','IKZF2','IKZF4'),idx,drop=FALSE]>=th)>0))
  stopifnot(truth==joint[i,paste0('n_positive_raw',th)])
 }
}
cts<-read.delim(file.path(out,'cluster_counts.tsv'));stopifnot(sum(cts$n_cells)==ncol(raw),ncol(raw)==356)
coverage<-read.delim(file.path(out,'cluster_source_counts.tsv'),colClasses=c(cluster='character'))
ill<-e[e$gene=='IL1RL1',]
for(i in seq_len(nrow(ill)))stopifnot(ill$n_measured[i]==sum(coverage$Freq[coverage$cluster==ill$cluster[i]&coverage$source=='sourceA']))
write_json(list(status='PASS',scope='synthetic only',retained_cells=ncol(raw),expected_QC_excluded=4,n_clusters=nrow(cts),expression_rows_checked=nrow(e),max_abs_cp10k_error=max_err,same_cell_Treg_rows_checked=nrow(joint),source_missingness_checked=TRUE,final_cell_labels_assigned=0),file.path(root,'independent_arithmetic_audit.json'),auto_unbox=TRUE,pretty=TRUE)
