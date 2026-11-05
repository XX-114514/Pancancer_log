suppressPackageStartupMessages({library(Matrix);library(jsonlite)})
a<-commandArgs(TRUE);out<-a[1];reg<-fromJSON(file.path(out,'registry.json'),simplifyVector=FALSE)
set.seed(191);gs<-unique(unlist(lapply(reg$panels,function(p)c(unlist(p$positive_genes),unlist(p$expected_lower_genes)))))
gs<-c(gs,paste0('SYNTH',seq_len(2100-length(gs))));cc<-paste0('synthetic',seq_len(360));src<-rep(c('sourceA','sourceB'),each=180);kind<-rep(rep(1:3,each=60),2)
x<-matrix(rpois(length(gs)*360,.15),nrow=length(gs),dimnames=list(gs,cc))
for(k in 1:3){sig<-list(c('CD3D','CD3E','TRAC','CD4','FOXP3','IL2RA','IKZF2'),c('CD3D','CD3E','TRAC','CD8A','CD8B','GZMK','NKG7','GNLY'),c('MS4A1','CD79A','CD79B','CD19','CD22','TCL1A'))[[k]];x[sig,kind==k]<-matrix(rpois(length(sig)*sum(kind==k),6),nrow=length(sig))
 block<-paste0('SYNTH',seq((k-1)*60+1,k*60));x[block,kind==k]<-matrix(rpois(60*sum(kind==k),3),nrow=60)
}
av<-matrix(TRUE,nrow=length(gs),ncol=2,dimnames=list(gs,unique(src)));av['IL1RL1','sourceB']<-FALSE;x['IL1RL1',src=='sourceB']<-0
md<-data.frame(source=src,truth_for_synthetic_test_only=kind,row.names=cc)
saveRDS(list(counts=as(x,'dgCMatrix'),metadata=md,availability=av,full_library_counts=TRUE,source_column='source'),file.path(out,'input.rds'))
q<-data.frame(cell_id=cc,keep=!(seq_along(cc)%in%c(179,180,359,360)),reason='synthetic_keep')
q$reason[!q$keep]<-'synthetic_QC_exclusion_for_test'
write.table(q,file.path(out,'synthetic_QC.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write_json(list(registry=file.path(out,'registry.json'),seed=191,n_hvg=500,n_pcs=20,graph_pcs=10,resolutions=c(.2,.4),marker_resolution=.4,min_model_source_coverage=1,run_ucell=FALSE,reviewed_keep_file=file.path(out,'synthetic_QC.tsv')),file.path(out,'config.json'),auto_unbox=TRUE)
