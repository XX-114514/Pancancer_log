# Literature-guided evidence producer. Input: full raw-count RDS contract.
# This rebuilds a reviewed lineage subset and outputs evidence, never final labels.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3)stop('Usage: lineage_workflow.R input_contract.rds config.json NEW_OUTPUT_DIR')
input<-readRDS(args[1]);cfg<-fromJSON(args[2],simplifyVector=FALSE);out<-args[3]
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE,mode='0700')
need<-c('counts','metadata','availability','full_library_counts','source_column')
stopifnot(all(need%in%names(input)),isTRUE(input$full_library_counts))
counts<-as(input$counts,'dgCMatrix');md<-input$metadata;avail<-input$availability
stopifnot(identical(colnames(counts),rownames(md)),!anyDuplicated(rownames(counts)),!anyDuplicated(colnames(counts)))
stopifnot(all(is.finite(counts@x)),all(counts@x>=0),all(counts@x==floor(counts@x)))
stopifnot(identical(rownames(avail),rownames(counts)),all(!is.na(avail)),all(avail%in%c(TRUE,FALSE,0,1)))
src<-as.character(md[[input$source_column]]);stopifnot(!anyNA(src),all(src%in%colnames(avail)))
stopifnot(!any(grepl('_|\\|',rownames(counts)))) # require explicit safe feature map upstream
required<-c('registry','seed','n_hvg','n_pcs','graph_pcs','resolutions','marker_resolution','min_model_source_coverage','run_ucell')
stopifnot(all(required%in%names(cfg)))
registry<-fromJSON(cfg$registry,simplifyVector=FALSE)
# Exact versioned aliases only, with collision refusal. No suffix stripping.
old_names<-rownames(counts);new_names<-vapply(old_names,function(g)if(!is.null(registry$aliases[[g]]))registry$aliases[[g]] else g,character(1))
if(anyDuplicated(new_names))stop('Alias collision requires stable-ID adjudication')
rownames(counts)<-new_names;rownames(avail)<-new_names
set.seed(cfg$seed);start<-format(Sys.time(),tz='UTC',usetz=TRUE)
qc_status<-'initial_review_all_cells_retained'
if(!is.null(cfg$reviewed_keep_file)){
 q<-read.delim(cfg$reviewed_keep_file,check.names=FALSE,stringsAsFactors=FALSE)
 stopifnot(all(c('cell_id','keep','reason')%in%names(q)),!anyDuplicated(q$cell_id),setequal(q$cell_id,colnames(counts)))
 q<-q[match(colnames(counts),q$cell_id),];stopifnot(is.logical(q$keep),!anyNA(q$keep),all(nchar(q$reason[!q$keep])>0))
 write.table(as.data.frame(table(keep=q$keep,reason=q$reason)),file.path(out,'reviewed_QC_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 saveRDS(q,file.path(out,'PRIVATE_reviewed_QC_manifest.rds'))
 counts<-counts[,q$keep,drop=FALSE];md<-md[q$keep,,drop=FALSE];src<-src[q$keep]
 qc_status<-'reviewed_keep_manifest_applied_rebuild_from_counts'
}
stopifnot(ncol(counts)>cfg$n_pcs+2)
sources<-sort(unique(src));avail<-avail[,sources,drop=FALSE]
# Validate measurement mask against nonzero counts, no structural-zero inference.
for(s in sources)if(any(Matrix::rowSums(counts[,src==s,drop=FALSE])[!avail[,s]]>0))stop('Nonzero counts contradict source mask')
obj<-CreateSeuratObject(counts=counts,meta.data=md,min.cells=0,min.features=0)
obj<-NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
eligible<-rownames(counts)[rowMeans(avail)>=cfg$min_model_source_coverage]
model<-obj[eligible,]
model<-FindVariableFeatures(model,selection.method='vst',nfeatures=min(cfg$n_hvg,length(eligible)),verbose=FALSE)
hvg<-VariableFeatures(model);stopifnot(length(hvg)>cfg$n_pcs)
obj<-ScaleData(obj,features=hvg,verbose=FALSE) # no automatic cancer/IFN/stress regression
obj<-RunPCA(obj,features=hvg,npcs=cfg$n_pcs,seed.use=cfg$seed,verbose=FALSE)
write.table(data.frame(gene=hvg),file.path(out,'model_HVG.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(PC=seq_along(Stdev(obj)),sd=Stdev(obj)),file.path(out,'PCA_stdev.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
reduction<-'pca'
if(!is.null(cfg$verified_batch_column)){
 if(!requireNamespace('harmony',quietly=TRUE))stop('harmony unavailable')
 stopifnot(cfg$verified_batch_column%in%names(md),!anyNA(md[[cfg$verified_batch_column]]))
 obj<-harmony::RunHarmony(obj,group.by.vars=cfg$verified_batch_column,verbose=FALSE)
 reduction<-'harmony'
}
pcs<-seq_len(cfg$graph_pcs);stopifnot(max(pcs)<=cfg$n_pcs)
obj<-FindNeighbors(obj,reduction=reduction,dims=pcs,verbose=FALSE)
obj<-FindClusters(obj,resolution=unlist(cfg$resolutions),algorithm=1,random.seed=cfg$seed,verbose=FALSE)
obj<-RunUMAP(obj,reduction=reduction,dims=pcs,seed.use=cfg$seed,verbose=FALSE)
key<-paste0('RNA_snn_res.',cfg$marker_resolution);stopifnot(key%in%names(obj[[]]));Idents(obj)<-obj[[key,drop=TRUE]]
cluster<-as.character(Idents(obj));tab<-table(cluster)
write.table(data.frame(cluster=names(tab),n_cells=as.integer(tab)),file.path(out,'cluster_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(as.data.frame(table(cluster=cluster,source=src)),file.path(out,'cluster_source_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Preserve all resolutions as diagnostics; cluster number is never targetted.
reskeys<-grep('^RNA_snn_res',names(obj[[]]),value=TRUE)
rs<-do.call(rbind,lapply(reskeys,function(k){z<-table(obj[[k,drop=TRUE]]);data.frame(resolution=k,cluster=names(z),n_cells=as.integer(z))}))
write.table(rs,file.path(out,'resolution_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
de_one<-function(o,features){
 if(length(unique(Idents(o)))<2)return(data.frame())
 FindAllMarkers(o,assay='RNA',slot='data',features=features,test.use='wilcox',
    min.pct=.10,min.diff.pct=.25,logfc.threshold=.25,return.thresh=.01,
    only.pos=FALSE,min.cells.group=3,verbose=FALSE)
}
# Pooled DE only on genes measured in every represented source. Other genes are
# tested within source. Counts/normalisation retain full-library denominator.
common<-rownames(avail)[rowSums(avail)==ncol(avail)]
pooled<-de_one(obj,common)
write.table(pooled,file.path(out,'RNA_DEG_common_measured.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
source_de<-list()
for(s in sources){
 so<-subset(obj,cells=colnames(obj)[src==s]);g<-rownames(avail)[avail[,s]]
 d<-de_one(so,g);if(nrow(d)){d$source<-s;source_de[[s]]<-d}
}
if(length(source_de))write.table(do.call(rbind,source_de),file.path(out,'RNA_DEG_by_source.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Coverage-aware positive-expression summaries for every registry gene.
gene_union<-unique(unlist(lapply(registry$panels,function(p)c(unlist(p$positive_genes),unlist(p$expected_lower_genes)))))
normal<-GetAssayData(obj,assay='RNA',layer='data');expr<-list();i<-0L
for(cl in names(tab))for(g in gene_union){
 ix<-which(cluster==cl);measured<-if(g%in%rownames(counts))ix[avail[g,src[ix]]] else integer()
 i<-i+1L
 expr[[i]]<-data.frame(cluster=cl,gene=g,n_cells=length(ix),n_measured=length(measured),
  pct_expressing_measured=if(length(measured))sum(counts[g,measured,drop=FALSE]>0)/length(measured) else NA_real_,
  mean_cp10k_measured=if(length(measured))mean(expm1(as.numeric(normal[g,measured]))) else NA_real_)
}
write.table(do.call(rbind,expr),file.path(out,'measured_marker_expression.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
# Same-cell positive evidence, not phenotype totals. Complete-core measurement
# is required for the denominator; failure to express never means identity-negative.
joint<-list()
for(p in registry$panels){
 core<-lapply(p$core_groups,unlist);all_core<-unique(unlist(core))
 for(s in sources){
  cells<-which(src==s);evaluable<-all(all_core%in%rownames(counts))
  if(evaluable)evaluable<-all(avail[all_core,s])
  for(cl in unique(cluster[cells])){
   ix<-cells[cluster[cells]==cl];n<-length(ix)
   vals<-c(NA_integer_,NA_integer_)
   if(evaluable)for(th in c(1L,2L)){
    component_hits<-lapply(core,function(gs)Matrix::colSums(counts[gs,ix,drop=FALSE]>=th)>0)
    vals[th]<-sum(Reduce('&',component_hits))
   }
   joint[[length(joint)+1L]]<-data.frame(panel=p$id,source=s,cluster=cl,n_cells=n,
     n_evaluable_complete_core=if(evaluable)n else 0L,n_positive_raw1=vals[1],n_positive_raw2=vals[2],identity_negative_inferred=FALSE)
  }
 }
}
write.table(do.call(rbind,joint),file.path(out,'same_cell_core_positive_evidence.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
ucell_state<-'not_requested';uagg<-list();coverage<-list()
if(isTRUE(cfg$run_ucell)){
 if(!is.null(cfg$r_library_paths)).libPaths(c(unlist(cfg$r_library_paths),.libPaths()))
 if(!requireNamespace('UCell',quietly=TRUE))stop('UCell requested but official package unavailable; no homemade replacement')
 rankmax<-if(is.null(cfg$ucell_maxRank))1500L else as.integer(cfg$ucell_maxRank)
 for(s in sources){
  genes<-rownames(avail)[avail[,s]];cells<-which(src==s)
  stopifnot(length(genes)>=rankmax)
  signatures<-list()
  for(p in registry$panels){
   signature<-unlist(p$positive_genes);n<-sum(signature%in%genes)
   coverage[[length(coverage)+1L]]<-data.frame(source=s,panel=p$id,n_genes=length(signature),n_measured=n,score_computed=n==length(signature))
   if(n==length(signature))signatures[[p$id]]<-signature
  }
  if(!length(signatures))next
  u<-UCell::ScoreSignatures_UCell(counts[genes,cells,drop=FALSE],features=signatures,maxRank=rankmax,ncores=1)
  saveRDS(u,file.path(out,paste0('PRIVATE_UCell_source_',match(s,sources),'.rds')))
  stopifnot(identical(rownames(u),colnames(counts)[cells]))
  for(cl in unique(cluster[cells]))for(n in colnames(u)){
   vals<-u[cluster[cells]==cl,n]
   uagg[[length(uagg)+1L]]<-data.frame(source=s,cluster=cl,signature=n,n_cells=length(vals),mean=mean(vals),median=median(vals))
  }
 }
 if(length(uagg))write.table(do.call(rbind,uagg),file.path(out,'UCell_within_source_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(do.call(rbind,coverage),file.path(out,'UCell_source_signature_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 ucell_state<-paste0('official_',packageVersion('UCell'),'_complete_signatures_only_no_labels')
}
saveRDS(obj,file.path(out,'PRIVATE_lineage_review_object.rds'))
write_json(list(status='evidence_complete_no_final_labels',start_UTC=start,end_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),
 n_cells=ncol(obj),n_features=nrow(obj),qc=qc_status,reduction=reduction,Seurat=as.character(packageVersion('Seurat')),
 UCell=ucell_state,config=cfg,final_label_changes=0,
 limitations=c('Within-source DE is not independent patient replication','Human-reviewed subtype decision remains necessary','Pooled DEG excludes structurally unmeasured genes','No automatic source/IFN/stress regression','Source masks are not individual-cell dropout masks')),
 file.path(out,'complete.json'),auto_unbox=TRUE,pretty=TRUE)
