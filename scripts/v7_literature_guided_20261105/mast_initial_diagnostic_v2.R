# Lodi/Qian-guided initial lineage inspection, NOT a final annotation assignment.
# Full source-library raw counts; no cancer, IFN, stress, cell-cycle regression.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(jsonlite);library(ggplot2)})
root <- Sys.getenv('PROJECT_ROOT');stopifnot(nchar(root)>0)
run <- file.path(root,'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1')
input <- file.path(run,'stage_b_mast_fullgene')
out <- file.path(run,'stage_c2_mast_initial_diagnostic')
stopifnot(!dir.exists(out));dir.create(out,recursive=FALSE)
started <- format(Sys.time(),tz='UTC',usetz=TRUE)
set.seed(20261105)
write_json(list(status='running',stage='initial_diagnostic',UTC_start=started,annotation_changes=0),file.path(out,'started.json'),pretty=TRUE)
options(future.globals.maxSize=8*1024^3)
counts <- Read10X_h5(file.path(input,'mast_v7_fullgene_raw_10x.h5'),use.names=TRUE,unique.features=TRUE)
stopifnot(ncol(counts)==45036,nrow(counts)==95038)
meta <- read.delim(gzfile(file.path(input,'mast_v7_metadata_PRIVATE.tsv.gz')),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(identical(colnames(counts),rownames(meta)))
# Keep the entire feature universe through normalization so the library-size
# denominator remains the complete source raw library. VST filters constant genes.
n_genes_detected_ge3 <- sum(Matrix::rowSums(counts>0)>=3)
original_gene <- rownames(counts)
# Seurat converts underscores/pipes and can create duplicate symbols. Preserve
# every raw feature and use explicit unique IDs only for incompatible names.
feature_id <- original_gene
unsafe <- grepl("_|\\|",original_gene)
feature_id[unsafe] <- paste0('V7FEATURE',which(unsafe))
stopifnot(!anyDuplicated(feature_id),!any(grepl("_|\\|",feature_id)))
rownames(counts) <- feature_id
gene_map <- data.frame(feature_id=feature_id,gene_symbol=original_gene)
write.table(gene_map,file.path(out,'feature_id_gene_symbol_map.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
obj <- CreateSeuratObject(counts=counts,meta.data=meta,min.cells=0,min.features=0,project='V7_Mast_initial_review')
obj <- NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
obj <- FindVariableFeatures(obj,selection.method='vst',nfeatures=3000,verbose=FALSE)
hvg <- VariableFeatures(obj)
hvg <- hvg[!grepl('^HB[AB]',hvg)]
hvg <- head(hvg,2000);stopifnot(length(hvg)==2000)
VariableFeatures(obj) <- hvg
write.table(data.frame(feature_id=hvg,gene_symbol=original_gene[match(hvg,feature_id)]),file.path(out,'mast_initial_HVG.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
obj <- ScaleData(obj,features=hvg,verbose=FALSE)
obj <- RunPCA(obj,features=hvg,npcs=50,verbose=FALSE,seed.use=20261105)
write.table(data.frame(PC=seq_along(Stdev(obj)),stdev=Stdev(obj)),file.path(out,'PCA_stdev.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# 30 PCs is a declared exploratory setting, not the paper's dataset-specific elbow choice.
# No technology covariate has been verified for this V7 subset: start unintegrated
# as in the lineage review branch of Qian, and inspect source segregation explicitly.
obj <- FindNeighbors(obj,reduction='pca',dims=1:30,verbose=FALSE)
res <- c(.2,.4,.6,.8,1,1.2,1.6,2)
obj <- FindClusters(obj,resolution=res,algorithm=1,random.seed=20261105,verbose=FALSE)
obj <- RunUMAP(obj,reduction='pca',dims=1:30,n.neighbors=30,n.components=2,seed.use=20261105,verbose=FALSE)
md <- obj[[]]
cols <- grep('^RNA_snn_res',colnames(md),value=TRUE)
summary <- do.call(rbind,lapply(cols,function(cn){z=table(md[[cn]]);data.frame(resolution=cn,cluster=names(z),n_cells=as.integer(z))}))
write.table(summary,file.path(out,'cluster_resolution_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
oldcross <- as.data.frame(table(cluster=md$RNA_snn_res.0.4,old_label=md$final_annotation_v7))
write.table(oldcross[oldcross$Freq>0,],file.path(out,'cluster_oldlabel_crosswalk.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Source tables remain private; public summaries below contain only source counts/fractions.
source_tab <- table(md$RNA_snn_res.0.4,md$gse_id)
sources <- data.frame(cluster=rownames(source_tab),n_gse=rowSums(source_tab>0),max_gse_fraction=apply(source_tab,1,max)/rowSums(source_tab))
write.table(sources,file.path(out,'cluster_source_concentration.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
Idents(obj)<-'RNA_snn_res.0.4';DefaultAssay(obj)<-'RNA'
markers <- FindAllMarkers(obj,assay='RNA',slot='data',min.pct=.10,min.diff.pct=.25,logfc.threshold=.25,return.thresh=.01,only.pos=FALSE,test.use='wilcox',verbose=FALSE)
markers$feature_id <- markers$gene
markers$gene <- original_gene[match(markers$feature_id,feature_id)]
write.table(markers,file.path(out,'mast_initial_RNA_FindAllMarkers_res0.4.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
panel <- c('TPSAB1','TPSB2','CPA3','KIT','MS4A2','HDC','LTC4S','IL1RL1','HPGDS','CMA1','CTSG','GATA2','GATA3','PTGDR2','IL7R','RORC','KLRB1','CD3D','CD3E','TRAC','NKG7','GNLY','LYZ','LST1','TYROBP','CD79A','MS4A1','LILRA4','CLEC4C','FCER1A','CD1C','CLC','SIGLEC8','EPX','PF4','PPBP','HBA1','HBB','GYPA','EPCAM','KRT18','MKI67','TOP2A','ISG15','IFIT1','FOS','JUN','DNAJB1')
panel <- intersect(panel,rownames(obj))
# Explicit expression detection/means for interpretation; availability is separately audited.
mat <- GetAssayData(obj,assay='RNA',layer='data')[panel,,drop=FALSE]
levels <- levels(Idents(obj))
expr <- do.call(rbind,lapply(levels,function(cl){sel=which(Idents(obj)==cl);data.frame(cluster=cl,gene=panel,n_cells=length(sel),pct_expressing=Matrix::rowMeans(mat[,sel,drop=FALSE]>0),mean_log1p_cp10k=Matrix::rowMeans(mat[,sel,drop=FALSE]))}))
write.table(expr,file.path(out,'mast_initial_marker_expression_res0.4.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
qc <- do.call(rbind,lapply(levels,function(cl){ii=which(Idents(obj)==cl);row=data.frame(cluster=cl,n_cells=length(ii));for(key in intersect(c('nCount_RNA','nFeature_RNA','pct_counts_mt','pct_counts_hb','doublet_score','decontx_contamination'),colnames(md))){v=suppressWarnings(as.numeric(md[ii,key]));row[[paste0(key,'_n_available')]]=sum(!is.na(v));row[[paste0(key,'_median')]]=if(all(is.na(v))) NA_real_ else median(v,na.rm=TRUE)};row}))
write.table(qc,file.path(out,'mast_initial_cluster_QC_res0.4.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
pdf(file.path(out,'mast_initial_diagnostic.pdf'),width=13,height=8)
print(ElbowPlot(obj,ndims=50))
print(DimPlot(obj,group.by='RNA_snn_res.0.4',label=TRUE,raster=TRUE)+ggtitle('V7 Mast candidates: initial full-gene lineage model'))
print(DimPlot(obj,group.by='final_annotation_v7',raster=TRUE))
print(DimPlot(obj,group.by='gse_id',raster=TRUE,cols=NULL)+NoLegend())
print(DotPlot(obj,features=panel,group.by='RNA_snn_res.0.4')+RotatedAxis()+ggtitle('Parent, competing lineage and state markers; exploratory'))
print(FeaturePlot(obj,features=intersect(c('TPSAB1','CPA3','KIT','MS4A2','CD3D','LYZ','NKG7','MKI67'),rownames(obj)),ncol=4,raster=TRUE))
dev.off()
saveRDS(obj,file.path(out,'mast_initial_seurat_PRIVATE.rds'),compress=FALSE)
write.table(cbind(cell_id=rownames(md),md),gzfile(file.path(out,'mast_initial_metadata_PRIVATE.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(status='completed_initial_diagnostic_only',UTC_start=started,UTC_end=format(Sys.time(),tz='UTC',usetz=TRUE),n_cells=ncol(obj),n_genes=nrow(obj),HVG=2000,PCA=50,graph_PCs=30,resolutions=res,marker_resolution=.4,normalization='full-library LogNormalize 10000',batch_correction='none; technology not yet verified',regressed_variables=character(),annotation_changes=0,Seurat=as.character(packageVersion('Seurat')),R=R.version.string,limits=c('Only V7 Mast parent included; no outside recruitment yet','Uncleaned initial QC/identity diagnostic, not accepted subtypes','No final subtype assignment from cluster IDs','Source-specific measurement mask required for negative-marker decisions')),file.path(out,'complete.json'),auto_unbox=TRUE,pretty=TRUE)
