# Alternative algorithm on the same balanced global cells; no invented paper parameters.
args<-commandArgs(trailingOnly=TRUE);r<-normalizePath(args[1]);out<-file.path(r,'results','CellChat_global');stopifnot(!dir.exists(out));dir.create(out)
suppressPackageStartupMessages({library(CellChat);library(Matrix)})
if(requireNamespace('future',quietly=TRUE))future::plan('sequential')
options(mc.cores=4);set.seed(2026110702);start<-Sys.time()
x<-Matrix::readMM(file.path(r,'inputs/CellChat_logexpression_PRIVATE.mtx'));genes<-readLines(file.path(r,'inputs/CellChat_genes_PRIVATE.txt'));meta<-read.delim(file.path(r,'inputs/CellChat_meta_PRIVATE.tsv'),row.names=1,check.names=FALSE)
rownames(x)<-genes;colnames(x)<-rownames(meta);stopifnot(nrow(x)==length(genes),ncol(x)==nrow(meta));meta$subCluster<-factor(meta$subCluster)
data('CellChatDB.human',package='CellChat');db<-CellChatDB.human
saveRDS(db,file.path(out,'database_snapshot_PRIVATE.rds'))
cc<-createCellChat(object=x,meta=meta,group.by='subCluster');cc@DB<-db
cc<-subsetData(cc);cc<-identifyOverExpressedGenes(cc);cc<-identifyOverExpressedInteractions(cc)
cc<-computeCommunProb(cc,type='triMean',raw.use=TRUE,population.size=FALSE,nboot=100,seed.use=2026110702)
cc<-filterCommunication(cc,min.cells=10);cc<-computeCommunProbPathway(cc);cc<-aggregateNet(cc)
saveRDS(cc,file.path(out,'model_PRIVATE.rds'))
edges<-subsetCommunication(cc);write.table(edges,file.path(out,'candidate_LR_by_identity.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
count<-cc@net$count;pair<-expand.grid(source=rownames(count),target=colnames(count),stringsAsFactors=FALSE);pair$n_candidate_CCI<-as.vector(count);write.table(pair,file.path(out,'identity_pair_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
annotation<-table(edges$annotation);write.table(data.frame(annotation=names(annotation),n_identity_LR_events=as.integer(annotation)),file.path(out,'signaling_category_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='CELLCHAT_ALTERNATIVE_ALGORITHM_FINISHED_PENDING_REVIEW',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),seconds=as.numeric(difftime(Sys.time(),start,units='secs')),package=as.character(packageVersion('CellChat')),n_cells=ncol(x),n_labels=nlevels(meta$subCluster),nboot=100,min_cells=10,raw_use=TRUE,population_size=FALSE,no_spatial_contact_or_patient_inference=TRUE,paper_exact_version_or_parameters_claim=FALSE),file.path(out,'complete.json'),pretty=TRUE,auto_unbox=TRUE)
