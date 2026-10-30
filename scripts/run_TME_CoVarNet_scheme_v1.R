# Frozen author-CoVarNet baseline; inputs are already complete-denominator views.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
for(pkg in c('NMF','dplyr','reshape2','igraph','psych','jsonlite'))if(!requireNamespace(pkg,quietly=TRUE))stop('Missing package: ',pkg)
suppressPackageStartupMessages({library(NMF);library(dplyr);library(reshape2);library(igraph)})
cfg <- jsonlite::fromJSON(args[[1]],simplifyVector=TRUE)
out <- cfg$output
stopifnot(dir.exists(out))
definitions <- parse(cfg$author_code)
author <- new.env(parent=.GlobalEnv)
for(name in c('freq_normalize','pair_correlation','cm_network')) {
  index <- which(vapply(definitions,function(e)is.call(e)&&length(e)>=3L&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name(name)),logical(1)))
  stopifnot(length(index)==1L)
  expression <- definitions[[index]]
  if(name=='cm_network') {
    text <- paste(deparse(expression,width.cutoff=500L),collapse='\n')
    old <- 'rownames(weight) = as.list(1:nrow(weight))'
    stopifnot(length(gregexpr(old,text,fixed=TRUE)[[1]])==1L,grepl(old,text,fixed=TRUE))
    text <- sub(old,'rownames(weight) = as.character(seq_len(nrow(weight)))',text,fixed=TRUE)
    expression <- parse(text=text)[[1]]
  }
  eval(expression,envir=author)
}
# API-name compatibility only; keep the author's edge selection intact.
author$subgraph <- igraph::induced_subgraph
author$delete.vertices <- igraph::delete_vertices
empty_network_schema <- function(network) {
  if(!is.data.frame(network$global$node)) {
    stopifnot(length(network$global$node)==0L)
    network$global$node<-data.frame(subCluster=character(),majorCluster=character())
  }
  stopifnot(all(c('subCluster','majorCluster')%in%colnames(network$global$node)))
  network
}
fq <- read.delim(cfg$frequency,check.names=FALSE,row.names=1,stringsAsFactors=FALSE)
meta <- read.delim(cfg$metadata,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot('majorCluster'%in%colnames(fq),identical(colnames(fq)[1],'majorCluster'),!anyDuplicated(rownames(fq)),!anyDuplicated(meta$sampleID))
meta <- meta[match(colnames(fq)[-1],meta$sampleID),,drop=FALSE]
stopifnot(!anyNA(meta$sampleID),identical(meta$sampleID,colnames(fq)[-1]))
raw <- as.matrix(fq[,-1,drop=FALSE]);storage.mode(raw)<-'double'
stopifnot(all(is.finite(raw)),all(raw>=0),all(raw<=1+1e-12),nrow(raw)>=4L,ncol(raw)>=12L)
minimum <- apply(raw,1,min);span <- apply(raw,1,max)-minimum
stopifnot(all(span>0))
X <- author$freq_normalize(fq,normalize='minmax')
independent <- sweep(sweep(raw,1,minimum,'-'),1,span,'/')
stopifnot(identical(dim(X),dim(independent)),max(abs(X-independent))<1e-12)
mode <- cfg$weight_mode
weights <- rep(1,ncol(X))
if(mode=='patient_candidate_equal') {
  stopifnot('patient_candidate_block'%in%colnames(meta),!anyNA(meta$patient_candidate_block))
  weights <- 1/as.numeric(table(meta$patient_candidate_block)[meta$patient_candidate_block])
} else if(mode=='source_equal') {
  weights <- 1/as.numeric(table(meta$source)[meta$source])
} else stopifnot(mode=='unit')
weights <- weights/mean(weights)
weighted <- sweep(X,2,weights,'*')
write.table(data.frame(sampleID=meta$sampleID,weight=weights),file.path(out,'private_analysis_weights_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(isTRUE(cfg$synthetic_counterexamples)) {
  stopifnot(cfg$analysis_role=='synthetic_validation',all(startsWith(meta$source,'ARTIFICIAL_')))
  renamed<-fq;rownames(renamed)<-paste0('RENAMED_',seq_len(nrow(fq)))
  stopifnot(max(abs(unname(author$freq_normalize(renamed,normalize='minmax'))-unname(X)))<1e-12)
  KL<-function(a,b){z<-b;positive<-a>0;z[positive]<-a[positive]*log(a[positive]/b[positive])-a[positive]+b[positive];sum(z)}
  target<-X+.1
  stopifnot(abs(KL(weighted,sweep(target,2,weights,'*'))-sum(vapply(seq_len(ncol(X)),function(j)weights[j]*KL(X[,j,drop=FALSE],target[,j,drop=FALSE]),numeric(1))))<1e-8)
}
max_rank <- min(as.integer(cfg$maximum_rank),nrow(X)-1L,floor(ncol(X)/4))
ranks <- seq.int(as.integer(cfg$minimum_rank),max_rank)
stopifnot(length(ranks)>0L,min(ranks)>=2L)
profiles <- list()
for(k in ranks) {
  model <- nmf(weighted,k,method='nsNMF',nrun=as.integer(cfg$survey_nrun),seed=as.integer(cfg$survey_seed+k),.options='-p',maxIter=as.integer(cfg$survey_max_iter))
  co <- as.numeric(cophcor(model));stopifnot(length(co)==1L,is.finite(co))
  profiles[[as.character(k)]] <- data.frame(rank=k,cophenetic=co,RSS=sum((weighted-fitted(model))^2))
}
survey <- do.call(rbind,profiles)
write.table(survey,file.path(out,'rank_survey_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
drops <- which(diff(survey$cophenetic)<(-cfg$cophenetic_drop))
selected <- if(length(drops))survey$rank[min(drops)]else max(survey$rank)
if(!is.null(cfg$rank_override)) {
  stopifnot(cfg$rank_override%in%ranks);selected <- as.integer(cfg$rank_override)
}
corres <- author$pair_correlation(fq,method='pearson')
stopifnot(nrow(corres)==choose(nrow(raw),2))
# The exact author psych defaults may already adjust its upper triangle.
# Preserve that compatibility result and separately export raw Pearson/BH.
corres$author_pval <- corres$pval;corres$author_pval_fdr <- corres$pval_fdr
r <- pmin(1,pmax(-1,corres$correlation));df <- ncol(raw)-2
statistic <- abs(r)*sqrt(df/pmax(0,1-r*r))
corres$raw_Pearson_p <- 2*pt(-statistic,df=df)
corres$raw_Pearson_BH <- p.adjust(corres$raw_Pearson_p,'BH')
write.table(corres,file.path(out,'private_all_pairs_author_and_raw_BH_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
parts <- function(model) {
  W<-as.matrix(basis(model));Hw<-as.matrix(coef(model));S<-smoothing(NMF::fit(model));A<-W%*%S
  stopifnot(max(abs(A%*%Hw-as.matrix(fitted(model))))<1e-8)
  H<-sweep(Hw,2,weights,'/');list(W=W,S=S,A=A,H=H,H_weighted=Hw)
}
for(seed in as.integer(cfg$final_seeds)) {
  dir <- file.path(out,paste0('seed_',seed));stopifnot(!dir.exists(dir));dir.create(dir)
  model <- nmf(weighted,selected,method='nsNMF',nrun=as.integer(cfg$final_nrun),seed=seed,.options='-p',maxIter=as.integer(cfg$final_max_iter))
  labels <- sprintf('CM%02d',seq_len(selected));colnames(basis(model))<-labels;rownames(coef(model))<-labels
  p <- parts(model);saveRDS(model,file.path(dir,'private_NMF_model_v1.rds'))
  if(isTRUE(cfg$synthetic_counterexamples)) {
    empty<-empty_network_schema(author$cm_network(model,corres,top_n=min(as.integer(cfg$top_n),nrow(X)),corr=1,pval_fdr=0))
    stopifnot(nrow(empty$global$edge)==0L,nrow(empty$global$node)==0L,nrow(empty$each$node)==0L)
  }
  for(name in names(p))write.table(p[[name]],file.path(dir,paste0('private_',name,'_v1.tsv')),sep='\t',quote=FALSE,col.names=NA)
  for(flavor in c('author_compatibility','raw_Pearson_BH_sensitivity')) {
    cor <- corres
    if(flavor=='raw_Pearson_BH_sensitivity'){cor$pval<-cor$raw_Pearson_p;cor$pval_fdr<-cor$raw_Pearson_BH}
    network <- empty_network_schema(author$cm_network(model,cor,top_n=min(as.integer(cfg$top_n),nrow(X)),corr=cfg$minimum_positive_correlation,pval_fdr=cfg$network_FDR))
    saveRDS(network,file.path(dir,paste0('private_',flavor,'_network_v1.rds')))
    for(level in c('global','each'))for(kind in c('node','edge'))write.table(network[[level]][[kind]],file.path(dir,paste0(flavor,'_',level,'_',kind,'_v1.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  }
  jsonlite::write_json(list(seed=seed,rank=selected,effective_WSH_reconstruction_checked=TRUE),file.path(dir,'seed_terminal_v1.json'),auto_unbox=TRUE,pretty=TRUE)
}
report <- list(status='COVARNET_SCHEME_NMF_AND_AUTHOR_NETWORK_EXECUTION_COMPLETE',analysis_role=cfg$analysis_role,selected_rank=selected,rank_at_survey_boundary=selected==max_rank,rank_selection='first_0.02_cophenetic_drop_or_boundary_unless_frozen_common_rank',n_nodes=nrow(X),n_analysis_observations=ncol(X),n_sources=length(unique(meta$source)),weight_mode=mode,weight_sum=sum(weights),n_final_seed_models=length(cfg$final_seeds),NMF_version=as.character(packageVersion('NMF')),psych_version=as.character(packageVersion('psych')),psych_default_adjust=as.character(formals(psych::corr.test)$adjust),network_spe_cutoff_is_author_top_n_dependent=TRUE,author_empty_table_serialization_only_fix=TRUE,raw_W_and_effective_WS_separately_saved=TRUE,clinical_outcomes_used=FALSE,rows_are_not_automatically_independent_patients=TRUE,author_pval_fdr_not_renamed_as_unadjusted_raw_BH=TRUE,independent_source_patient_validation_not_claimed=TRUE)
jsonlite::write_json(report,file.path(out,'CoVarNet_execution_report_v1.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('COVARNET_SCHEME_COMPLETE\n');flush.console()
