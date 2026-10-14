#!/usr/bin/env Rscript
# Original CoVarNet cellular composition/NMF/network method with bounded rank
# survey, source-centered correlations and one independent patient block per
# network observation. All row-level results remain on the private server.
options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly=TRUE)
if(!(length(args)%in%c(3L,4L,5L)))stop('Usage: run_immune_CM.R SCOPE_DIR COVARNET_R_DIR CODE_SHA [ATTEMPT [RESUME_DISCOVERY_DIR]]')
scope<-args[[1]];rdir<-args[[2]];code_sha<-args[[3]]
attempt<-if(length(args)>=4L)args[[4]]else 'discovery_v1'
if(!grepl('^discovery_v[0-9]+$',attempt))stop('Invalid immutable discovery attempt')
out<-file.path(scope,attempt)
if(dir.exists(out))stop('Existing discovery preserved')
dir.create(out,recursive=TRUE)
for(pkg in c('NMF','dplyr','reshape2','igraph','jsonlite'))if(!requireNamespace(pkg,quietly=TRUE))stop('Missing ',pkg)
suppressPackageStartupMessages({library(NMF);library(dplyr);library(reshape2);library(igraph)})
lines<-readLines(file.path(rdir,'R','CoVarNet.R'),warn=FALSE)
lines<-lines[!grepl('^\\s*usethis::use_package\\(',lines)]
eval(parse(text=paste(lines,collapse='\n')),envir=.GlobalEnv)
fq<-read.delim(file.path(scope,'inputs','frequency_private.tsv'),check.names=FALSE,row.names=1)
meta<-read.delim(file.path(scope,'inputs','sample_metadata_private.tsv'),check.names=FALSE)
meta$patient_verified<-tolower(as.character(meta$patient_verified))%in%c('true','1')
meta<-meta[match(colnames(fq)[-1],meta$sampleID),,drop=FALSE]
if(anyNA(meta$sampleID))stop('CM metadata coverage incomplete')
original<-as.matrix(fq[,-1,drop=FALSE]);storage.mode(original)<-'double'
if(any(!is.finite(original))||any(original<0))stop('Invalid complete-denominator frequencies')
matrix_norm<-freq_normalize(fq,normalize='minmax')
if(any(!is.finite(matrix_norm)))stop('Invalid minmax normalization')
# Balance discovery over verified patients; unverified libraries are explicitly
# exploratory units, never falsely described as independent patients.
blocks<-unique(meta$independent_unit)
balanced<-sapply(blocks,function(b)rowMeans(matrix_norm[,meta$independent_unit==b,drop=FALSE]))
rownames(balanced)<-rownames(matrix_norm)
meta_block<-meta[match(blocks,meta$independent_unit),,drop=FALSE]
max_rank<-min(12L,nrow(balanced)-1L,floor(ncol(balanced)/4))
if(max_rank<2L)stop('Insufficient independent composition blocks')
profiles<-list();models<-list()
resume<-if(length(args)>=5L)args[[5]]else ''
if(nzchar(resume)){
 profile<-read.delim(file.path(resume,'rank_survey.tsv'))
 if(!identical(as.integer(profile$rank),as.integer(2:max_rank)))stop('Resume rank survey scope differs')
}else for(k in 2:max_rank){
 message('rank survey K=',k,' blocks=',ncol(balanced))
 fit<-nmf(balanced,k,method='nsNMF',nrun=10,seed=20261013+k,.options='v',maxIter=1000)
 summ<-summary(fit)
 cophenetic<-if('cophenetic'%in%names(summ))as.numeric(summ[['cophenetic']])else as.numeric(cophcor(fit))
 profiles[[as.character(k)]]<-data.frame(rank=k,cophenetic=cophenetic,rss=sum((balanced-fitted(fit))^2))
 models[[as.character(k)]]<-fit
}
if(!nzchar(resume))profile<-do.call(rbind,profiles)
write.table(profile,file.path(out,'rank_survey.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(any(!is.finite(profile$cophenetic)))stop('NMF consensus stability unavailable')
# Pre-specified stability heuristic: largest rank before the first material
# consensus drop (>0.02); if no drop, the survey boundary is reported explicitly.
drops<-which(diff(profile$cophenetic)<(-0.02))
selected<-if(length(drops))profile$rank[min(drops)]else max(profile$rank)
message('selected discovery rank K=',selected)
if(nzchar(resume)){
 fit<-readRDS(file.path(resume,'NMF_discovery.rds'))
 oldW<-as.matrix(read.delim(file.path(resume,'W_private.tsv'),row.names=1,check.names=FALSE))
 oldH<-as.matrix(read.delim(file.path(resume,'H_patient_block_private.tsv'),row.names=1,check.names=FALSE))
 if(!identical(rownames(basis(fit)),rownames(balanced))||!identical(colnames(coef(fit)),colnames(balanced))||ncol(basis(fit))!=selected)stop('Resume model identity or selected rank differs')
 if(!identical(dim(oldW),dim(basis(fit)))||!identical(dim(oldH),dim(coef(fit)))||max(abs(oldW-basis(fit)))>1e-8||max(abs(oldH-coef(fit)))>1e-8)stop('Resume saved model/table parity failed')
 if(any(!is.finite(basis(fit)))||any(!is.finite(coef(fit)))||any(basis(fit)<0)||any(coef(fit)<0))stop('Invalid resume factors')
 message('Reusing completed immutable main fit; rerun diagnostics only')
}else fit<-nmf(balanced,selected,method='nsNMF',nrun=30,seed=20261013,.options='v',maxIter=2000)
colnames(basis(fit))<-sprintf('CM%02d',seq_len(selected));rownames(coef(fit))<-colnames(basis(fit))
saveRDS(fit,file.path(out,'NMF_discovery.rds'))
write.table(basis(fit),file.path(out,'W_private.tsv'),sep='\t',quote=FALSE,col.names=NA)
write.table(coef(fit),file.path(out,'H_patient_block_private.tsv'),sep='\t',quote=FALSE,col.names=NA)
# Average frequencies within verified patient blocks before source-adjusted
# network correlations: repeated sampling does not inflate nominal n.
block_freq<-sapply(blocks,function(b)rowMeans(original[,meta$independent_unit==b,drop=FALSE]))
residual<-t(block_freq)
design<-model.matrix(~0+factor(meta_block$dataset))
residual<-qr.resid(qr(design),residual)
cor_mat<-cor(residual);N<-ncol(cor_mat);df<-nrow(residual)-qr(design)$rank-1L
if(df<3L)stop('Insufficient source-adjusted network degrees of freedom')
pair<-which(upper.tri(cor_mat),arr.ind=TRUE)
r<-cor_mat[pair];r<-pmax(-.999999,pmin(.999999,r))
p<-2*pt(-abs(r*sqrt(df/(1-r^2))),df=df)
corres<-data.frame(subCluster1=rownames(cor_mat)[pair[,1]],subCluster2=colnames(cor_mat)[pair[,2]],
                  correlation=r,pval=p,pval_fdr=p.adjust(p,'BH'))
major<-setNames(as.character(fq$majorCluster),rownames(fq))
corres$majorCluster1<-major[corres$subCluster1];corres$majorCluster2<-major[corres$subCluster2]
corres$spe<-vapply(seq_len(nrow(pair)),function(i){a<-pair[i,1];b<-pair[i,2];bg<-c(cor_mat[a,-a],cor_mat[-c(a,b),b]);mean(bg<=cor_mat[a,b],na.rm=TRUE)},numeric(1))
# Require observable positive support in at least two studies for pooled edges.
# Single-study scopes retain explicitly within-cohort exploratory networks.
source_effects<-list();corres$n_positive_sources<-0L
for(study in unique(meta_block$dataset)){
 keep<-meta_block$dataset==study
 if(sum(keep)<5L)next
 cm<-cor(t(block_freq[,keep,drop=FALSE]));values<-cm[pair]
 corres$n_positive_sources<-corres$n_positive_sources+as.integer(is.finite(values)&values>0)
 source_effects[[study]]<-data.frame(source=study,subCluster1=corres$subCluster1,subCluster2=corres$subCluster2,correlation=values,n_blocks=sum(keep))
}
corres$source_support_eligible<-corres$n_positive_sources>=if(length(unique(meta_block$dataset))>1L)2L else 1L
corres$pval_fdr[!corres$source_support_eligible|!is.finite(corres$pval_fdr)]<-1
network<-cm_network(fit,corres,top_n=min(10L,nrow(balanced)),corr=.2,pval_fdr=.05)
saveRDS(network,file.path(out,'CM_network.rds'))
write.table(corres,file.path(out,'source_adjusted_pair_correlations_private.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(source_effects))write.table(do.call(rbind,source_effects),file.path(out,'per_source_effects_private.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
for(name in c('global','each'))for(type in c('node','edge'))write.table(network[[name]][[type]],file.path(out,paste0(name,'_',type,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
write.table(network$filter,file.path(out,'CM_members.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Source hold-out stability refits learn new W for diagnostics, never count as
# fixed-reference external validation. Compare all components by cosine matrix.
stability<-list();holdout_status<-list();W<-basis(fit)
cosine<-function(A,B){A<-sweep(A,2,sqrt(colSums(A^2)),'/');B<-sweep(B,2,sqrt(colSums(B^2)),'/');t(A)%*%B}
for(study in unique(meta_block$dataset)){
 keep<-meta_block$dataset!=study
 null_rows<-rowSums(balanced[,keep,drop=FALSE])<=0
 reason<-if(sum(keep)<max(8L,selected+1L))'insufficient_training_blocks'else if(any(null_rows))'unobservable_nodes_after_source_holdout'else 'estimable'
 holdout_status[[study]]<-data.frame(heldout_source=study,n_training_blocks=sum(keep),n_unobservable_nodes=sum(null_rows),status=reason)
 if(reason!='estimable'){
   message('LOSO not estimable for ',study,': ',reason,'; missing nodes are not biological zero')
   next
 }
 message('LOSO refit ',study)
 held<-nmf(balanced[,keep,drop=FALSE],selected,method='nsNMF',nrun=10,seed=20261013,.options='v',maxIter=1000)
 sim<-cosine(W,basis(held));write.table(sim,file.path(out,paste0('LOSO_cosine_',make.names(study),'.tsv')),sep='\t',quote=FALSE,col.names=NA)
 stability[[study]]<-data.frame(heldout_source=study,cm=colnames(W),best_cosine=apply(sim,1,max),one_to_one_assignment_not_assumed=TRUE)
}
if(length(stability))write.table(do.call(rbind,stability),file.path(out,'LOSO_stability.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,holdout_status),file.path(out,'LOSO_scope_status.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
png(file.path(out,'rank_survey.png'),width=1400,height=1000,res=150)
plot(profile$rank,profile$cophenetic,type='b',xlab='Rank',ylab='Consensus cophenetic');abline(v=selected,col='red',lty=2);dev.off()
png(file.path(out,'CM_W_heatmap.png'),width=1800,height=1600,res=160)
heatmap(W,scale='none',Colv=NA,cexRow=.35);dev.off()
report<-list(observed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),code_sha=code_sha,
 scope=basename(scope),selected_rank=selected,rank_selection='largest_before_first_0.02_consensus_drop_or_survey_boundary',
 rank_at_survey_boundary=selected==max_rank,n_nodes=nrow(W),n_complete_denominator_units=ncol(original),
 n_composition_blocks=length(blocks),n_verified_patient_blocks=length(unique(meta$patient_block[meta$patient_verified])),
 n_unverified_library_units=sum(!meta$patient_verified),n_sources=length(unique(meta$dataset)),
 n_network_edges=nrow(network$global$edge),n_network_nodes=nrow(network$global$node),
 n_modules_with_network_members=length(unique(network$each$node$cm)),
 new_discovery_not_fixed_reference_validation=TRUE,clinical_outcomes_used=FALSE,
 repeated_patient_network_averaged=TRUE,source_centered_correlations=TRUE,
 LOSO_refit_count=length(stability),LOSO_matching='full_cosine_matrix;row_best_is_diagnostic_not_one_to_one_proof',
 LOSO_unestimable_count=sum(vapply(holdout_status,function(x)x$status!='estimable',logical(1))),
 resumed_completed_main_fit=nzchar(resume),
 authority='candidate_discovery_modules;clinical_validation_and_biological_review_pending')
jsonlite::write_json(report,file.path(out,'CM_public.json'),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(out,'session_info_private.txt'))
message('CM scope complete: ',basename(scope))
