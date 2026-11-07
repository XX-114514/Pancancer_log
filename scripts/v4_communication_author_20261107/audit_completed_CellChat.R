# Recompute stored global CellChat counts from model probabilities and permutation p.
args<-commandArgs(trailingOnly=TRUE);r<-normalizePath(args[1]);p<-file.path(r,'results/CellChat_global');out<-file.path(r,'audit/CellChat_global_v1');stopifnot(!dir.exists(out));dir.create(out)
suppressPackageStartupMessages(library(CellChat));cc<-readRDS(file.path(p,'model_PRIVATE.rds'));tab<-read.delim(file.path(p,'candidate_LR_by_identity.tsv'),check.names=FALSE);pairs<-read.delim(file.path(p,'identity_pair_counts.tsv'),check.names=FALSE)
stopifnot(all(is.finite(cc@net$prob)),all(cc@net$prob>=0),all(is.finite(cc@net$pval)),all(cc@net$pval>=0 & cc@net$pval<=1))
expected<-apply(cc@net$prob>0 & cc@net$pval<0.05,c(1,2),sum)
stopifnot(identical(dim(expected),dim(cc@net$count)),max(abs(expected-cc@net$count))==0)
stopifnot(!anyDuplicated(pairs[,c('source','target')]),nrow(pairs)==45*45)
for(i in seq_len(nrow(pairs))){s<-pairs$source[i];t<-pairs$target[i];stopifnot(pairs$n_candidate_CCI[i]==expected[s,t]);stopifnot(sum(tab$source==s & tab$target==t)==expected[s,t])}
write.table(pairs,file.path(out,'audited_identity_pair_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
category<-as.data.frame(table(tab$annotation));colnames(category)<-c('category','n_identity_pair_LR_events');write.table(category,file.path(out,'signaling_category_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='PASS_independent_CellChat_saved_model_prob_p_and_complete45_by45_candidate_counts',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),CellChat=as.character(packageVersion('CellChat')),n_identity_pairs=nrow(pairs),n_identity_pair_LR_events=nrow(tab),n_LR_types_in_fitted_model=dim(cc@net$prob)[3],n_eligible_database_interactions=nrow(cc@DB$interaction),p_threshold=0.05,nboot=100,nonnegative_finite_model_and_p_range=TRUE,no_independent_patient_or_contact_or_causal_claim=TRUE),file.path(out,'complete.json'),pretty=TRUE,auto_unbox=TRUE)
