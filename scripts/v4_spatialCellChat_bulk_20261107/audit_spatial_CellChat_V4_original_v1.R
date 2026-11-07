args<-commandArgs(trailingOnly=TRUE);r<-normalizePath(args[1]);unit<-args[2];p<-file.path(r,'imaging',unit,'CellChat_v1');dirs<-list.dirs(p,recursive=FALSE);checked<-0;maxerr<-0
for(d in dirs){
 f<-file.path(d,'probability_arrays_PRIVATE.rds');if(!file.exists(f))next
 z<-readRDS(f);t<-read.delim(file.path(d,'identity_pair_counts.tsv'),check.names=FALSE);edges<-read.delim(file.path(d,'candidate_LR.tsv'),check.names=FALSE);stopifnot(!anyDuplicated(t[,c('source','target')]),all(is.finite(z$prob)),all(z$prob>=0),all(is.finite(z$pval)),all(z$pval>=0 & z$pval<=1),nrow(t)==dim(z$prob)[1]*dim(z$prob)[2])
 q<-array(p.adjust(as.vector(z$pval),'BH'),dim=dim(z$pval),dimnames=dimnames(z$pval));for(i in seq_len(nrow(t))){s<-as.character(t$source[i]);a<-as.character(t$target[i]);ex<-sum(z$prob[s,a,]>0 & z$pval[s,a,]<.05);eq<-sum(z$prob[s,a,]>0 & q[s,a,]<.05);stopifnot(t$n_candidate_CCI[i]==ex,t$n_BH_CCI[i]==eq,sum(edges$source==s & edges$target==a)==ex)}
 if(nrow(edges))for(i in seq_len(nrow(edges))){ix<-cbind(match(edges$source[i],dimnames(z$prob)[[1]]),match(edges$target[i],dimnames(z$prob)[[2]]),match(edges$interaction[i],dimnames(z$prob)[[3]]));stopifnot(abs(z$prob[ix]-edges$prob[i])<1e-10,abs(z$pval[ix]-edges$p[i])<1e-12,abs(q[ix]-edges$BH_q[i])<1e-12)}
 checked<-checked+1
}
jsonlite::write_json(list(status=if(checked)'PASS_all_available_native_spatial_CellChat_saved_arrays_counts_and_full_BH' else 'UNASSESSABLE_NO_FITTED_NATIVE_SPATIAL_CELLCHAT_VARIANT',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),n_variants_checked=checked,patient_contact_or_causality_claim=FALSE),file.path(p,'independent_audit.json'),pretty=TRUE,auto_unbox=TRUE)
