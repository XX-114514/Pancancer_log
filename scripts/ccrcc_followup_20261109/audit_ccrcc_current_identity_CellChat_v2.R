args<-commandArgs(trailingOnly=TRUE);r<-normalizePath(args[1]);unit<-args[2];p<-file.path(r,'imaging',unit,'CellChat_existing_identity');dirs<-list.dirs(p,recursive=FALSE);checked<-0;maxerr<-0
for(d in dirs){
 f<-file.path(d,'probability_arrays_PRIVATE.rds');if(!file.exists(f))next
 z<-readRDS(f);t<-read.delim(file.path(d,'identity_pair_counts.tsv'),check.names=FALSE);edges<-read.delim(file.path(d,'candidate_LR.tsv'),check.names=FALSE);stopifnot(!anyDuplicated(t[,c('source','target')]),all(is.finite(z$prob)),all(z$prob>=0),all(is.finite(z$pval)),all(z$pval>=0 & z$pval<=1),nrow(t)==dim(z$prob)[1]*dim(z$prob)[2])
 q<-array(p.adjust(as.vector(z$pval),'BH'),dim=dim(z$pval),dimnames=dimnames(z$pval));for(i in seq_len(nrow(t))){s<-as.character(t$source[i]);a<-as.character(t$target[i]);ex<-sum(z$prob[s,a,]>0 & z$pval[s,a,]<.05);eq<-sum(z$prob[s,a,]>0 & q[s,a,]<.05);stopifnot(t$n_candidate_CCI[i]==ex,t$n_BH_CCI[i]==eq,sum(edges$source==s & edges$target==a)==ex)}
 if(nrow(edges))for(i in seq_len(nrow(edges))){ix<-cbind(match(edges$source[i],dimnames(z$prob)[[1]]),match(edges$target[i],dimnames(z$prob)[[2]]),match(edges$interaction[i],dimnames(z$prob)[[3]]));stopifnot(abs(z$prob[ix]-edges$prob[i])<1e-10,abs(z$pval[ix]-edges$p[i])<1e-12,abs(q[ix]-edges$BH_q[i])<1e-12)}
 checked<-checked+1
}
jsonlite::write_json(list(status=if(checked)'PASS_all_available_native_spatial_CellChat_saved_arrays_counts_and_full_BH' else 'UNASSESSABLE_NO_FITTED_NATIVE_SPATIAL_CELLCHAT_VARIANT',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),n_variants_checked=checked,patient_contact_or_causality_claim=FALSE),file.path(p,'independent_audit_v2.json'),pretty=TRUE,auto_unbox=TRUE)

f<-file.path(p,'all_variants_finite_MC_sensitivity.tsv')
stopifnot(file.exists(f),checked==4)
a<-read.delim(f);tol<-1e-12
p_error<-max(abs(a$plus_one_MC_sensitivity_P-(round(a$raw_native_P*100)+1)/101))
bh_error<-max(abs(a$BH_joint_all_variants-p.adjust(a$plus_one_MC_sensitivity_P,'BH')))
by_error<-max(abs(a$BY_joint_all_variants-p.adjust(a$plus_one_MC_sensitivity_P,'BY')))
lower_bound_shortfall<-max(0,1/101-min(a$plus_one_MC_sensitivity_P))
stopifnot(all(is.finite(a$plus_one_MC_sensitivity_P)),all(a$plus_one_MC_sensitivity_P>=1/101-tol),all(a$plus_one_MC_sensitivity_P<=1+tol),p_error<tol,bh_error<tol,by_error<tol)
jsonlite::write_json(list(status='PASS_finite_MC_saved_count_sensitivity_and_joint_adjustments',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),rows=nrow(a),tolerance=tol,minimum_P=min(a$plus_one_MC_sensitivity_P),lower_bound_roundtrip_shortfall=lower_bound_shortfall,max_P_error=p_error,max_BH_error=bh_error,max_BY_error=by_error,model_refitted=FALSE,identity_changed=FALSE),file.path(p,'finite_MC_independent_audit_v2.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
