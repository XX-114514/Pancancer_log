args<-commandArgs(TRUE);r<-normalizePath(args[1]);library(jsonlite)
d<-read.csv(file.path(r,'analysis_inputs_PRIVATE.csv'),row.names=1,check.names=FALSE)
stopifnot(!anyDuplicated(rownames(d)));out<-file.path(r,'models_PRIVATE');stopifnot(!dir.exists(out));dir.create(out)
scores<-c('CM03','CM04','CM07','CM09','CM11','Vascular_B_NK_equal_mean','Vascular_B_NK_min_component');records<-list();designs<-list();i<-0
fit_one<-function(score,scope,incremental=FALSE){
 cols<-c(score,'panT','HR','pCR','treatment');if(incremental)cols<-unique(c(cols,'EC','mean_B_NK'))
 z<-d[,cols,drop=FALSE];z<-z[complete.cases(z),,drop=FALSE];if(scope!='interaction')z<-z[z$treatment==if(scope=='pembro')1 else 0,,drop=FALSE]
 x<-cbind(Intercept=rep(1,nrow(z)),score=z[[score]],panT=z$panT,HR=z$HR)
 target<-'score'
 if(incremental)x<-cbind(x,EC=z$EC,mean_B_NK=z$mean_B_NK)
 if(scope=='interaction'){
  if(incremental){tx<-x[,-1,drop=FALSE]*z$treatment;colnames(tx)<-paste0(colnames(tx),'_treatment');x<-cbind(x,treatment=z$treatment,tx);target<-'score_treatment'}
  else{x<-cbind(x,treatment=z$treatment,score_treatment=z[[score]]*z$treatment);target<-'score_treatment'}
 }
 y<-z$pCR;n<-length(y);minority<-if(n)min(sum(y),sum(1-y)) else 0
 rec<-list(family=if(incremental)'secondary_incremental3' else 'primary21',score=score,target=scope,n=n,events=sum(y),minority=minority,n_parameters=ncol(x),status='unavailable',reason='',beta=NA_real_,SE=NA_real_,OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,LRT_P_diagnostic=NA_real_,condition_number=NA_real_)
 if(n<30||minority<10){rec$reason<-'sample_or_minority_gate';return(rec)}
 if(qr(x)$rank<ncol(x)){rec$reason<-'rank_deficient';return(rec)}
 cond<-kappa(x,exact=TRUE);rec$condition_number<-cond
 if(!is.finite(cond)||cond>1e6){rec$reason<-'conditioning_gate';return(rec)}
 warnings<-character();fit<-withCallingHandlers(glm.fit(x=x,y=y,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 info<-crossprod(x,x*as.numeric(fit$fitted.values*(1-fit$fitted.values)));cov<-tryCatch(solve(info),error=function(e)NULL)
 if(!fit$converged||is.null(cov)||any(!is.finite(fit$coefficients))||max(abs(fit$coefficients))>20||length(warnings)){rec$reason<-'nonconvergence_separation_or_warning';return(rec)}
 j<-match(target,colnames(x));b<-fit$coefficients[j];se<-sqrt(cov[j,j]);if(!is.finite(se)||se>10){rec$reason<-'unstable_SE';return(rec)}
 reduced<-glm.fit(x=x[,-j,drop=FALSE],y=y,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100))
 rec$status<-'estimable';rec$beta<-unname(b);rec$SE<-unname(se);rec$OR<-exp(b);rec$CI_low<-exp(b-1.959963984540054*se);rec$CI_high<-exp(b+1.959963984540054*se);rec$P<-2*pnorm(-abs(b/se));rec$LRT_P_diagnostic<-pchisq(max(0,reduced$deviance-fit$deviance),df=1,lower.tail=FALSE)
 tag<-paste(rec$family,score,scope,sep='__');write.csv(data.frame(patient_key=rownames(z),outcome=y,x,check.names=FALSE),file.path(out,paste0(tag,'.csv')),row.names=FALSE)
 list_record<-list(tag=tag,target_index=j,n=n,p=ncol(x),beta=unname(fit$coefficients),covariance=unname(cov),deviance=fit$deviance,max_abs_gradient=max(abs(crossprod(x,y-fit$fitted.values))))
 write_json(list_record,file.path(out,paste0(tag,'.json')),digits=17,pretty=TRUE,auto_unbox=TRUE);rec
}
for(s in scores)for(scope in c('pembro','control','interaction')){i<-i+1;records[[i]]<-fit_one(s,scope)}
for(scope in c('pembro','control','interaction')){i<-i+1;records[[i]]<-fit_one('Vascular_B_NK_min_component',scope,TRUE)}
t<-do.call(rbind,lapply(records,as.data.frame));t$BH<-t$BY<-NA_real_
for(f in unique(t$family)){ix<-which(t$family==f);p<-ifelse(is.na(t$P[ix]),1,t$P[ix]);t$BH[ix]<-p.adjust(p,'BH');t$BY[ix]<-p.adjust(p,'BY');t$BH[ix[is.na(t$P[ix])]]<-t$BY[ix[is.na(t$P[ix])]]<-NA_real_}
write.table(t,file.path(r,'all24_prespecified_target_terms.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write_json(list(UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),status='PRODUCTION_ENDED_INDEPENDENT_AUDIT_PENDING',R=R.version.string,n_slots=nrow(t),n_estimable=sum(t$status=='estimable'),n_primary_BH=sum(t$family=='primary21'&t$BH<.05,na.rm=TRUE),n_secondary_BH=sum(t$family=='secondary_incremental3'&t$BH<.05,na.rm=TRUE)),file.path(r,'production_complete.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Produced ',nrow(t),' target slots; ',sum(t$status=='estimable'),' estimable.\n',sep='')
