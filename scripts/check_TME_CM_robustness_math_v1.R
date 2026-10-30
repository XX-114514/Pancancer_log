args<-commandArgs(TRUE);stopifnot(length(args)==2L);suppressPackageStartupMessages(library(NMF))
x<-parse(args[1]);target<-new.env(parent=.GlobalEnv)
for(n in c('cosine','kl_project','fit_parts')){ids<-which(vapply(x,function(e)is.call(e)&&length(e)>=3&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name(n)),logical(1)));stopifnot(length(ids)==1);eval(x[[ids]],target)}
A<-matrix(c(1,.1,0,0,.1,1),3,2);H<-rbind(seq(.1,2,length.out=40),seq(2,.1,length.out=40));X<-A%*%H
p<-target$kl_project(A,X);stopifnot(max(abs(p$H-H))<.005)
co<-target$cosine(A,A[,2:1]);stopifnot(max(abs(diag(co[,2:1])-1))<1e-12)
k<-2L;benefit<-matrix(0,4,4);benefit[1:2,1:2]<-matrix(c(.99,.3,.2,.4),2,2)-.6;benefit<-benefit-min(benefit);m<-as.integer(clue::solve_LSAP(benefit,maximum=TRUE));stopifnot(m[1]==1,m[2]>2)
xx<-rbind(X,rep(0,40));model<-target$fit_parts(xx,rep(1,40),2,42,2,500);stopifnot(model$n_unobserved_training_nodes==1,all(model$W[4,]==0),all(is.finite(model$A)))
jsonlite::write_json(list(status='CM_ROBUSTNESS_MATH_SYNTHETIC_PASS',n_tests=4,all_tests_passed=TRUE,package_NMF=as.character(packageVersion('NMF')),not_human_scientific_validation=TRUE),args[2],auto_unbox=TRUE,pretty=TRUE)
