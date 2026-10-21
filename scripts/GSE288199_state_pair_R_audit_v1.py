"""Independent QR/rank/HC3 replay through RAM pipes; never writes input arrays."""
import json
import os
import subprocess

R_CODE = r'''
suppressPackageStartupMessages(library(jsonlite))
input<-fromJSON(file('stdin'),simplifyVector=FALSE)
v<-function(x)vapply(x,function(z)if(is.null(z))NA_real_ else as.numeric(z),numeric(1))
mat<-function(x)do.call(rbind,lapply(x,v))
std<-function(x){s<-sqrt(mean((x-mean(x))^2));if(s==0)s<-1;(x-mean(x))/s}
prepare<-function(x,C){
 n<-length(x);if(n<20L)return(NULL)
 z<-std(rank(x,ties.method='average'))
 if(ncol(C)){rc<-sapply(seq_len(ncol(C)),function(j)std(rank(C[,j],ties.method='average')));raw<-cbind(1,rc)}else raw<-matrix(1,n,1)
 s<-svd(raw,nu=0,nv=0)$d;r<-sum(s>s[1]*1e-12)
 if(s[r]/s[1]<1e-8)return(NULL)
 qr0<-qr(raw,tol=1e-12,LAPACK=FALSE);stopifnot(qr0$rank==r)
 Q<-qr.Q(qr0)[,seq_len(r),drop=FALSE]
 stopifnot(max(abs(raw-Q%*%crossprod(Q,raw)))<1e-9)
 rx<-as.numeric(z-Q%*%crossprod(Q,z))
 if(sqrt(sum(rx^2))<=1e-8*max(1,sqrt(sum(z^2)))||n-r-1L<10L)return(NULL)
 D<-cbind(Q,z);ds<-svd(D,nu=0,nv=0)$d;if(ds[1]/tail(ds,1)>=1e8)return(NULL)
 inverse<-qr.solve(D,diag(n),tol=1e-12);h<-rowSums(D*t(inverse))
 if(max(h)>=.999999)return(NULL)
 list(Q=Q,rx=rx,D=D,inverse=inverse,h=h,h0=rowSums(Q^2))
}
fit<-function(m,y){
 if(is.null(m))return(NULL)
 y<-std(rank(y,ties.method='average'));ry<-as.numeric(y-m$Q%*%crossprod(m$Q,y))
 if(sqrt(sum(ry^2))<=1e-8*max(1,sqrt(sum(y^2))))return(NULL)
 op<-m$inverse[nrow(m$inverse),];res<-as.numeric(y-m$D%*%(m$inverse%*%y))
 se<-sqrt(sum((op*res/(1-m$h))^2));if(!is.finite(se)||se^2<=1e-20)return(NULL)
 beta<-sum(op*y);rho<-sum(m$rx*ry)/sqrt(sum(m$rx^2)*sum(ry^2))
 list(beta=beta,se=se,t=beta/se,rho=max(-1,min(1,rho)),y=y)
}
answers<-list()
for(b in input$blocks){
 x<-v(b$x);y<-v(b$y);C<-mat(b$C);m<-prepare(x,C);f<-fit(m,y)
 expected<-b$expected;stopifnot(is.null(f)==is.null(expected))
 ans<-list(label=b$label,estimable=!is.null(f))
 if(!is.null(f)){
  err<-max(abs(c(f$beta,f$se,f$t,f$rho)-v(expected)));stopifnot(err<1e-8)
  ans$coefficient_SE_t_rho_max_error<-err
  if(!is.null(b$signs)){
   S<-t(mat(b$signs));op<-m$inverse[nrow(m$inverse),]
   e0<-as.numeric(f$y-m$Q%*%crossprod(m$Q,f$y))/(1-m$h0)
   noise<-e0*S;res<-noise-m$D%*%(m$inverse%*%noise)
   se<-sqrt(colSums((op/(1-m$h)*res)^2));wt<-as.numeric(op%*%noise)/se
   delta<-max(abs(wt-v(b$wild)));stopifnot(delta<1e-8)
   ans$all_wild_statistics_max_error<-delta
  }
 }
 if(!is.null(b$indices)){
  values<-vapply(b$indices,function(a){
   ix<-as.integer(v(a))+1L;mm<-prepare(x[ix],C[ix,,drop=FALSE]);ff<-fit(mm,y[ix]);if(is.null(ff))NA_real_ else ff$rho
  },numeric(1))
  exp<-v(b$bootstrap);stopifnot(identical(is.na(values),is.na(exp)))
  delta<-if(any(is.finite(values)))max(abs(values-exp),na.rm=TRUE)else 0
  stopifnot(delta<1e-8);ans$all_bootstrap_rho_max_error<-delta
 }
 answers[[length(answers)+1L]]<-ans
}
cat(toJSON(list(status='PASS_ALL_SUBMITTED_SOURCE_DESIGNS_AND_FULL_RANDOM_STREAM_QR_REPLAY',blocks=answers),auto_unbox=TRUE,digits=16,na='null'))
'''


def audit(blocks):
    payload = json.dumps(dict(blocks=blocks), allow_nan=False)
    result = subprocess.run([os.environ['PANCANCER_R_BLOCK_REFERENCE'], '--vanilla', '-e', R_CODE],
                            input=payload, capture_output=True, text=True, timeout=600)
    if result.returncode:
        raise RuntimeError('Independent_R_QR_replay_failed')
    answer = json.loads(result.stdout)
    assert answer['status'].startswith('PASS_ALL_')
    return answer
