# Conditional-on-frozen-input robustness for completed composition CM fits.
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
suppressPackageStartupMessages(library(NMF));cfg<-jsonlite::fromJSON(args[1],simplifyVector=FALSE)
stopifnot(as.character(packageVersion('NMF'))==cfg$expected_NMF_version)
write_json<-function(x,p)jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,digits=16,na='null')
cosine<-function(a,b){an<-sqrt(colSums(a*a));bn<-sqrt(colSums(b*b));z<-crossprod(a,b)/outer(an,bn);z[!is.finite(z)]<-0;z}
kl_project<-function(A,X,maxiter=2000L){
 H<-matrix(1/ncol(A),ncol(A),ncol(X));den<-pmax(colSums(A),1e-12);old<-Inf
 for(it in seq_len(maxiter)){
  H<-H*(crossprod(A,X/pmax(A%*%H,1e-12))/den)
  if(it%%20L==0L){Y<-pmax(A%*%H,1e-12);pos<-X>0;loss<-sum(Y-X)+sum(X[pos]*log(X[pos]/Y[pos]));if(is.finite(old)&&abs(old-loss)<1e-8*(1+abs(old)))break;old<-loss}
 }
 stopifnot(all(is.finite(H)),all(H>=0));list(H=H,iterations=it)
}
fit_parts<-function(X,weights,k,seed,nrun,maxiter){
 keep<-rowSums(X)>0;stopifnot(sum(keep)>k,sum(weights>0)>k+4L)
 Y<-sweep(X[keep,weights>0,drop=FALSE],2,weights[weights>0], '*');model<-nmf(Y,k,method='nsNMF',seed=seed,nrun=nrun,maxIter=maxiter,.options='-p')
 w<-matrix(0,nrow(X),k);w[keep,]<-as.matrix(basis(model));s<-smoothing(NMF::fit(model));a<-w%*%s
 stopifnot(max(abs(a[keep,,drop=FALSE]%*%as.matrix(coef(model))-as.matrix(fitted(model))))<1e-8)
 list(W=w,A=a,n_unobserved_training_nodes=sum(!keep))
}
all_results<-list();nr<-0L
for(scheme in cfg$schemes){
 root<-scheme$root;rc<-jsonlite::fromJSON(file.path(root,'config/spec_v1.json'));acc<-jsonlite::fromJSON(rc$input_acceptance$path)
 fq<-read.delim(acc$frequency$path,row.names=1,check.names=FALSE,stringsAsFactors=FALSE);meta<-read.delim(acc$metadata$path,check.names=FALSE,stringsAsFactors=FALSE)
 stopifnot(identical(colnames(fq)[-1],meta$sampleID),!anyNA(meta$patient_candidate_block),!anyNA(meta$source))
 X<-as.matrix(fq[,-1,drop=FALSE]);span<-apply(X,1,max)-apply(X,1,min);stopifnot(all(span>0));X<-sweep(sweep(X,1,apply(X,1,min),'-'),1,span,'/')
 k<-as.integer(cfg$rank);D<-file.path(root,'outputs',paste0('seed_',cfg$baseline_seed));base<-list(W=as.matrix(read.delim(file.path(D,'private_W_v1.tsv'),row.names=1,check.names=FALSE)),A=as.matrix(read.delim(file.path(D,'private_A_v1.tsv'),row.names=1,check.names=FALSE)))
 stopifnot(ncol(base$W)==k,identical(rownames(base$W),rownames(fq)))
 baselineH<-kl_project(base$A,X)$H
 blocks<-unique(meta$patient_candidate_block);blockid<-match(meta$patient_candidate_block,blocks);bsize<-as.numeric(table(blockid)[as.character(blockid)]);sources<-sort(unique(meta$source))
 tasks<-list(list(kind='source_equal',identifier='all_sources',weights=1/as.numeric(table(meta$source)[meta$source])))
 set.seed(as.integer(cfg$seed))
 for(b in seq_len(cfg$n_bootstrap)){
  draws<-sample(seq_along(blocks),length(blocks),replace=TRUE);multi<-tabulate(draws,nbins=length(blocks));tasks[[length(tasks)+1L]]<-list(kind='candidate_block_bootstrap',identifier=paste0('repeat_',b),weights=multi[blockid]/bsize)
 }
 for(s in sources)tasks[[length(tasks)+1L]]<-list(kind='leave_one_source',identifier=s,weights=ifelse(meta$source==s,0,1/bsize))
 # Matching must be one-to-one; unmatched factors are retained explicitly.
 for(ti in seq_along(tasks)){
  task<-tasks[[ti]];weights<-task$weights;if(sum(weights>0)<12L){nr<-nr+1L;all_results[[nr]]<-data.frame(scheme=scheme$id,kind=task$kind,identifier=task$identifier,status='insufficient_training_units');next}
  weights<-weights/mean(weights[weights>0]);p<-fit_parts(X,weights,k,as.integer(cfg$seed+ti),as.integer(cfg$nrun),as.integer(cfg$max_iter));projection<-kl_project(p$A,X)
  for(basis in c('A','W')){
   co<-cosine(base[[basis]],p[[basis]]);# padding allows every baseline factor to remain unmatched
   benefit<-matrix(0,2*k,2*k);benefit[seq_len(k),seq_len(k)]<-co-cfg$minimum_cosine;benefit<-benefit-min(benefit)
   match<-as.integer(clue::solve_LSAP(benefit,maximum=TRUE))
   for(i in seq_len(k)){
    j<-match[i];ok<-j<=k&&co[i,j]>=cfg$minimum_cosine
    hc<-if(ok&&sd(baselineH[i,])>0&&sd(projection$H[j,])>0)cor(baselineH[i,],projection$H[j,],method='spearman')else NA_real_
    held<-meta$source==task$identifier;heldrho<-if(ok&&task$kind=='leave_one_source'&&sum(held)>=10&&sd(baselineH[i,held])>0&&sd(projection$H[j,held])>0)cor(baselineH[i,held],projection$H[j,held],method='spearman')else NA_real_
    nr<-nr+1L;all_results[[nr]]<-data.frame(scheme=scheme$id,kind=task$kind,identifier=task$identifier,status=if(ok)'matched'else'unmatched',basis=basis,baseline_module=i,refit_module=if(ok)j else NA_integer_,loading_cosine=if(ok)co[i,j]else NA_real_,projected_activity_Spearman=hc,held_source_activity_Spearman=heldrho,n_training_units=sum(weights>0),n_held_source_units=if(task$kind=='leave_one_source')sum(held)else 0,n_unobserved_training_nodes=p$n_unobserved_training_nodes,projection_iterations=projection$iterations)
   }
  }
  write_json(list(scheme=scheme$id,n_tasks_completed=ti,n_tasks_planned=length(tasks)),file.path(cfg$output,paste0('checkpoint_',scheme$id,'.json')))
 }
}
# Fill schema columns for deliberately nonestimable tasks.
keys<-unique(unlist(lapply(all_results,names)));all_results<-lapply(all_results,function(x){for(n in setdiff(keys,names(x)))x[[n]]<-NA;x[,keys,drop=FALSE]});result<-do.call(rbind,all_results)
write.table(result,file.path(cfg$output,'private_CM_block_source_robustness_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(status='FROZEN_MATRIX_CM_BLOCK_AND_SOURCE_REFITS_COMPLETE',n_schemes=length(cfg$schemes),rank=k,n_bootstrap=cfg$n_bootstrap,minimum_cosine=cfg$minimum_cosine,matching='one_to_one_Hungarian_with_unmatched',NMF_version=as.character(packageVersion('NMF')),normalization_frozen_before_refits=TRUE,candidate_blocks_not_certified_independent_patients=TRUE,heldout_atlas_and_annotation_not_refitted=TRUE,clinical_outcomes_used=FALSE),file.path(cfg$output,'CM_block_source_robustness_report_v1.json'))
