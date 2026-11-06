args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);cfg<-jsonlite::fromJSON(file.path(root,'config/paper_config.json'));suppressPackageStartupMessages(library(NMF));out<-file.path(root,'CM_results')
fq<-read.delim(file.path(root,'inputs/frequency_PRIVATE.tsv'),row.names=1,check.names=FALSE);x<-as.matrix(fq[,-1]);mins<-apply(x,1,min);span<-apply(x,1,max)-mins;X<-sweep(sweep(x,1,mins,'-'),1,span,'/');actual<-list()
for(k in cfg$ranks){
 folder<-file.path(root,'rank_survey',sprintf('K%02d',k));f<-readRDS(file.path(folder,'model_PRIVATE.rds'));z<-read.delim(file.path(folder,'rank_metrics.tsv'));C<-consensus(f);D<-as.dist(1-C);co<-cor(as.vector(D),as.vector(cophenetic(hclust(D,method='average'))));stopifnot(abs(co-z$cophenetic)<1e-10,NMF::nrun(f)==30,ncol(basis(f))==k,ncol(coef(f))==ncol(X),identical(rownames(basis(f)),rownames(X)),max(abs(basis(f)%*%smoothing(NMF::fit(f))%*%coef(f)-as.matrix(fitted(f))))<1e-7,abs(z$RSS-sum((X-as.matrix(fitted(f)))^2))<1e-7);actual[[as.character(k)]]<-co
}
rho<-unlist(actual);eligible<-cfg$ranks[vapply(cfg$ranks,function(k)all(as.character(c(k-2,k-1,k,k+1))%in%names(rho))&&rho[as.character(k-2)]<rho[as.character(k-1)]&&rho[as.character(k-1)]<rho[as.character(k)]&&rho[as.character(k)]>rho[as.character(k+1)],logical(1))];receipt<-jsonlite::fromJSON(file.path(out,'complete.json'));n<-length(eligible)
if(!n){stopifnot(receipt$status=='NO_ADMISSIBLE_RANK_UNDER_PAPER_CCC_RULE');selected<-NULL}else{
 selected<-eligible[order(-rho[as.character(eligible)],eligible)][1];stopifnot(selected==receipt$selected_rank);f<-readRDS(file.path(out,'selected_model_PRIVATE.rds'));W<-basis(f);H<-coef(f);S<-smoothing(NMF::fit(f));A<-W%*%S
 for(name in c('W','H','S','A')){saved<-as.matrix(read.delim(file.path(out,paste0(name,'_PRIVATE.tsv')),row.names=1,check.names=FALSE));stopifnot(identical(dim(saved),dim(get(name))),max(abs(saved-get(name)))<1e-8)}
 net<-readRDS(file.path(out,'author_network_PRIVATE.rds'));nodes<-rownames(W);adj<-matrix(0,length(nodes),length(nodes),dimnames=list(nodes,nodes));edges<-net$global$edge
 if(nrow(edges)){a<-match(edges$subCluster1,nodes);b<-match(edges$subCluster2,nodes);adj[cbind(a,b)]<-1;adj[cbind(b,a)]<-1}
 topology<-read.delim(file.path(out,'CM_network_topology_and_permutation.tsv'));set.seed(cfg$permutation_seed)
 for(cm in colnames(W)){
  selectednodes<-as.character(net$each$node$subCluster[net$each$node$cm==cm]);size<-length(selectednodes);observed<-if(size>=2)sum(adj[selectednodes,selectednodes,drop=FALSE])/2/choose(size,2)else 0
  null<-numeric(cfg$permutations)
  if(size>=2)for(j in seq_len(cfg$permutations)){draw<-sample.int(length(nodes),size);null[j]<-sum(adj[draw,draw,drop=FALSE])/2/choose(size,2)}
  z<-topology[topology$cm==cm,,drop=FALSE];p<-(1+sum(null>=observed))/(1+cfg$permutations);stopifnot(nrow(z)==1,z$n_network_nodes==size,abs(z$connectivity-observed)<1e-12,abs(z$permutation_P-p)<1e-12,abs(z$null_mean-mean(null))<1e-12)
 }
}
jsonlite::write_json(list(status='PASS_independent_R_saved_models_CCC_WSH_and_permutations',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),n_rank_models=length(cfg$ranks),n_starts_per_rank=30,paper_admissible_ranks=as.integer(eligible),selected_rank=selected,no_rank_forced=TRUE,all_CCC_recomputed_from_saved_consensus=TRUE,all_WS_H_reconstructions_verified=TRUE,all_topology_permutations_recomputed=TRUE,all_atlas_or_biological_validation_claim=FALSE),file.path(out,'independent_R_audit.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
