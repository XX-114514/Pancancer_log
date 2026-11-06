# Author network functions and paper CCC rule; no result-dependent rank or edge relaxation.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L);root<-normalizePath(args[1]);cfg<-jsonlite::fromJSON(file.path(root,'config','paper_config.json'))
suppressPackageStartupMessages({library(NMF);library(dplyr);library(reshape2);library(igraph)})
out<-file.path(root,'CM_results');stopifnot(!dir.exists(out));dir.create(out)
definition<-parse(file.path(root,'author','CoVarNet.R'));author<-new.env(parent=globalenv())
for(name in c('pair_correlation','cm_network')){
 i<-which(vapply(definition,function(x)is.call(x)&&length(x)>=3L&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name)),logical(1)));stopifnot(length(i)==1L);x<-definition[[i]]
 if(name=='cm_network'){
  s<-paste(deparse(x,width.cutoff=500L),collapse='\n');old<-'rownames(weight) = as.list(1:nrow(weight))';stopifnot(grepl(old,s,fixed=TRUE));s<-sub(old,'rownames(weight) = as.character(seq_len(nrow(weight)))',s,fixed=TRUE);x<-parse(text=s)[[1]]
 }
 eval(x,author)
}
author$subgraph<-igraph::induced_subgraph;author$delete.vertices<-igraph::delete_vertices
survey<-do.call(rbind,lapply(cfg$ranks,function(k)read.delim(file.path(root,'rank_survey',sprintf('K%02d',k),'rank_metrics.tsv'))));stopifnot(identical(as.integer(survey$rank),as.integer(cfg$ranks)));rho<-setNames(survey$cophenetic,survey$rank)
survey$paper_admissible<-vapply(survey$rank,function(k)all(as.character(c(k-2,k-1,k,k+1))%in%names(rho))&&rho[as.character(k-2)]<rho[as.character(k-1)]&&rho[as.character(k-1)]<rho[as.character(k)]&&rho[as.character(k)]>rho[as.character(k+1)],logical(1));write.table(survey,file.path(out,'rank_survey_and_paper_selection.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
eligible<-survey[survey$paper_admissible,,drop=FALSE]
if(!nrow(eligible)){
 jsonlite::write_json(list(status='NO_ADMISSIBLE_RANK_UNDER_PAPER_CCC_RULE',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),no_rank_forced=TRUE,no_primary_CM_claim=TRUE),file.path(out,'complete.json'),auto_unbox=TRUE,pretty=TRUE);quit(status=0)
}
best<-eligible[order(-eligible$cophenetic,eligible$rank),];k<-as.integer(best$rank[1]);fit<-readRDS(file.path(root,'rank_survey',sprintf('K%02d',k),'model_PRIVATE.rds'));labels<-sprintf('CM%02d',seq_len(k));colnames(basis(fit))<-labels;rownames(coef(fit))<-labels
fq<-read.delim(file.path(root,'inputs','frequency_PRIVATE.tsv'),row.names=1,check.names=FALSE,stringsAsFactors=FALSE);pairs<-author$pair_correlation(fq,method='pearson');pairs$author_pval<-pairs$pval;pairs$author_pval_fdr<-pairs$pval_fdr
rr<-pmin(1,pmax(-1,pairs$correlation));df<-ncol(fq)-3;stat<-abs(rr)*sqrt(df/pmax(0,1-rr^2));pairs$raw_Pearson_p<-2*pt(-stat,df=df);pairs$raw_Pearson_BH<-p.adjust(pairs$raw_Pearson_p,'BH');write.table(pairs,file.path(out,'all_pairs_author_and_raw_BH.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
network<-author$cm_network(fit,pairs,top_n=cfg$top_n,corr=cfg$correlation_cutoff,pval_fdr=cfg$FDR)
if(!is.data.frame(network$global$node))network$global$node<-data.frame(subCluster=character(),majorCluster=character())
saveRDS(network,file.path(out,'author_network_PRIVATE.rds'));saveRDS(fit,file.path(out,'selected_model_PRIVATE.rds'))
for(layer in c('global','each'))for(kind in c('node','edge'))write.table(network[[layer]][[kind]],file.path(out,paste0(layer,'_',kind,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
write.table(network$filter,file.path(out,'CM_members_after_author_network_filter.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
W<-basis(fit);H<-coef(fit);S<-smoothing(NMF::fit(fit));A<-W%*%S;stopifnot(max(abs(A%*%H-as.matrix(fitted(fit))))<1e-7)
for(name in c('W','H','S','A'))write.table(get(name),file.path(out,paste0(name,'_PRIVATE.tsv')),sep='\t',quote=FALSE,col.names=NA)
top<-do.call(rbind,lapply(seq_len(k),function(i){ord<-order(-W[,i],rownames(W));ix<-ord[seq_len(min(cfg$top_n,length(ord)))];data.frame(cm=labels[i],subCluster=rownames(W)[ix],majorCluster=fq[rownames(W)[ix],'majorCluster'],weight=W[ix,i],rank_in_factor=seq_along(ix))}));write.table(top,file.path(out,'factor_top15_candidates.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
allnodes<-rownames(W);adj<-matrix(0,length(allnodes),length(allnodes),dimnames=list(allnodes,allnodes));ed<-network$global$edge
if(nrow(ed)){ii<-match(ed$subCluster1,allnodes);jj<-match(ed$subCluster2,allnodes);adj[cbind(ii,jj)]<-1;adj[cbind(jj,ii)]<-1}
set.seed(cfg$permutation_seed);result<-lapply(labels,function(cm){
 nodes<-as.character(network$each$node$subCluster[network$each$node$cm==cm]);n<-length(nodes);edges<-network$each$edge[network$each$edge$cm==cm,,drop=FALSE];m<-nrow(edges);density<-if(n>=2)m/choose(n,2)else 0
 null<-if(n>=2)replicate(cfg$permutations,{ids<-sample.int(length(allnodes),n);sum(adj[ids,ids,drop=FALSE])/2/choose(n,2)})else rep(0,cfg$permutations)
 graph<-if(n)igraph::graph_from_adjacency_matrix(adj[nodes,nodes,drop=FALSE],mode='undirected')else igraph::make_empty_graph()
 data.frame(cm=cm,n_top_candidates=sum(top$cm==cm),n_network_nodes=n,n_network_edges=m,n_parent_types=length(unique(fq[nodes,'majorCluster'])),n_cross_parent_edges=if(m)sum(edges$majorCluster1!=edges$majorCluster2)else 0,n_components=if(n)components(graph)$no else 0,connectivity=density,permutation_P=(1+sum(null>=density))/(cfg$permutations+1),null_mean=mean(null))
});summary<-do.call(rbind,result);summary$permutation_BH_diagnostic<-p.adjust(summary$permutation_P,'BH');summary$multicellular_network<-summary$n_parent_types>=2&summary$n_cross_parent_edges>0;write.table(summary,file.path(out,'CM_network_topology_and_permutation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='AUTHOR_CONFIG_CM_IDENTIFICATION_COMPLETE_PENDING_INDEPENDENT_AUDIT',UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC'),selected_rank=k,paper_admissible_ranks=as.integer(eligible$rank),n_factors=k,n_nonempty_networks=sum(summary$n_network_edges>0),n_multicellular_networks=sum(summary$multicellular_network),n_nodes=nrow(W),n_samples=ncol(H),top_n=cfg$top_n,specificity_cutoff=1-(2*(cfg$top_n-1)-1)/(2*(nrow(W)-1)-1),permutations=cfg$permutations,permutation_definition='Sample the same number of nodes from all measured input nodes; preserve global author-selected graph. Report all factors including empty and single-parent networks.',author_P_values='psych default adjusted upper triangle followed by author BH; raw Pearson/BH separate diagnostics',rank_rule='rho[k-2]<rho[k-1]<rho[k]>rho[k+1]; maximize admissible rho; smaller k only for exact ties',runtime=lapply(c('NMF','dplyr','igraph','psych'),function(x)paste0(x,'=',packageVersion(x))),all_atlas_or_biological_validation_claim=FALSE),file.path(out,'complete.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
