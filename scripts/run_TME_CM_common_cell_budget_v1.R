# Author minmax is re-estimated in each cell draw, with identity/rank fixed.
args<-commandArgs(TRUE);stopifnot(length(args)==1)
cfg<-jsonlite::fromJSON(args[1],simplifyVector=FALSE);suppressPackageStartupMessages(library(NMF))
defs<-parse(cfg$math_source);api<-new.env(parent=.GlobalEnv)
for(n in c('cosine','kl_project','fit_parts')){i<-which(vapply(defs,function(e)is.call(e)&&length(e)>=3&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name(n)),logical(1)));stopifnot(length(i)==1);eval(defs[[i]],api)}
idx<-jsonlite::fromJSON(cfg$input_index,simplifyVector=FALSE);results<-list();nr<-0L
normalise<-function(x){span<-apply(x,1,max)-apply(x,1,min);z<-x;z[,]<-0;v<-span>0;if(any(v))z[v,]<-sweep(sweep(x[v,,drop=FALSE],1,apply(x[v,,drop=FALSE],1,min),'-'),1,span[v],'/');list(X=z,n_constant_rows=sum(!v))}
toy<-rbind(c(1,1,1),c(0,.5,1));tt<-normalise(toy);stopifnot(tt$n_constant_rows==1,all(tt$X[1,]==0),identical(as.numeric(tt$X[2,]),c(0,.5,1)))
for(s in idx$schemes){
 root<-s$root;sc<-jsonlite::fromJSON(file.path(root,'config/spec_v1.json'));input<-jsonlite::fromJSON(sc$input_acceptance$path);fq<-read.delim(input$frequency$path,row.names=1,check.names=FALSE,stringsAsFactors=FALSE);meta<-read.delim(input$metadata$path,stringsAsFactors=FALSE);stopifnot(identical(colnames(fq)[-1],meta$sampleID))
 raw<-as.matrix(fq[,-1,drop=FALSE]);baseX<-normalise(raw)$X;D<-file.path(root,'outputs',paste0('seed_',cfg$baseline_seed));base<-list(W=as.matrix(read.delim(file.path(D,'private_W_v1.tsv'),row.names=1,check.names=FALSE)),A=as.matrix(read.delim(file.path(D,'private_A_v1.tsv'),row.names=1,check.names=FALSE)));Hbase<-api$kl_project(base$A,baseX)$H
 weights<-1/as.numeric(table(meta$patient_candidate_block)[meta$patient_candidate_block]);weights<-weights/mean(weights);k<-as.integer(cfg$rank)
 for(rep in idx$repeats){
  file<-Filter(function(z)z$scheme==s$id,rep$schemes)[[1]]$frequency$path;sample<-read.delim(gzfile(file),row.names=1,check.names=FALSE,stringsAsFactors=FALSE);stopifnot(identical(rownames(sample),rownames(fq)),identical(colnames(sample),colnames(fq)))
  sampled<-as.matrix(sample[,-1,drop=FALSE]);norm<-normalise(sampled);X<-norm$X
  if(sum(rowSums(X)>0)<=k){nr<-nr+1L;results[[nr]]<-data.frame(scheme=s$id,repeat_index=rep[['repeat']],status='insufficient_variable_nodes_after_cell_draw',n_constant_rows=norm$n_constant_rows);next}
  p<-api$fit_parts(X,weights,k,as.integer(cfg$seed+rep[['repeat']]),as.integer(cfg$nrun),as.integer(cfg$max_iter));H<-api$kl_project(p$A,X)$H
  for(kind in c('W','A')){
   co<-api$cosine(base[[kind]],p[[kind]]);benefit<-matrix(0,2*k,2*k);benefit[seq_len(k),seq_len(k)]<-co-cfg$minimum_cosine;benefit<-benefit-min(benefit);mt<-as.integer(clue::solve_LSAP(benefit,maximum=TRUE))
   for(i in seq_len(k)){
    j<-mt[i];ok<-j<=k&&co[i,j]>=cfg$minimum_cosine;rho<-if(ok&&sd(Hbase[i,])>0&&sd(H[j,])>0)cor(Hbase[i,],H[j,],method='spearman')else NA_real_
    nr<-nr+1L;results[[nr]]<-data.frame(scheme=s$id,repeat_index=rep[['repeat']],status=if(ok)'matched'else'unmatched',basis=kind,baseline_module=i,refit_module=if(ok)j else NA_integer_,loading_cosine=if(ok)co[i,j]else NA_real_,same_material_activity_Spearman=rho,n_constant_rows=norm$n_constant_rows,n_analysis_units=ncol(X),cells_per_parent_unit=idx$cells_per_parent_unit,normalization='author_minmax_each_draw_fixed_identity_rank')
   }
  }
  jsonlite::write_json(list(scheme=s$id,n_repeats_completed=rep[['repeat']]+1,n_repeats_planned=idx$n_repeats),file.path(cfg$output,paste0('checkpoint_',s$id,'.json')),auto_unbox=TRUE,pretty=TRUE)
 }
}
keys<-unique(unlist(lapply(results,names)));results<-lapply(results,function(x){for(n in setdiff(keys,names(x)))x[[n]]<-NA;x[,keys,drop=FALSE]});out<-do.call(rbind,results);write.table(out,file.path(cfg$output,'private_CM_same_cell_budget_stability_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='SAME_CELL_DRAWS_ALL_SCHEMES_COVARNET_ROBUSTNESS_COMPLETE',n_schemes=length(idx$schemes),n_repeats=idx$n_repeats,n_analysis_units=idx$n_analysis_units,cells_per_parent_unit=idx$cells_per_parent_unit,normalization='author_minmax_reestimated_per_draw',same_cells_across_all_annotation_schemes=TRUE,rank_and_node_vocabulary_frozen=TRUE,nonvariable_sampled_rows_explicit_not_imputed=TRUE,not_independent_patient_validation=TRUE),file.path(cfg$output,'CM_same_cell_budget_report_v1.json'),auto_unbox=TRUE,pretty=TRUE)
