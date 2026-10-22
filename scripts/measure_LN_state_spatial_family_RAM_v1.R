# Patient-level source summaries only; every graph, membership and case vector stays RAM.
LN_spatial_image_metrics_RAM <- function(neighbors, states, parents, testable, minimum_parent=20L, minimum_state=10L) {
  n<-nrow(neighbors);k<-ncol(neighbors);stopifnot(all(neighbors>=1L&neighbors<=n),all(neighbors!=row(neighbors)))
  src<-rep(seq_len(n),each=k);dst<-as.vector(t(neighbors));rows<-list()
  rnames<-c("R0_FOXP3_CTLA4","R1_FOXP3_CTLA4_MKI67","R2_general_CD4_MKI67");enames<-c("E0_PD1_TIM3","E1_CD8_MKI67","E2_CCR7_CD127")
  for(r in rnames)for(e in enames) {
    A<-states[[r]];B<-states[[e]];P<-parents$CD4;Q<-parents$CD8a;stopifnot(all(!A|P),all(!B|Q))
    eligible<-sum(P)>=minimum_parent&&sum(Q)>=minimum_parent&&isTRUE(testable[[r]])&&isTRUE(testable[[e]])&&sum(A)>=minimum_state&&sum(B)>=minimum_state
    if(!eligible){reason<-character();if(sum(P)<minimum_parent||sum(Q)<minimum_parent)reason<-c(reason,"PARENT_COVERAGE");if(!isTRUE(testable[[r]])||!isTRUE(testable[[e]]))reason<-c(reason,"STATE_MARKER_MEASUREMENT_UNINFORMATIVE");if(sum(A)<minimum_state||sum(B)<minimum_state)reason<-c(reason,"STATE_CELL_COVERAGE")
      rows[[length(rows)+1L]]<-list(pair=paste(r,e,sep="__"),whole=NA_real_,conditional=NA_real_,reason=paste(reason,collapse=";"),counts=NULL);next}
    obs<-sum(A[src]&B[dst]);expected<-sum(A)*k*sum(B)/(n-1)
    take<-P[src]&Q[dst];as<-A[src[take]];bs<-B[dst[take]]
    tab<-c(a=sum(as&bs),b=sum(as&!bs),c=sum(!as&bs),d=sum(!as&!bs))
    margins_ok<-all(c(tab[1]+tab[2],tab[3]+tab[4],tab[1]+tab[3],tab[2]+tab[4])>0)
    conditional<-if(margins_ok)log((tab[1]+.5)*(tab[4]+.5)/((tab[2]+.5)*(tab[3]+.5))) else NA_real_
    rows[[length(rows)+1L]]<-list(pair=paste(r,e,sep="__"),whole=log((obs+.5)/(expected+.5)),conditional=unname(conditional),
      reason=if(margins_ok)"ESTIMABLE" else "ZERO_CROSS_PARENT_EDGE_MARGIN",counts=c(observed=obs,expected=expected,tab))
  }
  rows
}

LN_patient_effects_RAM <- function(records,meta,metric,station=NULL) {
  v<-vapply(records,function(z)z[[metric]],numeric(1));keep<-is.finite(v)
  if(!is.null(station))keep<-keep&!is.na(meta$station)&meta$station==station
  if(!any(keep))return(setNames(numeric(),character()))
  node_key<-paste0(nchar(meta$patient,type="bytes"),":",meta$patient,nchar(meta$node,type="bytes"),":",meta$node);groups<-split(which(keep),node_key[keep]);node_means<-vapply(groups,function(i)mean(v[i]),numeric(1))
  pat<-vapply(groups,function(i)meta$patient[i[1]],character(1))
  vapply(split(node_means,pat),mean,numeric(1))
}

LN_exact_direction_summary_RAM <- function(x,seed,B=2000L,minimum_cases=5L) {
  x<-x[is.finite(x)];n<-length(x);positive<-sum(x>1e-12);negative<-sum(x< -1e-12);zero<-n-positive-negative
  direction_n<-positive+negative;p<-if(n>=minimum_cases){if(direction_n)min(1,2*pbinom(min(positive,negative),direction_n,.5)) else 1} else NA_real_
  if(n){set.seed(seed);indices<-matrix(sample.int(n,n*B,replace=TRUE),nrow=B);boot<-apply(matrix(x[indices],nrow=B),1,median);ci<-as.numeric(quantile(boot,c(.025,.975),type=7,names=FALSE))}else{indices<-matrix(integer(),0,0);ci<-c(NA_real_,NA_real_)}
  list(summary=list(cases=n,positive_cases=positive,negative_cases=negative,zero_cases=zero,median_effect=if(n)median(x) else NA_real_,
    bootstrap_lower=ci[1],bootstrap_upper=ci[2],exact_sign_p=p,minimum_case_gate=minimum_cases,
    leave_one_case_median_min=if(n>1)min(vapply(seq_len(n),function(i)median(x[-i]),numeric(1))) else NA_real_,
    leave_one_case_median_max=if(n>1)max(vapply(seq_len(n),function(i)median(x[-i]),numeric(1))) else NA_real_,
    status=if(n>=minimum_cases)"ESTIMABLE" else "INSUFFICIENT_CASES"),
    internal=list(effects=x,bootstrap_indices=indices))
}

measure_LN_state_spatial_family_RAM <- function(x,phenotypes,config) {
  av<-function(v,k)attr(v,k,exact=TRUE);dl<-function(v){if(is.list(v)&&!isS4(v))return(v);attr(v,"listData",exact=TRUE)}
  cd<-dl(av(x,"colData"));coords<-av(x,"spatialCoords");if(is.null(coords))coords<-dl(av(x,"int_colData"))[["spatialCoords"]]
  image<-as.character(cd[[config$image_field]]);patient<-as.character(cd[[config$patient_field]]);node<-as.character(cd[[config$node_field]]);station<-toupper(trimws(as.character(cd[[config$station_field]])))
  groups<-split(seq_along(image),image);meta<-data.frame(image=names(groups),patient=vapply(groups,function(i)patient[i[1]],character(1)),node=vapply(groups,function(i)node[i[1]],character(1)),station=vapply(groups,function(i)station[i[1]],character(1)),stringsAsFactors=FALSE)
  stopifnot(all(vapply(groups,function(i)length(unique(patient[i]))==1L&&length(unique(node[i]))==1L&&length(unique(station[i]))==1L,logical(1))))
  graphs<-list();metrics<-list()
  for(im in names(groups)) {
    if(Sys.time()>as.POSIXct(config$absolute_deadline_UTC,tz="UTC"))stop("FIXED_DEADLINE")
    i<-groups[[im]];xy<-coords[i,1:2,drop=FALSE];stopifnot(all(is.finite(xy)))
    has_parents<-any(vapply(phenotypes,function(p)sum(p$internal$parents$CD4[i])>=config$minimum_parent_cells&&sum(p$internal$parents$CD8a[i])>=config$minimum_parent_cells,logical(1)))
    if(!has_parents) {
      graphs[[im]]<-list(cell_rows=i,neighbors=NULL,distances=NULL,reason="INSUFFICIENT_PARENT_COVERAGE")
      for(mode in names(phenotypes))for(k in unlist(config$k_values)) {
        key<-paste(mode,k,sep="_");if(is.null(metrics[[key]]))metrics[[key]]<-list()
        metrics[[key]][[im]]<-unlist(lapply(c("R0_FOXP3_CTLA4","R1_FOXP3_CTLA4_MKI67","R2_general_CD4_MKI67"),function(r)lapply(c("E0_PD1_TIM3","E1_CD8_MKI67","E2_CCR7_CD127"),function(e)list(pair=paste(r,e,sep="__"),whole=NA_real_,conditional=NA_real_,reason="PARENT_COVERAGE",counts=NULL))),recursive=FALSE)
      }
      next
    }
    stopifnot(length(i)>max(unlist(config$k_values)))
    nn<-RANN::nn2(xy,query=xy,k=max(unlist(config$k_values))+1L,eps=0)
    ids<-matrix(NA_integer_,length(i),max(unlist(config$k_values)))
    distances<-matrix(NA_real_,nrow(ids),ncol(ids))
    for(j in seq_along(i)){take<-which(nn$nn.idx[j,]!=j)[seq_len(ncol(ids))];ids[j,]<-nn$nn.idx[j,take];distances[j,]<-nn$nn.dists[j,take]}
    graphs[[im]]<-list(cell_rows=i,neighbors=ids,distances=distances)
    for(mode in names(phenotypes))for(k in unlist(config$k_values)) {
      p<-phenotypes[[mode]]$internal;st<-lapply(p$states,function(z)z[i]);par<-lapply(p$parents,function(z)z[i]);testable<-lapply(p$testable_images,function(z)im%in%z)
      key<-paste(mode,k,sep="_");if(is.null(metrics[[key]]))metrics[[key]]<-list()
      metrics[[key]][[im]]<-LN_spatial_image_metrics_RAM(ids[,seq_len(k),drop=FALSE],st,par,testable,config$minimum_parent_cells,config$minimum_state_cells)
    }
  }
  pair_names<-vapply(metrics[[paste("primary",config$primary_k,sep="_")]][[1]],function(z)z$pair,character(1))
  rows<-list();internal_tests<-list();slot<-0L;sensitivities<-list()
  for(pair in pair_names) {
    get_records<-function(mode,k)lapply(metrics[[paste(mode,k,sep="_")]][meta$image],function(z)z[[match(pair,vapply(z,function(v)v$pair,character(1)))]])
    rec<-get_records("primary",config$primary_k)
    whole<-LN_patient_effects_RAM(rec,meta,"whole");conditional<-LN_patient_effects_RAM(rec,meta,"conditional")
    n1<-LN_patient_effects_RAM(rec,meta,"conditional","N1");n2<-LN_patient_effects_RAM(rec,meta,"conditional","N2");common<-intersect(names(n1),names(n2));paired<-n1[common]-n2[common]
    for(estimand in c("whole","conditional","N1_minus_N2")) {
      slot<-slot+1L;effects<-switch(estimand,whole=whole,conditional=conditional,N1_minus_N2=paired)
      t<-LN_exact_direction_summary_RAM(effects,config$seed+slot,config$bootstrap_replicates,config$minimum_cases)
      rows[[slot]]<-c(list(pair=pair,estimand=estimand,eligible_images=sum(vapply(rec,function(z)is.finite(z[[if(estimand=="whole")"whole" else "conditional"]]),logical(1)) & if(estimand=="N1_minus_N2")meta$patient%in%common else TRUE)),t$summary)
      internal_tests[[slot]]<-c(list(pair=pair,estimand=estimand),t$internal)
    }
    for(k in unlist(config$k_values))for(mode in names(phenotypes))for(metric in c("whole","conditional")) {
      rr<-get_records(mode,k);effects<-LN_patient_effects_RAM(rr,meta,metric)
      primary_records<-get_records("primary",config$primary_k);joint<-vapply(seq_along(rr),function(i)is.finite(rr[[i]][[metric]])&&is.finite(primary_records[[i]][[metric]]),logical(1))
      ref<-LN_patient_effects_RAM(primary_records[joint],meta[joint,,drop=FALSE],metric);matched<-LN_patient_effects_RAM(rr[joint],meta[joint,,drop=FALSE],metric);common<-intersect(names(ref),names(matched))
      sensitivities[[length(sensitivities)+1L]]<-list(pair=pair,mode=mode,k=k,estimand=metric,cases=length(effects),median_effect=if(length(effects))median(effects) else NA_real_,
        common_cases=length(common),primary_median_on_common_cases=if(length(common))median(ref[common]) else NA_real_,common_eligible_images=sum(joint),sensitivity_median_on_common_cases=if(length(common))median(matched[common]) else NA_real_,no_additional_p_test=TRUE)
    }
  }
  p<-vapply(rows,function(z)z$exact_sign_p,numeric(1));pad<-p;pad[!is.finite(pad)]<-1
  bh<-p.adjust(pad,"BH");by<-p.adjust(pad,"BY")
  for(j in seq_along(rows)){rows[[j]]$BH<-if(is.finite(p[j]))bh[j] else NA_real_;rows[[j]]$BY<-if(is.finite(p[j]))by[j] else NA_real_}
  list(summary=list(status="COMPLETE_FROZEN_SPATIAL_FAMILY_PENDING_INDEPENDENT_AUDIT",primary_rows=rows,sensitivities=sensitivities,
    fixed_planned_primary_tests=length(rows),all_planned_targets_retained=TRUE,
    test_unit="author-registered patient; equal image then equal node aggregation",no_cell_or_edge_level_p_values=TRUE,
    exact_direction_test_not_causal_or_independent_cohort_validation=TRUE,clinical_stage_or_ICI_outcomes_not_tested=TRUE),
    internal=list(graphs=graphs,metrics=metrics,meta=meta,pairs=pair_names,tests=internal_tests,phenotypes=phenotypes))
}
