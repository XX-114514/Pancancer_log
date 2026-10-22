# Data-informed secondary family. Native objects, graphs and resamples stay RAM.
LN_spec_table <- function(reg_by_cycle, totals, endpoint) {
  a<-reg_by_cycle[,endpoint];b<-reg_by_cycle[,1]-a
  c<-totals[,endpoint]-a;d<-totals[,1]-a-b-c
  out<-cbind(a,b,c,d);stopifnot(all(out>= -1e-10));pmax(out,0)
}
LN_spec_effects <- function(tab) {
  z<-tab+.5;n<-rowSums(z)
  lor<-log(z[,1]*z[,4]/(z[,2]*z[,3]))
  mh<-log(sum(z[,1]*z[,4]/n)/sum(z[,2]*z[,3]/n))
  c(noncycling=lor[1],cycling=lor[2],combined_MH=mh,cycling_minus_noncycling=lor[2]-lor[1])
}
LN_spec_blocks <- function(xy,anchor,cycle,degree,grid) {
  ranges<-apply(xy,2,range);width<-ranges[2,]-ranges[1,];stopifnot(all(width>0))
  scaled<-sweep(sweep(xy[anchor,,drop=FALSE],2,ranges[1,],"-"),2,width,"/")
  gx<-pmin(grid-1L,floor(scaled[,1]*grid));gy<-pmin(grid-1L,floor(scaled[,2]*grid))
  db<-cut(degree,breaks=c(-1,0,1,3,6,Inf),labels=FALSE)
  split(seq_along(anchor),paste(gx,gy,as.integer(cycle),db,sep=":"))
}
LN_spec_shuffle <- function(reg,blocks) {
  r<-reg
  for(i in blocks)if(length(i)>1L)r[i]<-reg[i[sample.int(length(i))]]
  r
}
LN_spec_cases_matrix <- function(values,meta) {
  groups<-split(seq_len(nrow(meta)),meta$patient)
  out<-lapply(groups,function(i){nodes<-split(i,meta$node[i]);Reduce("+",lapply(nodes,function(j)colMeans(values[j,,drop=FALSE])))/length(nodes)})
  if(!length(out))return(matrix(numeric(),0,ncol(values)))
  z<-do.call(rbind,out);rownames(z)<-names(groups);z
}
LN_spec_null_summary <- function(obs,perms,seed,Bboot=2000L,minimum_cases=5L,infer=TRUE) {
  stopifnot(nrow(perms)==length(obs),all(is.finite(perms)),all(is.finite(obs)))
  n<-length(obs);B<-ncol(perms)
  # Symmetric centering includes the observed column, preserving column exchangeability.
  center<-rowMeans(cbind(obs,perms));corrected<-obs-center
  nt<-if(n)apply(sweep(perms,1,center,"-"),2,median) else rep(NA_real_,B)
  point<-if(n)median(corrected) else NA_real_
  low<-if(n)sum(nt<=point+1e-12) else 0L;high<-if(n)sum(nt>=point-1e-12) else 0L
  variable<-n&&diff(range(nt))>1e-12
  p<-if(infer&&n>=minimum_cases&&variable)min(1,2*(1+min(low,high))/(B+1)) else NA_real_
  ci<-c(NA_real_,NA_real_);indices<-matrix(integer(),0,0)
  if(n&&infer){set.seed(seed);indices<-matrix(sample.int(n,n*Bboot,replace=TRUE),nrow=Bboot);boots<-apply(matrix(corrected[indices],nrow=Bboot),1,median);ci<-as.numeric(quantile(boots,c(.025,.975),names=FALSE,type=7))}
  half<-B%/%2L;split_difference<-if(n&&half>0L){a<-obs-rowMeans(perms[,seq_len(half),drop=FALSE]);b<-obs-rowMeans(perms[,(half+1L):B,drop=FALSE]);abs(median(a)-median(b))}else NA_real_
  list(summary=list(cases=n,median_effect=point,bootstrap_lower=ci[1],bootstrap_upper=ci[2],
    positive_cases=sum(corrected>1e-12),negative_cases=sum(corrected< -1e-12),zero_cases=sum(abs(corrected)<=1e-12),
    permutation_p=p,permutations=B,lower_tail_count=low,upper_tail_count=high,
    source_null_varies=isTRUE(variable),split_stream_centering_difference=split_difference,
    leave_one_case_median_min=if(n>1)min(vapply(seq_len(n),function(i)median(corrected[-i]),numeric(1))) else NA_real_,
    leave_one_case_median_max=if(n>1)max(vapply(seq_len(n),function(i)median(corrected[-i]),numeric(1))) else NA_real_,
    status=if(n<minimum_cases)"INSUFFICIENT_CASES" else if(!variable)"DEGENERATE_CONDITIONAL_NULL" else "ESTIMABLE",
    inference=if(infer)"two-sided source-level conditional randomization; patient bootstrap conditional on fitted labels and null center" else "point and coverage only; no additional hypothesis test"),
    internal=list(corrected_cases=corrected,center=center,null_statistics=nt,bootstrap_indices=indices))
}

measure_LN_regulatory_specificity_RAM <- function(x,phenotype,config) {
  av<-function(v,k)attr(v,k,exact=TRUE);dl<-function(v){if(is.list(v)&&!isS4(v))return(v);attr(v,"listData",exact=TRUE)}
  cd<-dl(av(x,"colData"));xy<-av(x,"spatialCoords");if(is.null(xy))xy<-dl(av(x,"int_colData"))[["spatialCoords"]]
  image<-as.character(cd[[config$image_field]]);patient<-as.character(cd[[config$patient_field]]);node<-as.character(cd[[config$node_field]])
  groups<-split(seq_along(image),image);meta<-data.frame(image=names(groups),patient=vapply(groups,function(i)patient[i[1]],character(1)),node=vapply(groups,function(i)node[i[1]],character(1)),stringsAsFactors=FALSE)
  stopifnot(all(vapply(groups,function(i)length(unique(patient[i]))==1L&&length(unique(node[i]))==1L,logical(1))))
  endpoints<-c("E0_PD1_TIM3","E1_CD8_MKI67","E2_CCR7_CD127");estimands<-c("noncycling","cycling","combined_MH","cycling_minus_noncycling")
  direct<-array(NA_real_,c(length(groups),3,4),dimnames=list(names(groups),endpoints,estimands));internals<-list()
  nulls<-setNames(lapply(config$grids,function(g)array(NA_real_,c(length(groups),config$permutations,3),dimnames=list(names(groups),NULL,endpoints))),as.character(config$grids))
  mobility_records<-list();missing<-list();edge_margin_flags<-matrix(FALSE,length(groups),3);p<-phenotype$internal
  for(imno in seq_along(groups)) {
    if(Sys.time()>as.POSIXct(config$absolute_deadline_UTC,tz="UTC"))stop("FIXED_DEADLINE")
    im<-names(groups)[imno];i<-groups[[imno]];coord<-xy[i,1:2,drop=FALSE];stopifnot(all(is.finite(coord)))
    anchor<-which(p$parents$CD4[i]);Q<-p$parents$CD8a[i];reg<-p$states$R0_FOXP3_CTLA4[i][anchor];cycle<-p$states$R2_general_CD4_MKI67[i][anchor]
    base_support<-im%in%p$testable_images$R0_FOXP3_CTLA4&&im%in%p$testable_images$R2_general_CD4_MKI67
    cell_gate<-base_support&&length(anchor)>=config$minimum_parent&&sum(Q)>=config$minimum_parent&&all(vapply(c(FALSE,TRUE),function(cy)sum(reg&cycle==cy)>=config$minimum_class&&sum(!reg&cycle==cy)>=config$minimum_class,logical(1)))
    reasons<-rep("ESTIMABLE",3);eligible<-vapply(endpoints,function(e)cell_gate&&im%in%p$testable_images[[e]]&&sum(p$states[[e]][i])>=config$minimum_class&&sum(Q&!p$states[[e]][i])>=config$minimum_class,logical(1))
    reasons[!eligible]<-if(!base_support)"REQUIRED_CD4_MARKER_MEASUREMENT" else if(!cell_gate)"DISJOINT_CD4_CLASS_OR_PARENT_COVERAGE" else "CD8_ENDPOINT_MEASUREMENT_OR_COMPLEMENT_COVERAGE"
    missing[[imno]]<-reasons
    if(!any(eligible)){internals[[im]]<-list(eligible=eligible,reason=reasons);next}
    stopifnot(length(i)>config$k);nn<-RANN::nn2(coord,query=coord,k=config$k+1L,eps=0)
    nb<-matrix(NA_integer_,length(i),config$k);distances<-nb*NA_real_
    for(j in seq_along(i)){take<-which(nn$nn.idx[j,]!=j)[seq_len(config$k)];nb[j,]<-nn$nn.idx[j,take];distances[j,]<-nn$nn.dists[j,take]}
    D<-vapply(c("CD8_parent",endpoints),function(e){mask<-if(e=="CD8_parent")Q else p$states[[e]][i];rowSums(matrix(mask[nb[anchor,,drop=FALSE]],nrow=length(anchor)))},numeric(length(anchor)))
    bycycle<-cbind(!cycle,cycle);totals<-crossprod(bycycle,D);regtotals<-crossprod(bycycle*as.numeric(reg),D)
    tables<-lapply(seq_len(3),function(e)LN_spec_table(regtotals,totals,e+1L))
    for(e in which(eligible)){
      direct[imno,e,]<-LN_spec_effects(tables[[e]])
      edge_margin_flags[imno,e]<-any(apply(tables[[e]],1,function(z)any(c(z[1]+z[2],z[3]+z[4],z[1]+z[3],z[2]+z[4])==0)))
    }
    blocksets<-list();permutation_seeds<-list()
    for(grid in config$grids) {
      blocks<-LN_spec_blocks(coord,anchor,cycle,D[,1],grid);mixed<-blocks[vapply(blocks,function(j)any(reg[j])&&any(!reg[j]),logical(1))]
      movable<-unique(unlist(mixed,use.names=FALSE));movable<-movable[D[movable,1]>0]
      can_shuffle<-length(mixed)>=config$minimum_mixed_blocks&&sum(reg[movable])>=config$minimum_movable_class&&sum(!reg[movable])>=config$minimum_movable_class
      mobility_records[[length(mobility_records)+1L]]<-list(grid=grid,mixed_blocks=length(mixed),edge_relevant_movable_reg=sum(reg[movable]),edge_relevant_movable_other=sum(!reg[movable]),eligible=can_shuffle)
      seed<-config$seed+imno*100L+grid;blocksets[[as.character(grid)]]<-blocks;permutation_seeds[[as.character(grid)]]<-seed
      if(can_shuffle) {
        set.seed(seed)
        for(b in seq_len(config$permutations)) {
          shuffled<-LN_spec_shuffle(reg,blocks);r<-crossprod(bycycle*as.numeric(shuffled),D)
          for(e in which(eligible))nulls[[as.character(grid)]][imno,b,e]<-LN_spec_effects(LN_spec_table(r,totals,e+1L))["combined_MH"]
        }
      }
    }
    internals[[im]]<-list(eligible=eligible,cell_rows=i,neighbors=nb,distances=distances,anchor=anchor,reg=reg,cycle=cycle,Q=Q,D=D,totals=totals,tables=tables,blocks=blocksets,seeds=permutation_seeds)
  }
  rows<-list();tests<-list();slot<-0L
  for(e in seq_along(endpoints))for(est in seq_along(estimands)) {
    slot<-slot+1L;ok<-is.finite(direct[,e,est]);case<-LN_spec_cases_matrix(matrix(direct[ok,e,est],ncol=1),meta[ok,,drop=FALSE])[,1]
    z<-LN_exact_direction_summary_RAM(case,config$seed+10000L+slot,config$bootstrap_replicates,config$minimum_cases)
    rows[[slot]]<-c(list(endpoint=endpoints[e],estimand=estimands[est],eligible_images=sum(ok),p_method="exact two-sided patient direction test"),z$summary,list(primary_p=z$summary$exact_sign_p))
    tests[[slot]]<-c(list(kind="direct",case_values=case),z$internal)
  }
  sensitivities<-list();null_tests<-list()
  for(grid in config$grids)for(e in seq_along(endpoints)) {
    vals<-matrix(nulls[[as.character(grid)]][,,e],nrow=nrow(meta),ncol=config$permutations);ok<-is.finite(direct[,e,"combined_MH"])&rowSums(is.finite(vals))==config$permutations
    cases<-LN_spec_cases_matrix(cbind(observed=direct[ok,e,"combined_MH"],vals[ok,,drop=FALSE]),meta[ok,,drop=FALSE]);infer<-grid==config$primary_grid
    z<-LN_spec_null_summary(cases[,1],cases[,-1,drop=FALSE],config$seed+20000L+e,config$bootstrap_replicates,config$minimum_cases,infer)
    item<-c(list(endpoint=endpoints[e],estimand="local_composition_degree_and_cycle_conditional_null",grid=grid,eligible_images=sum(ok)),z$summary)
    if(infer){slot<-slot+1L;rows[[slot]]<-c(item,list(primary_p=z$summary$permutation_p));tests[[slot]]<-c(list(kind="null",case_matrix=cases),z$internal)}else sensitivities[[length(sensitivities)+1L]]<-item
    null_tests[[paste(grid,e,sep="_")]]<-list(case_matrix=cases,summary=z)
  }
  stopifnot(length(rows)==15L);pv<-vapply(rows,function(r)r$primary_p,numeric(1));pad<-pv;pad[!is.finite(pad)]<-1
  for(j in seq_along(rows)){rows[[j]]$BH<-if(is.finite(pv[j]))p.adjust(pad,"BH")[j] else NA_real_;rows[[j]]$BY<-if(is.finite(pv[j]))p.adjust(pad,"BY")[j] else NA_real_}
  reason_totals<-as.list(table(unlist(missing,use.names=FALSE)))
  mobility_summary<-lapply(config$grids,function(grid){z<-Filter(function(x)x$grid==grid,mobility_records);list(grid=grid,images_evaluated=length(z),images_with_mobility=sum(vapply(z,function(x)x$eligible,logical(1))),mixed_blocks=sum(vapply(z,function(x)x$mixed_blocks,numeric(1))),edge_relevant_movable_reg=sum(vapply(z,function(x)x$edge_relevant_movable_reg,numeric(1))),edge_relevant_movable_other=sum(vapply(z,function(x)x$edge_relevant_movable_other,numeric(1))))})
  edge_quality<-lapply(seq_along(endpoints),function(e)list(endpoint=endpoints[e],eligible_images=sum(is.finite(direct[,e,1])),images_with_any_zero_edge_margin_stratum=sum(edge_margin_flags[,e])))
  list(summary=list(status="COMPLETE_DATA_INFORMED_SECONDARY_SPECIFICITY_PENDING_AUDIT",primary_rows=rows,sensitivities=sensitivities,coverage_reasons=reason_totals,mobility=mobility_summary,edge_margin_quality=edge_quality,
    planned_primary_tests=15,states_frozen_from_previous_family=TRUE,disjoint_CD4_classes_and_cycle_strata=TRUE,
    local_null_preserves_parent_geometry_receiver_states_and_block_regulatory_counts=TRUE,
    local_null_does_not_fully_preserve_fine_scale_CD4_autocorrelation=TRUE,
    same_source_reanalysis_not_independent_cohort_or_direct_communication=TRUE),
    internal=list(images=internals,meta=meta,direct=direct,nulls=nulls,tests=tests,null_tests=null_tests,phenotype=phenotype))
}
