# Source-aggregate extension; masks, identifiers, graphs and resamples stay RAM.
# Requires unchanged LN table/null helpers, exact-direction helper and unit_key_RAM.
ETDLN_case_ROI_means <- function(values,meta) {
  groups<-split(seq_len(nrow(meta)),meta$case)
  if(!length(groups))return(matrix(numeric(),0,ncol(values)))
  out<-do.call(rbind,lapply(groups,function(i)colMeans(values[i,,drop=FALSE])))
  rownames(out)<-names(groups);out
}

measure_external_TDLN_fixed_spatial_RAM <- function(xy,case,roi,phenotype,config) {
  stopifnot(is.matrix(xy),is.numeric(xy),ncol(xy)==2L,nrow(xy)==length(case),length(case)==length(roi))
  case<-as.character(case);roi<-as.character(roi);image<-unit_key_RAM(list(case,roi))
  groups<-split(seq_along(image),image)
  meta<-data.frame(image=names(groups),case=vapply(groups,function(i)case[i[1]],character(1)),stringsAsFactors=FALSE)
  stopifnot(identical(image,phenotype$internal$image),all(vapply(groups,function(i)length(unique(case[i]))==1L,logical(1))))
  endpoints<-c("E0_PD1_TIM3","E1_CD8_MKI67","E2_CCR7_CD127");estimands<-c("noncycling","cycling","combined_MH","cycling_minus_noncycling")
  available<-setNames(vapply(endpoints,function(e)isTRUE(phenotype$internal$measurement_available[[e]]),logical(1)),endpoints)
  direct<-array(NA_real_,c(length(groups),3,4),dimnames=list(names(groups),endpoints,estimands));internals<-list()
  nulls<-setNames(lapply(config$grids,function(g)array(NA_real_,c(length(groups),config$permutations,3),dimnames=list(names(groups),NULL,endpoints))),as.character(config$grids))
  mobility_records<-list();missing<-list();edge_margin_flags<-matrix(FALSE,length(groups),3);p<-phenotype$internal
  for(imno in seq_along(groups)) {
    if(Sys.time()>as.POSIXct(config$absolute_deadline_UTC,tz="UTC"))stop("FIXED_DEADLINE")
    im<-names(groups)[imno];i<-groups[[imno]];coord<-xy[i,1:2,drop=FALSE];stopifnot(all(is.finite(coord)))
    if(nrow(coord)<2L||any(apply(coord,2,function(z)diff(range(z)))<=0)) {
      reasons<-ifelse(available,"DEGENERATE_IMAGE_GEOMETRY","STRUCTURAL_MARKER_ABSENCE")
      missing[[imno]]<-reasons;internals[[im]]<-list(eligible=rep(FALSE,3),reason=reasons);next
    }
    anchor<-which(p$parents$CD4[i]);Q<-p$parents$CD8a[i];reg<-p$states$R0_FOXP3_CTLA4[i][anchor];cycle<-p$states$R2_general_CD4_MKI67[i][anchor]
    base_support<-im%in%p$testable_images$R0_FOXP3_CTLA4&&im%in%p$testable_images$R2_general_CD4_MKI67
    cell_gate<-base_support&&length(anchor)>=config$minimum_parent&&sum(Q)>=config$minimum_parent&&all(vapply(c(FALSE,TRUE),function(cy)sum(reg&cycle==cy)>=config$minimum_class&&sum(!reg&cycle==cy)>=config$minimum_class,logical(1)))
    reasons<-rep("ESTIMABLE",3);eligible<-vapply(endpoints,function(e)available[[e]]&&cell_gate&&im%in%p$testable_images[[e]]&&sum(p$states[[e]][i])>=config$minimum_class&&sum(Q&!p$states[[e]][i])>=config$minimum_class,logical(1))
    reasons[available&!eligible]<-if(!base_support)"REQUIRED_CD4_MARKER_MEASUREMENT" else if(!cell_gate)"DISJOINT_CD4_CLASS_OR_PARENT_COVERAGE" else "CD8_ENDPOINT_MEASUREMENT_OR_COMPLEMENT_COVERAGE"
    reasons[!available]<-"STRUCTURAL_MARKER_ABSENCE"
    missing[[imno]]<-reasons
    if(!any(eligible)){internals[[im]]<-list(eligible=eligible,reason=reasons);next}
    stopifnot(length(i)>config$k);nn<-RANN::nn2(coord,query=coord,k=config$k+1L,eps=0)
    nb<-matrix(NA_integer_,length(i),config$k);distances<-nb*NA_real_
    for(j in seq_along(i)){take<-which(nn$nn.idx[j,]!=j)[seq_len(config$k)];nb[j,]<-nn$nn.idx[j,take];distances[j,]<-nn$nn.dists[j,take]}
    D<-matrix(NA_real_,length(anchor),4,dimnames=list(NULL,c("CD8_parent",endpoints)))
    D[,1]<-rowSums(matrix(Q[nb[anchor,,drop=FALSE]],nrow=length(anchor)))
    for(e in which(eligible)) {
      mask<-p$states[[endpoints[e]]][i];stopifnot(!anyNA(mask))
      D[,e+1L]<-rowSums(matrix(mask[nb[anchor,,drop=FALSE]],nrow=length(anchor)))
    }
    bycycle<-cbind(!cycle,cycle);totals<-crossprod(bycycle,D);regtotals<-crossprod(bycycle*as.numeric(reg),D)
    tables<-rep(list(NULL),3)
    for(e in which(eligible))tables[[e]]<-LN_spec_table(regtotals,totals,e+1L)
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
    slot<-slot+1L;ok<-is.finite(direct[,e,est]);case<-ETDLN_case_ROI_means(matrix(direct[ok,e,est],ncol=1),meta[ok,,drop=FALSE])[,1]
    z<-LN_exact_direction_summary_RAM(case,config$seed+10000L+slot,config$bootstrap_replicates,config$minimum_cases)
    rows[[slot]]<-c(list(endpoint=endpoints[e],estimand=estimands[est],eligible_images=sum(ok),p_method="exact two-sided registered-case direction test"),z$summary,list(primary_p=z$summary$exact_sign_p))
    tests[[slot]]<-c(list(kind="direct",case_values=case),z$internal)
  }
  sensitivities<-list();null_tests<-list()
  for(grid in config$grids)for(e in seq_along(endpoints)) {
    vals<-matrix(nulls[[as.character(grid)]][,,e],nrow=nrow(meta),ncol=config$permutations);ok<-is.finite(direct[,e,"combined_MH"])&rowSums(is.finite(vals))==config$permutations
    cases<-ETDLN_case_ROI_means(cbind(observed=direct[ok,e,"combined_MH"],vals[ok,,drop=FALSE]),meta[ok,,drop=FALSE]);infer<-grid==config$primary_grid
    z<-LN_spec_null_summary(cases[,1],cases[,-1,drop=FALSE],config$seed+20000L+e,config$bootstrap_replicates,config$minimum_cases,infer)
    item<-c(list(endpoint=endpoints[e],estimand="local_composition_degree_and_cycle_conditional_null",grid=grid,eligible_images=sum(ok)),z$summary)
    if(infer){slot<-slot+1L;rows[[slot]]<-c(item,list(primary_p=z$summary$permutation_p));tests[[slot]]<-c(list(kind="null",case_matrix=cases),z$internal)}else sensitivities[[length(sensitivities)+1L]]<-item
    null_tests[[paste(grid,e,sep="_")]]<-list(case_matrix=cases,summary=z)
  }
  stopifnot(length(rows)==15L);pv<-vapply(rows,function(r)r$primary_p,numeric(1));pad<-pv;pad[!is.finite(pad)]<-1
  for(j in seq_along(rows)){rows[[j]]$BH<-if(is.finite(pv[j]))p.adjust(pad,"BH")[j] else NA_real_;rows[[j]]$BY<-if(is.finite(pv[j]))p.adjust(pad,"BY")[j] else NA_real_}
  for(j in seq_along(rows)) {
    if(!available[[rows[[j]]$endpoint]])rows[[j]]$status<-"STRUCTURAL_MARKER_ABSENCE"
    if(!is.null(rows[[j]]$inference))rows[[j]]$inference<-"source-level conditional randomization; registered-case bootstrap conditional on fixed labels and null center"
  }
  reason_totals<-as.list(table(unlist(missing,use.names=FALSE)))
  mobility_summary<-lapply(config$grids,function(grid){z<-Filter(function(x)x$grid==grid,mobility_records);list(grid=grid,images_evaluated=length(z),images_with_mobility=sum(vapply(z,function(x)x$eligible,logical(1))),mixed_blocks=sum(vapply(z,function(x)x$mixed_blocks,numeric(1))),edge_relevant_movable_reg=sum(vapply(z,function(x)x$edge_relevant_movable_reg,numeric(1))),edge_relevant_movable_other=sum(vapply(z,function(x)x$edge_relevant_movable_other,numeric(1))))})
  edge_quality<-lapply(seq_along(endpoints),function(e)list(endpoint=endpoints[e],eligible_images=sum(is.finite(direct[,e,1])),images_with_any_zero_edge_margin_stratum=sum(edge_margin_flags[,e])))
  list(summary=list(status="COMPLETE_EXTERNAL_FIXED_STATE_EXTENSION_PENDING_AUDIT",primary_rows=rows,sensitivities=sensitivities,coverage_reasons=reason_totals,mobility=mobility_summary,edge_margin_quality=edge_quality,
    planned_primary_tests=15,relative_state_cutoffs_frozen_from_previous_construct=TRUE,disjoint_CD4_classes_and_cycle_strata=TRUE,
    local_null_preserves_parent_geometry_receiver_states_and_block_regulatory_counts=TRUE,
    local_null_does_not_fully_preserve_fine_scale_CD4_autocorrelation=TRUE,
    new_source_relative_construct_not_exact_RNA_CM_projection_or_direct_communication=TRUE,
    case_ROI_aggregation="ROI equal within registered case; no inferred lymph-node level",
    literal_case_keys_not_accepted_independent_patients=TRUE,
    structural_missing_endpoints_are_NA_not_negative=TRUE),
    internal=list(images=internals,meta=meta,direct=direct,nulls=nulls,tests=tests,null_tests=null_tests,phenotype=phenotype))
}
