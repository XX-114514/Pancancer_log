# Independent numeric rank/OLS, sparse graph counts, patient aggregation and statistics.
LN_audit_rank_RAM <- function(z) {
  o<-order(z);r<-rle(z[o]);mid<-cumsum(r$lengths)-(r$lengths-1)/2;out<-numeric(length(z));out[o]<-rep(mid,r$lengths);out
}
LN_audit_residual_RAM <- function(y,x) {xc<-x-mean(x);yc<-y-mean(y);yc-sum(xc*yc)/sum(xc^2)*xc}
LN_audit_median_RAM <- function(z) {z<-sort(z);n<-length(z);if(!n)NA_real_ else (z[(n+1L)%/%2L]+z[(n+2L)%/%2L])/2}
LN_audit_cases_RAM <- function(records,meta,metric,station=NULL) {
  out<-numeric()
  for(p in unique(meta$patient)) {
    values<-numeric()
    for(node in unique(meta$node[meta$patient==p])) {
      i<-which(meta$patient==p&meta$node==node)
      if(!is.null(station))i<-i[!is.na(meta$station[i])&meta$station[i]==station]
      z<-vapply(records[i],function(a)a[[metric]],numeric(1));z<-z[is.finite(z)]
      if(length(z))values<-c(values,sum(z)/length(z))
    }
    if(length(values))out[p]<-sum(values)/length(values)
  }
  out
}
LN_audit_test_RAM <- function(x,indices,minimum_cases=5L) {
  n<-length(x);pos<-sum(x>1e-12);neg<-sum(x< -1e-12);d<-pos+neg;p<-NA_real_
  if(n>=minimum_cases&&d>0L) {
    patterns<-0:(2^d-1);count<-integer(length(patterns))
    for(j in seq_len(d))count<-count+bitwAnd(bitwShiftR(patterns,j-1L),1L)
    p<-sum(count<=min(pos,neg)|count>=d-min(pos,neg))/length(patterns)
  }
  if(n>=minimum_cases&&d==0L)p<-1
  if(n) {
    o<-order(x);sorted<-x[o]
    boot<-vapply(seq_len(nrow(indices)),function(b){w<-tabulate(indices[b,],nbins=n)[o];cw<-cumsum(w);(sorted[which(cw>=(n+1L)%/%2L)[1]]+sorted[which(cw>=(n+2L)%/%2L)[1]])/2},numeric(1));boot<-sort(boot)
    ci<-vapply(c(.025,.975),function(q){h<-1+(length(boot)-1)*q;l<-floor(h);u<-ceiling(h);boot[l]+(h-l)*(boot[u]-boot[l])},numeric(1))
  }else ci<-c(NA_real_,NA_real_)
  list(cases=n,positive_cases=pos,negative_cases=neg,zero_cases=n-d,median_effect=LN_audit_median_RAM(x),bootstrap_lower=ci[1],bootstrap_upper=ci[2],exact_sign_p=p)
}
LN_audit_BH_RAM <- function(p,method="BH") {
  q<-p;q[!is.finite(q)]<-1;m<-length(q);o<-order(q);rank<-seq_len(m);raw<-q[o]*m/rank
  if(method=="BY")raw<-raw*sum(1/rank)
  for(j in rev(seq_len(max(0,m-1L))))if(j>=1L)raw[j]<-min(raw[j],raw[j+1L])
  result<-numeric(m);result[o]<-pmin(1,raw);result[!is.finite(p)]<-NA_real_;result
}

audit_LN_state_spatial_family_RAM <- function(x,production,config,phenotype_config) {
  av<-function(v,k)attr(v,k,exact=TRUE);dl<-function(v){if(is.list(v)&&!isS4(v))return(v);attr(v,"listData",exact=TRUE)}
  coords<-av(x,"spatialCoords");if(is.null(coords))coords<-dl(av(x,"int_colData"))[["spatialCoords"]]
  audited_phenotypes<-list(primary=derive_LN_hierarchical_protein_states_RAM(x,phenotype_config,"counts",LN_audit_rank_RAM),
    area_adjusted=derive_LN_hierarchical_protein_states_RAM(x,phenotype_config,"counts",LN_audit_rank_RAM,area_adjust=TRUE,residual_fun=LN_audit_residual_RAM))
  phenotype_differences<-0L
  for(mode in names(audited_phenotypes))for(name in names(audited_phenotypes[[mode]]$internal$states))
    phenotype_differences<-phenotype_differences+sum(audited_phenotypes[[mode]]$internal$states[[name]]!=production$internal$phenotypes[[mode]]$internal$states[[name]])
  for(mode in names(audited_phenotypes))for(name in names(audited_phenotypes[[mode]]$internal$parents))
    phenotype_differences<-phenotype_differences+sum(audited_phenotypes[[mode]]$internal$parents[[name]]!=production$internal$phenotypes[[mode]]$internal$parents[[name]])
  metrics<-list();max_count_error<-0;max_geometry_error<-0;queries<-0L
  for(im in names(production$internal$graphs)) {
    if(Sys.time()>as.POSIXct(config$absolute_deadline_UTC,tz="UTC"))stop("FIXED_DEADLINE")
    g<-production$internal$graphs[[im]];i<-g$cell_rows;xy<-coords[i,1:2,drop=FALSE];n<-length(i)
    if(is.null(g$neighbors)) {
      stopifnot(all(vapply(audited_phenotypes,function(p)sum(p$internal$parents$CD4[i])<config$minimum_parent_cells||sum(p$internal$parents$CD8a[i])<config$minimum_parent_cells,logical(1))))
      for(mode in names(audited_phenotypes))for(k in unlist(config$k_values)) {
        key<-paste(mode,k,sep="_");if(is.null(metrics[[key]]))metrics[[key]]<-list()
        metrics[[key]][[im]]<-lapply(production$internal$pairs,function(pair)list(pair=pair,whole=NA_real_,conditional=NA_real_,counts=NULL))
      }
      next
    }
    queries_i<-unique(as.integer(round(seq(1,n,length.out=min(10L,n)))))
    for(j in queries_i) {
      delta<-sweep(xy,2,xy[j,],"-");distances<-sqrt(rowSums(delta^2));distances[j]<-Inf;truth<-sort(distances)[seq_len(ncol(g$neighbors))]
      actual<-distances[g$neighbors[j,]];stopifnot(!anyDuplicated(g$neighbors[j,]),all(g$neighbors[j,]!=j))
      max_geometry_error<-max(max_geometry_error,max(abs(actual-truth)),max(abs(actual-g$distances[j,])));queries<-queries+1L
    }
    for(k in unlist(config$k_values)) {
      neighbors<-g$neighbors[,seq_len(k),drop=FALSE]
      adj<-Matrix::sparseMatrix(i=rep(seq_len(n),times=k),j=as.vector(neighbors),x=1,dims=c(n,n))
      stopifnot(sum(Matrix::diag(adj))==0)
      for(mode in names(audited_phenotypes)) {
        p<-audited_phenotypes[[mode]]$internal;rn<-c("R0_FOXP3_CTLA4","R1_FOXP3_CTLA4_MKI67","R2_general_CD4_MKI67");en<-c("E0_PD1_TIM3","E1_CD8_MKI67","E2_CCR7_CD127")
        S<-cbind(vapply(p$states[rn],function(z)as.numeric(z[i]),numeric(n)),parent=as.numeric(p$parents$CD4[i]))
        T<-cbind(vapply(p$states[en],function(z)as.numeric(z[i]),numeric(n)),parent=as.numeric(p$parents$CD8a[i]))
        C<-as.matrix(Matrix::crossprod(S,adj%*%T));records<-list()
        for(r in rn)for(e in en) {
          eligible<-sum(S[,"parent"])>=config$minimum_parent_cells&&sum(T[,"parent"])>=config$minimum_parent_cells&&im%in%p$testable_images[[r]]&&im%in%p$testable_images[[e]]&&sum(S[,r])>=config$minimum_state_cells&&sum(T[,e])>=config$minimum_state_cells
          if(!eligible){records[[length(records)+1L]]<-list(pair=paste(r,e,sep="__"),whole=NA_real_,conditional=NA_real_,counts=NULL);next}
          a<-C[r,e];b<-C[r,"parent"]-a;c<-C["parent",e]-a;d<-C["parent","parent"]-a-b-c
          E<-sum(S[,r])*sum(T[,e])*k/(n-1);margin<-all(c(a+b,c+d,a+c,b+d)>0)
          records[[length(records)+1L]]<-list(pair=paste(r,e,sep="__"),whole=unname(log((a+.5)/(E+.5))),conditional=if(margin)unname(log((a+.5)*(d+.5)/((b+.5)*(c+.5)))) else NA_real_,counts=c(observed=unname(a),expected=E,a=unname(a),b=unname(b),c=unname(c),d=unname(d)))
        }
        key<-paste(mode,k,sep="_");if(is.null(metrics[[key]]))metrics[[key]]<-list();metrics[[key]][[im]]<-records
        prod<-production$internal$metrics[[key]][[im]]
        for(j in seq_along(records)) {
          for(field in c("whole","conditional")){v<-records[[j]][[field]];w<-prod[[j]][[field]];stopifnot(is.finite(v)==is.finite(w));if(is.finite(v))max_count_error<-max(max_count_error,abs(v-w))}
          if(!is.null(records[[j]]$counts)){stopifnot(!is.null(prod[[j]]$counts));max_count_error<-max(max_count_error,max(abs(records[[j]]$counts-prod[[j]]$counts)))}
        }
      }
    }
  }
  meta<-production$internal$meta;rows<-list();max_stat_error<-0;slot<-0L
  for(pair in production$internal$pairs) {
    rec<-lapply(metrics[[paste("primary",config$primary_k,sep="_")]][meta$image],function(z)z[[match(pair,vapply(z,function(v)v$pair,character(1)))]])
    whole<-LN_audit_cases_RAM(rec,meta,"whole");conditional<-LN_audit_cases_RAM(rec,meta,"conditional");n1<-LN_audit_cases_RAM(rec,meta,"conditional","N1");n2<-LN_audit_cases_RAM(rec,meta,"conditional","N2");paired<-intersect(names(n1),names(n2));location<-n1[paired]-n2[paired]
    for(estimand in c("whole","conditional","N1_minus_N2")) {
      slot<-slot+1L;v<-switch(estimand,whole=whole,conditional=conditional,N1_minus_N2=location);pt<-production$internal$tests[[slot]]
      stopifnot(setequal(names(v),names(pt$effects)));v<-v[names(pt$effects)];if(length(v))max_stat_error<-max(max_stat_error,max(abs(v-pt$effects)))
      if(length(v)){set.seed(config$seed+slot);expected<-matrix(sample.int(length(v),length(v)*config$bootstrap_replicates,replace=TRUE),nrow=config$bootstrap_replicates);stopifnot(identical(expected,pt$bootstrap_indices))}
      a<-LN_audit_test_RAM(v,pt$bootstrap_indices,config$minimum_cases);rows[[slot]]<-c(list(pair=pair,estimand=estimand),a)
      prod<-production$summary$primary_rows[[slot]]
      for(name in names(a)){z<-a[[name]];w<-prod[[name]];stopifnot(is.finite(z)==is.finite(w));if(is.finite(z))max_stat_error<-max(max_stat_error,abs(z-w))}
    }
  }
  p<-vapply(rows,function(z)z$exact_sign_p,numeric(1));bh<-LN_audit_BH_RAM(p);by<-LN_audit_BH_RAM(p,"BY")
  for(j in seq_along(rows)){rows[[j]]$BH<-bh[j];rows[[j]]$BY<-by[j];for(name in c("BH","BY")){a<-rows[[j]][[name]];b<-production$summary$primary_rows[[j]][[name]];stopifnot(is.finite(a)==is.finite(b));if(is.finite(a))max_stat_error<-max(max_stat_error,abs(a-b))}}
  max_sensitivity_error<-0
  for(item in production$summary$sensitivities) {
    get<-function(mode,k)lapply(metrics[[paste(mode,k,sep="_")]][meta$image],function(z)z[[match(item$pair,vapply(z,function(v)v$pair,character(1)))]])
    rr<-get(item$mode,item$k);ref<-get("primary",config$primary_k);metric<-item$estimand
    effects<-LN_audit_cases_RAM(rr,meta,metric);joint<-vapply(seq_along(rr),function(i)is.finite(rr[[i]][[metric]])&&is.finite(ref[[i]][[metric]]),logical(1))
    base<-LN_audit_cases_RAM(ref[joint],meta[joint,,drop=FALSE],metric);alt<-LN_audit_cases_RAM(rr[joint],meta[joint,,drop=FALSE],metric);common<-intersect(names(base),names(alt))
    expected<-list(cases=length(effects),median_effect=LN_audit_median_RAM(effects),common_cases=length(common),common_eligible_images=sum(joint),
      primary_median_on_common_cases=LN_audit_median_RAM(base[common]),sensitivity_median_on_common_cases=LN_audit_median_RAM(alt[common]))
    for(name in names(expected)){a<-expected[[name]];b<-item[[name]];stopifnot(is.finite(a)==is.finite(b));if(is.finite(a))max_sensitivity_error<-max(max_sensitivity_error,abs(a-b))}
  }
  list(status=if(phenotype_differences==0L&&max_count_error<1e-9&&max_stat_error<1e-9&&max_sensitivity_error<1e-9&&max_geometry_error<1e-7)"PASS_FULL_FROZEN_FAMILY_NUMERIC_AND_GRAPH_AUDIT" else "AUDIT_DIFFERENCE_REQUIRES_REVIEW",
    phenotype_membership_differences=phenotype_differences,max_sparse_count_or_image_metric_error=max_count_error,
    max_case_statistic_or_BH_error=max_stat_error,max_sensitivity_point_or_coverage_error=max_sensitivity_error,exact_distance_query_checks=queries,max_exact_distance_error=max_geometry_error,
    complete_reference_primary_rows=rows,all_graph_images_and_all_pair_counts_recalculated=TRUE,
    phenotype_definition_shared_but_numeric_rank_and_OLS_independent=TRUE,
    patient_aggregation_sign_test_bootstrap_and_BH_independent_algorithms=TRUE,
    no_individual_identifiers_case_vectors_state_arrays_or_rank_matrices_returned=TRUE)
}
