# Full secondary table/null replay; sampled exact geometry, not biological identity validation.
LN_spec_audit_cases <- function(values,meta) {
  patients<-unique(meta$patient);w<-matrix(0,nrow(meta),length(patients),dimnames=list(NULL,patients))
  for(p in patients){nodes<-unique(meta$node[meta$patient==p]);for(nd in nodes){i<-which(meta$patient==p&meta$node==nd);w[i,p]<-1/(length(nodes)*length(i))}}
  if(!length(patients))return(matrix(numeric(),0,ncol(values)))
  crossprod(w,values)
}
LN_spec_audit_effects <- function(tab) {
  z<-tab+.5;den<-apply(z,1,sum)
  a<-log(z[1,1])+log(z[1,4])-log(z[1,2])-log(z[1,3]);b<-log(z[2,1])+log(z[2,4])-log(z[2,2])-log(z[2,3])
  c(a,b,log(sum(z[,1]*z[,4]/den))-log(sum(z[,2]*z[,3]/den)),b-a)
}
LN_spec_audit_null <- function(cases,indices,minimum_cases,infer) {
  n<-nrow(cases);B<-ncol(cases)-1L;center<-if(n)apply(cases,1,sum)/(B+1) else numeric()
  corrected<-cases[,1]-center;nt<-if(n)vapply(seq_len(B),function(b)LN_audit_median_RAM(cases[,b+1]-center),numeric(1)) else rep(NA_real_,B)
  pt<-LN_audit_median_RAM(corrected);low<-if(n)sum(nt<=pt+1e-12) else 0;high<-if(n)sum(nt>=pt-1e-12) else 0;variable<-n&&max(nt)-min(nt)>1e-12
  pv<-if(infer&&n>=minimum_cases&&variable)min(1,2*(min(low,high)+1)/(B+1)) else NA_real_
  ci<-if(n&&infer)LN_audit_test_RAM(corrected,indices,minimum_cases) else list(bootstrap_lower=NA_real_,bootstrap_upper=NA_real_)
  list(cases=n,median_effect=pt,bootstrap_lower=ci$bootstrap_lower,bootstrap_upper=ci$bootstrap_upper,
    positive_cases=sum(corrected>1e-12),negative_cases=sum(corrected< -1e-12),zero_cases=sum(abs(corrected)<=1e-12),permutation_p=pv,lower_tail_count=low,upper_tail_count=high)
}
audit_LN_regulatory_specificity_RAM <- function(x,production,config) {
  av<-function(v,k)attr(v,k,exact=TRUE);dl<-function(v){if(is.list(v)&&!isS4(v))return(v);attr(v,"listData",exact=TRUE)}
  xy<-av(x,"spatialCoords");if(is.null(xy))xy<-dl(av(x,"int_colData"))[["spatialCoords"]]
  p<-production$internal$phenotype$internal;images<-production$internal$images;direct<-production$internal$direct;nulls<-production$internal$nulls
  ends<-dimnames(direct)[[2]];maxerr<-0;distanceerr<-0;queries<-0L;replayed<-0L;block_checks<-0L;margin_flags<-integer(3)
  compare<-function(a,b){stopifnot(identical(dim(a),dim(b)),identical(as.vector(is.finite(a)),as.vector(is.finite(b))));if(any(is.finite(a)))max(abs(a[is.finite(a)]-b[is.finite(b)])) else 0}
  for(im in names(images)) {
    if(Sys.time()>as.POSIXct(config$absolute_deadline_UTC,tz="UTC"))stop("FIXED_DEADLINE")
    z<-images[[im]];if(is.null(z$neighbors))next
    i<-z$cell_rows;coord<-xy[i,1:2,drop=FALSE];n<-length(i);nb<-z$neighbors;anchor<-z$anchor;reg<-z$reg;cycle<-z$cycle
    for(j in unique(round(seq(1,n,length.out=min(n,10L))))) {
      d<-sqrt((coord[,1]-coord[j,1])^2+(coord[,2]-coord[j,2])^2);d[j]<-Inf
      stopifnot(!anyDuplicated(nb[j,]),all(nb[j,]!=j));truth<-sort(d)[seq_len(config$k)]
      distanceerr<-max(distanceerr,max(abs(d[nb[j,]]-truth)),max(abs(z$distances[j,]-truth)));queries<-queries+1L
    }
    masks<-cbind(z$Q,vapply(ends,function(e)p$states[[e]][i],logical(n)))
    D<-matrix(0,length(anchor),4)
    for(a in seq_along(anchor))for(e in seq_len(4))D[a,e]<-sum(masks[nb[anchor[a],],e])
    maxerr<-max(maxerr,compare(D,unname(z$D)))
    totals<-rbind(colSums(D[!cycle,,drop=FALSE]),colSums(D[cycle,,drop=FALSE]))
    make_tables<-function(rr,e){tab<-matrix(0,2,4);for(cy in 0:1){k<-which(cycle==as.logical(cy));a<-sum(D[k,e]*rr[k]);b<-sum((D[k,1]-D[k,e])*rr[k]);c<-sum(D[k,e]*(!rr[k]));d<-sum((D[k,1]-D[k,e])*(!rr[k]));tab[cy+1,]<-c(a,b,c,d)};tab}
    for(e in which(z$eligible)){
      tab<-make_tables(reg,e+1L);maxerr<-max(maxerr,compare(tab,unname(z$tables[[e]])),compare(unname(LN_spec_audit_effects(tab)),unname(direct[im,e,])))
      bad<-FALSE;for(row in seq_len(2)){v<-tab[row,];if(any(c(v[1]+v[2],v[3]+v[4],v[1]+v[3],v[2]+v[4])==0))bad<-TRUE};margin_flags[e]<-margin_flags[e]+bad
    }
    for(grid in config$grids) {
      gkey<-as.character(grid);blocks<-z$blocks[[gkey]]
      box<-apply(coord,2,range);sx<-floor((coord[anchor,1]-box[1,1])/(box[2,1]-box[1,1])*grid);sy<-floor((coord[anchor,2]-box[1,2])/(box[2,2]-box[1,2])*grid)
      deg<-D[,1];db<-1L+(deg>0)+(deg>1)+(deg>3)+(deg>6)
      truth<-split(seq_along(anchor),paste(pmin(sx,grid-1L),pmin(sy,grid-1L),as.integer(cycle),db,sep=":"));stopifnot(identical(truth,blocks))
      if(!any(is.finite(nulls[[gkey]][im,,])))next
      set.seed(z$seeds[[gkey]])
      for(b in seq_len(config$permutations)) {
        rr<-reg;for(j in blocks)if(length(j)>1L)rr[j]<-sample(reg[j],length(j),replace=FALSE)
        stopifnot(sum(rr&cycle)==sum(reg&cycle),sum(rr&!cycle)==sum(reg&!cycle))
        if(b%in%c(1L,(config$permutations+1L)%/%2L,config$permutations)){for(j in blocks)stopifnot(sum(rr[j])==sum(reg[j]));block_checks<-block_checks+length(blocks)}
        for(e in which(z$eligible)){v<-LN_spec_audit_effects(make_tables(rr,e+1L))[3];maxerr<-max(maxerr,abs(v-nulls[[gkey]][im,b,e]));replayed<-replayed+1L}
      }
    }
  }
  for(e in seq_len(3))stopifnot(margin_flags[e]==production$summary$edge_margin_quality[[e]]$images_with_any_zero_edge_margin_stratum)
  meta<-production$internal$meta;reference<-list();slot<-0L
  for(e in seq_len(3))for(est in seq_len(4)) {
    slot<-slot+1L;ok<-is.finite(direct[,e,est]);cases<-LN_spec_audit_cases(matrix(direct[ok,e,est],ncol=1),meta[ok,,drop=FALSE])[,1];pt<-production$internal$tests[[slot]]
    stopifnot(setequal(names(cases),names(pt$effects)));cases<-cases[names(pt$effects)]
    maxerr<-max(maxerr,compare(unname(cases),unname(pt$effects)))
    a<-LN_audit_test_RAM(cases,pt$bootstrap_indices,config$minimum_cases);ref<-production$summary$primary_rows[[slot]]
    for(k in names(a))maxerr<-max(maxerr,compare(a[[k]],ref[[k]]))
    reference[[slot]]<-c(a,list(primary_p=a$exact_sign_p))
  }
  for(grid in config$grids)for(e in seq_len(3)) {
    vals<-matrix(nulls[[as.character(grid)]][,,e],nrow=nrow(meta),ncol=config$permutations);ok<-is.finite(direct[,e,"combined_MH"])&rowSums(is.finite(vals))==config$permutations
    cases<-LN_spec_audit_cases(cbind(direct[ok,e,"combined_MH"],vals[ok,,drop=FALSE]),meta[ok,,drop=FALSE]);key<-paste(grid,e,sep="_")
    target<-production$internal$null_tests[[key]];stopifnot(setequal(rownames(cases),rownames(target$case_matrix)));cases<-cases[rownames(target$case_matrix),,drop=FALSE]
    maxerr<-max(maxerr,compare(unname(cases),unname(target$case_matrix)));infer<-grid==config$primary_grid
    a<-LN_spec_audit_null(cases,target$summary$internal$bootstrap_indices,config$minimum_cases,infer);ref<-target$summary$summary
    for(k in names(a))maxerr<-max(maxerr,compare(a[[k]],ref[[k]]))
    if(infer){slot<-slot+1L;reference[[slot]]<-c(a,list(primary_p=a$permutation_p))}
  }
  pv<-vapply(reference,function(x)x$primary_p,numeric(1));bh<-LN_audit_BH_RAM(pv);by<-LN_audit_BH_RAM(pv,"BY")
  for(j in seq_along(reference)){maxerr<-max(maxerr,compare(bh[j],production$summary$primary_rows[[j]]$BH),compare(by[j],production$summary$primary_rows[[j]]$BY))}
  list(status=if(maxerr<1e-8&&distanceerr<1e-7)"PASS_FULL_SECONDARY_SPECIFICITY_NUMERIC_AUDIT" else "AUDIT_DIFFERENCE_REQUIRES_REVIEW",
    maximum_numeric_error=maxerr,maximum_exact_distance_error=distanceerr,exact_distance_queries=queries,
    all_eligible_null_image_endpoint_statistics_replayed=replayed,selected_permutation_block_count_checks=block_checks,
    direct_tables_all_recalculated=TRUE,all_patient_aggregation_bootstrap_tail_probabilities_and_BH_recalculated=TRUE,
    no_individual_or_graph_case_resample_arrays_returned=TRUE)
}
