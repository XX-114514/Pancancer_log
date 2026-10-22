# Current numeric qualifier: small finite sets do not release quantiles.
# Legacy v1 was synthetic only and was never adopted for a real object.
# A later reviewed native schema contract must provide the assay and row map.
# No native object lookup, external backing-file access or output files here.

LN_predeclared_marker_names <- c("CD3","CD4","CD8","FOXP3","PD1","Ki67",
  "GZMB","GZMK","TCF1","CCR7","CD25","LAG3","TIM3","CTLA4","CCR8",
  "CXCL13","LAMP3","CD1c","CD141","CD68","CD163","CD20","HLA_DR",
  "CD11c","CD45","CD45RA","CD45RO","TIA1")

profile_LN_marker_assay_RAM <- function(assay, marker_rows, ROI_key=NULL) {
  stopifnot((is.matrix(assay)&&is.numeric(assay))||inherits(assay,"dgCMatrix"))
  stopifnot(nrow(assay)>0L,nrow(assay)<=256L,ncol(assay)>0L,
    is.list(marker_rows),identical(names(marker_rows),LN_predeclared_marker_names))
  for (v in marker_rows) stopifnot(is.numeric(v),!anyNA(v),
    all(is.finite(v)),all(v==as.integer(v)),all(v>=1L & v<=nrow(assay)),!anyDuplicated(v))
  assigned <- unlist(marker_rows,use.names=FALSE)
  stopifnot(!anyDuplicated(assigned))
  units <- if(is.null(ROI_key)) list() else {
    stopifnot(is.character(ROI_key),length(ROI_key)==ncol(assay))
    ok <- !is.na(ROI_key) & !tolower(trimws(ROI_key)) %in%
      c("","na","nan","none","null","unknown","n/a")
    split(which(ok),ROI_key[ok])
  }
  records <- lapply(names(marker_rows),function(name) {
    rows <- marker_rows[[name]]
    if(length(rows)!=1L) return(list(marker=name,matched_feature_rows=length(rows),
      status=if(length(rows))"AMBIGUOUS_FEATURE_ROWS_NO_AUTOMATIC_COMBINATION" else "MISSING_FEATURE_NOT_ZERO_EXPRESSION",
      numeric_profile=NULL))
    values <- as.numeric(assay[rows,,drop=TRUE])
    stopifnot(length(values)==ncol(assay))
    finite <- is.finite(values)
    x <- values[finite]
    within <- lapply(units,function(i) {
      y <- values[i];y<-y[is.finite(y)]
      c(two_finite=length(y)>=2L,variable=length(y)>=2L && min(y)<max(y))
    })
    summary <- list(total_cell_rows=length(values),finite_rows=length(x),
      nonfinite_rows=sum(!finite),negative_numeric_rows=sum(x<0),
      zero_numeric_rows=sum(x==0),strictly_positive_numeric_rows=sum(x>0),
      fractional_numeric_rows=sum(x!=trunc(x)),finite_unique_values=length(unique(x)),
      nonconstant_finite_values=length(x)>1L && min(x)<max(x),
      quantile_release_minimum_finite_rows=20L,
      quantiles_withheld_for_small_finite_set=length(x)<20L,
      quantiles=if(length(x)>=20L)as.list(setNames(as.numeric(quantile(x,
        c(.05,.25,.5,.75,.95),names=FALSE,type=7)),c("q05","q25","median","q75","q95"))) else NULL,
      explicitly_supplied_ROI_keys=length(units),
      ROIs_with_two_finite_values=if(length(within))sum(vapply(within,function(v)v[1],logical(1))) else 0L,
      ROIs_with_nonconstant_values=if(length(within))sum(vapply(within,function(v)v[2],logical(1))) else 0L,
      positive_number_not_accepted_marker_positivity=TRUE)
    list(marker=name,matched_feature_rows=1L,
      status=if(length(x))"MEASURED_NUMERIC_PROFILE_ONLY" else "MEASURED_NO_FINITE_VALUES",numeric_profile=summary)
  })
  list(status="COMPLETE_PREDECLARED_MARKER_NUMERIC_QUALIFICATION_NO_STATE_OR_EFFECT",
    recognized_assay_storage=if(is.matrix(assay))"IN_MEMORY_DENSE_MATRIX" else "IN_MEMORY_DGCMATRIX",
    features=nrow(assay),cells=ncol(assay),predeclared_markers=records,
    field_mapping_semantics_require_separate_native_review=TRUE,
    assay_transform_background_and_physical_units_not_inferred=TRUE,
    no_threshold_or_cell_state_or_clinical_effect_accepted=TRUE,
    duplicate_features_not_summed_and_missing_markers_not_zero_filled=TRUE,
    no_identifiers_ROI_records_or_per_cell_values_in_result=TRUE)
}

check_LN_marker_assay_synthetic_RAM <- function() {
  x <- rbind(c(-1,0,.5,1,NA,Inf),rep(0,6),rep(NA_real_,6),
    seq_len(6),seq_len(6)+1,rep(2,6))
  rows <- setNames(rep(list(integer()),length(LN_predeclared_marker_names)),LN_predeclared_marker_names)
  rows[c("CD3","CD4","CD8","PD1","Ki67")] <- list(1L,2L,3L,c(4L,5L),6L)
  u <- c("FAKE_ROI_A","FAKE_ROI_A","FAKE_ROI_A","FAKE_ROI_B","FAKE_ROI_B",NA_character_)
  a <- profile_LN_marker_assay_RAM(x,rows,u)
  get <- function(z,k)z$predeclared_markers[[match(k,LN_predeclared_marker_names)]]
  v <- get(a,"CD3")$numeric_profile
  stopifnot(v$finite_rows==4L,v$nonfinite_rows==2L,v$negative_numeric_rows==1L,
    v$zero_numeric_rows==1L,v$strictly_positive_numeric_rows==2L,v$fractional_numeric_rows==1L,
    v$ROIs_with_two_finite_values==1L,v$ROIs_with_nonconstant_values==1L)
  # Scalar sign/fraction counts and independently interpolated order statistics.
  counts <- c(finite=0,negative=0,zero=0,positive=0,fractional=0)
  finite <- numeric()
  for(y in x[1,]) if(is.finite(y)) {
    counts <- counts+c(1,y<0,y==0,y>0,y!=trunc(y));finite<-c(finite,y)
  }
  stopifnot(identical(as.numeric(counts),c(4,1,1,2,1)))
  stopifnot(is.null(v$quantiles),v$quantiles_withheld_for_small_finite_set)
  many <- matrix(seq(-3,3,length.out=31),nrow=1)
  many_rows <- setNames(rep(list(integer()),length(LN_predeclared_marker_names)),LN_predeclared_marker_names)
  many_rows$CD3 <- 1L
  qm <- get(profile_LN_marker_assay_RAM(many,many_rows),"CD3")$numeric_profile
  ordered <- sort(as.numeric(many))
  oracle <- vapply(c(.05,.25,.5,.75,.95),function(p) {
    h<-1+(length(ordered)-1)*p;lo<-floor(h);hi<-ceiling(h)
    ordered[lo]+(h-lo)*(ordered[hi]-ordered[lo])
  },numeric(1))
  stopifnot(max(abs(unlist(qm$quantiles,use.names=FALSE)-oracle))<1e-12)
  q19 <- get(profile_LN_marker_assay_RAM(many[,1:19,drop=FALSE],many_rows),"CD3")$numeric_profile
  q20 <- get(profile_LN_marker_assay_RAM(many[,1:20,drop=FALSE],many_rows),"CD3")$numeric_profile
  stopifnot(is.null(q19$quantiles),!is.null(q20$quantiles),
    q19$quantiles_withheld_for_small_finite_set,!q20$quantiles_withheld_for_small_finite_set)
  stopifnot(get(a,"CD4")$numeric_profile$zero_numeric_rows==6L,
    !get(a,"CD4")$numeric_profile$nonconstant_finite_values,
    get(a,"CD8")$status=="MEASURED_NO_FINITE_VALUES",
    is.null(get(a,"CD8")$numeric_profile$quantiles),
    get(a,"FOXP3")$status=="MISSING_FEATURE_NOT_ZERO_EXPRESSION",
    is.null(get(a,"FOXP3")$numeric_profile),
    get(a,"PD1")$status=="AMBIGUOUS_FEATURE_ROWS_NO_AUTOMATIC_COMBINATION",
    is.null(get(a,"PD1")$numeric_profile),
    get(a,"Ki67")$numeric_profile$strictly_positive_numeric_rows==6L,
    !get(a,"Ki67")$numeric_profile$nonconstant_finite_values)
  sparse <- as(Matrix::Matrix(x,sparse=TRUE),"dgCMatrix")
  b <- profile_LN_marker_assay_RAM(sparse,rows,u)
  b$recognized_assay_storage <- a$recognized_assay_storage
  stopifnot(identical(a,b))
  perm <- c(6,3,1,5,4,2)
  stopifnot(identical(a,profile_LN_marker_assay_RAM(x[,perm,drop=FALSE],rows,u[perm])))
  bad <- rows;bad$FOXP3<-1L
  stopifnot(inherits(try(profile_LN_marker_assay_RAM(x,bad,u),silent=TRUE),"try-error"))
  stopifnot(!grepl("FAKE_ROI",jsonlite::toJSON(a,auto_unbox=TRUE,null="null")))
  list(status="PASS_SYNTHETIC_MARKER_MISSING_ZERO_NONFINITE_AMBIGUOUS_AND_TRANSFORMED_VALUES",
    scalar_counts_and_order_statistics_independently_checked=TRUE,
    dense_sparse_and_column_permutation_equivalent=TRUE,
    small_finite_sets_withheld_and_release_boundary_checked=TRUE,
    marker_row_collisions_rejected=TRUE,no_real_native_object_read=TRUE,
    no_state_thresholds_or_effects_accepted=TRUE)
}
