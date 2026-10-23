# Explicit vectors only; no clinical fields, native accessors or individual output.
# Requires the unchanged unit_key_RAM helper from qualify_LN_unit_geometry_RAM_v1.R.
qualify_external_TDLN_units_RAM <- function(case, roi, image, sample, object, xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), ncol(xy)==2L, nrow(xy)>0L,
    is.numeric(object), length(object)==nrow(xy))
  fields<-list(case=case,roi=roi,image=image,sample=sample)
  stopifnot(all(lengths(fields)==nrow(xy)),
    all(vapply(fields,function(z)is.character(z)||is.factor(z),logical(1))))
  fields<-lapply(fields,as.character)
  valid<-function(z)!is.na(z)&!tolower(trimws(z))%in%
    c("","na","nan","none","null","unknown","n/a")
  good<-lapply(fields,valid);core<-good$case&good$roi
  key<-unit_key_RAM(fields[c("case","roi")])
  count_distinct<-function(z,ok)length(unique(z[ok]))
  nesting<-function(child,parent,ok) {
    edges<-unique(data.frame(child=child[ok],parent=parent[ok]))
    n<-table(edges$child)
    list(child_labels=length(n),child_labels_with_multiple_parents=sum(n>1L),
      maximum_parents_per_child=if(length(n))as.integer(max(n)) else NULL)
  }
  aliases<-lapply(c("image","sample"),function(name) {
    alias<-unit_key_RAM(list(fields$case,fields[[name]]))
    ok<-core&good[[name]]
    a<-nesting(key,alias,ok);b<-nesting(alias,key,ok)
    complete<-sum(core)>0L&&all(good[[name]][core])
    partition<-sum(ok)>0L&&a$child_labels_with_multiple_parents==0L&&
      b$child_labels_with_multiple_parents==0L
    list(alias_field_role=name,jointly_complete_rows=sum(ok),
      missing_alias_rows_among_primary_complete=sum(core&!good[[name]]),
      primary_to_alias=a,alias_to_primary=b,
      same_partition_on_jointly_complete_rows=partition,
      full_primary_coverage_and_same_partition=complete&&partition)
  })
  object_valid<-is.finite(object)&object>=0&object==floor(object)
  ok<-core&object_valid
  duplicated_objects<-sum(duplicated(data.frame(key=key[ok],object=object[ok])))
  finite_xy<-is.finite(xy[,1])&is.finite(xy[,2])
  eligible<-which(core&finite_xy)
  groups<-split(eligible,key[eligible])
  sizes<-lengths(groups)
  duplicate_xy<-sum(duplicated(data.frame(key=key[eligible],x=xy[eligible,1],y=xy[eligible,2])))
  degenerate<-vapply(groups,function(i)length(i)<2L||
    diff(range(xy[i,1]))<=0||diff(range(xy[i,2]))<=0,logical(1))
  reasons<-character()
  if(!all(core))reasons<-c(reasons,"MISSING_PRIMARY_CASE_OR_ROI")
  if(!all(vapply(aliases,function(z)z$full_primary_coverage_and_same_partition,logical(1))))
    reasons<-c(reasons,"ALIAS_PARTITION_OR_COVERAGE_UNQUALIFIED")
  if(!all(object_valid))reasons<-c(reasons,"INVALID_NUMERIC_OBJECT_KEY")
  if(duplicated_objects>0L)reasons<-c(reasons,"DUPLICATE_CASE_ROI_OBJECT_KEY")
  if(!all(finite_xy))reasons<-c(reasons,"NONFINITE_COORDINATES")
  list(status=if(!length(reasons))"REGISTERED_CASE_ROI_PARTITION_CANDIDATE" else
      "UNIT_PARTITION_REQUIRES_REVIEW_NO_AUTOMATIC_REPAIR",
    failure_reasons=as.list(reasons),input_cell_rows=nrow(xy),
    missing_rows_by_field=lapply(good,function(z)sum(!z)),
    primary_complete_rows=sum(core),
    literal_case_labels=count_distinct(fields$case,good$case),
    literal_ROI_labels=count_distinct(fields$roi,good$roi),
    primary_complete_case_ROI_tuples=count_distinct(key,core),
    raw_ROI_to_case=nesting(fields$roi,fields$case,core),alias_checks=aliases,
    invalid_object_rows=sum(!object_valid),duplicate_case_ROI_object_rows=duplicated_objects,
    finite_coordinate_rows=sum(finite_xy),primary_complete_finite_coordinate_tuples=length(groups),
    repeated_coordinate_rows_within_primary_ROI=duplicate_xy,
    primary_ROIs_with_single_cell_or_zero_axis_span=sum(degenerate),
    primary_ROIs_with_nonzero_two_axis_span=sum(!degenerate),
    primary_ROI_cell_count_quantiles=if(length(sizes))as.list(setNames(
      as.numeric(quantile(sizes,c(0,.25,.5,.75,1),names=FALSE,type=7)),
      c("min","q25","median","q75","max"))) else NULL,
    registered_case_ROI_partition_gate_passed=!length(reasons),
    coordinate_ties_reported_not_deleted=TRUE,
    alias_comparison_uses_patient_qualified_keys_not_raw_labels=TRUE,
    literal_case_labels_not_accepted_independent_patients=TRUE,
    pathology_coordinate_units_lineage_and_endpoint_not_accepted=TRUE,
    no_individual_identifiers_keys_rows_or_vectors_in_result=TRUE,
    no_native_clinical_or_state_fields_selected=TRUE)
}
