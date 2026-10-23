# Same relative state definitions; new measured CD3/CD45 parent support.
# Requires LN_rank_profiles_RAM and unit_key_RAM. No clinical fields are inputs.
derive_external_TDLN_parents_states_RAM <- function(a, panel, fine, case, roi, config,
    rank_fun=function(z)rank(z,ties.method="average")) {
  stopifnot(is.matrix(a),is.numeric(a),identical(rownames(a),as.character(panel)),
    !anyDuplicated(panel),length(fine)==ncol(a),length(case)==ncol(a),length(roi)==ncol(a))
  fine<-as.character(fine);case<-as.character(case);roi<-as.character(roi)
  valid<-function(z)!is.na(z)&!tolower(trimws(z))%in%c("","na","nan","none","null","unknown","n/a")
  take<-valid(fine)&valid(case)&valid(roi)
  image<-unit_key_RAM(list(case,roi));n<-ncol(a)
  lm<-c("CD3","CD45","CD4","CD8a","CD14","CD68","CD163","CD20","SOX10","CD31","Coll-I")
  stopifnot(all(lm%in%panel),any(take),length(unique(fine[take]))<=128L)
  # A constant CD3/CD45 channel has rank .5 everywhere, not lineage evidence.
  image_rows<-split(which(take),image[take])
  supported<-vapply(image_rows,function(i)all(vapply(c("CD3","CD45"),function(marker) {
    z<-a[marker,i];z<-z[is.finite(z)]
    length(z)>=config$minimum_cells_per_image_group&&length(unique(z))>=2L
  },logical(1))),logical(1))
  parent_take<-take&image%in%names(image_rows)[supported]
  p<-if(any(parent_take))LN_rank_profiles_RAM(a,lm,fine,image,case,parent_take,
    config$minimum_cells_per_image_group,rank_fun) else list(fine=character())
  rule<-function(r,label) {
    if(any(!is.finite(r)))return(FALSE)
    other<-if(label=="CD4")"CD8a" else "CD4"
    opponent<-max(mean(r[c("CD14","CD68","CD163")]),r[c("CD20","SOX10","CD31","Coll-I",other)])
    r["CD3"]>=config$minimum_CD3_relative_rank&&
      r["CD45"]>=config$minimum_CD45_relative_rank&&r[label]-opponent>=config$lineage_margin
  }
  parents<-list(CD4=character(),CD8a=character());records<-list()
  for(k in seq_along(p$fine)) {
    g<-which(p$group_fine==k&p$eligible);cases<-unique(p$image_patient[p$group_image[g]])
    cells<-sum(parent_take&fine==p$fine[k])
    coverage<-length(g)>=config$minimum_images&&length(cases)>=config$minimum_cases&&cells>=config$minimum_cells
    bycase<-split(g,p$image_patient[p$group_image[g]])
    agreement<-setNames(lapply(names(parents),function(label) {
      z<-vapply(bycase,function(j)mean(vapply(j,function(i)rule(p$group_profiles[,i],label),logical(1))),numeric(1))
      if(length(z))mean(z) else 0
    }),names(parents))
    passed<-character()
    for(label in names(parents))if(coverage&&rule(p$profiles[,k],label)&&agreement[[label]]>=config$minimum_case_balanced_agreement) {
      parents[[label]]<-c(parents[[label]],p$fine[k]);passed<-c(passed,label)
    }
    records[[k]]<-list(author_group_ordinal=k,candidate_lineages=as.list(passed),
      eligible_images=length(g),registered_case_labels=length(cases),cell_rows=cells,
      case_balanced_marker_agreement=agreement,relative_marker_candidate_not_formal_identity=TRUE)
  }
  stopifnot(!length(intersect(parents$CD4,parents$CD8a)))
  masks<-list(CD4=parent_take&fine%in%parents$CD4,CD8a=parent_take&fine%in%parents$CD8a)
  definitions<-list(
    R0_FOXP3_CTLA4=list(parent="CD4",markers=c("FoxP3","CTLA-4"),cut=c(.75,.60)),
    R1_FOXP3_CTLA4_MKI67=list(parent="CD4",markers=c("FoxP3","CTLA-4","Ki-67"),cut=c(.75,.60,.75)),
    R2_general_CD4_MKI67=list(parent="CD4",markers="Ki-67",cut=.75),
    E0_PD1_TIM3=list(parent="CD8a",markers=c("PD-1","TIM-3"),cut=c(.70,.70)),
    E1_CD8_MKI67=list(parent="CD8a",markers="Ki-67",cut=.75),
    E2_CCR7_CD127=list(parent="CD8a",markers=c("CCR7","CD127"),cut=c(.70,.70)))
  available<-vapply(definitions,function(d)all(d$markers%in%panel),logical(1))
  membership<-setNames(rep(list(rep(NA,n)),length(definitions)),names(definitions))
  testable<-setNames(rep(list(character()),length(definitions)),names(definitions))
  for(lineage in names(parents)) {
    inside<-masks[[lineage]];if(!any(inside))next
    names_here<-names(definitions)[vapply(definitions,function(d)d$parent==lineage,logical(1))&available]
    markers<-unique(unlist(lapply(definitions[names_here],function(d)d$markers),use.names=FALSE))
    imgs<-split(which(inside),image[inside]);ranks<-setNames(rep(list(rep(NA_real_,n)),length(markers)),markers)
    varying<-matrix(FALSE,length(imgs),length(markers),dimnames=list(names(imgs),markers))
    for(marker in markers)for(im in names(imgs)) {
      i<-imgs[[im]];z<-a[marker,i];i<-i[is.finite(z)];z<-z[is.finite(z)]
      if(length(i)<config$minimum_parent_cells_per_image||length(unique(z))<2L)next
      ranks[[marker]][i]<-rank_fun(z)/(length(i)+1);varying[im,marker]<-TRUE
    }
    for(name in names_here) {
      d<-definitions[[name]];good<-rownames(varying)[rowSums(varying[,d$markers,drop=FALSE])==length(d$markers)]
      value<-rep(FALSE,n);value[inside]<-NA
      selected<-which(inside&image%in%good)
      keep<-rep(TRUE,length(selected))
      for(j in seq_along(d$markers))keep<-keep&ranks[[d$markers[j]]][selected]>=d$cut[j]
      value[selected]<-keep;membership[[name]]<-value;testable[[name]]<-good
    }
  }
  state_records<-lapply(names(definitions),function(name) {
    d<-definitions[[name]];m<-membership[[name]];positive<-which(!is.na(m)&m)
    parent_resolved<-any(masks[[d$parent]])
    estimable<-available[[name]]&&parent_resolved&&length(testable[[name]])>0L
    status<-if(!available[[name]])"STRUCTURAL_MARKER_ABSENCE" else if(!parent_resolved)
      "PARENT_NOT_QUALIFIED" else if(!length(testable[[name]]))
      "NO_IMAGE_WITH_SUFFICIENT_VARIABLE_MARKERS" else "RELATIVE_STATE_MEASUREMENT_AVAILABLE"
    list(state=name,parent=d$parent,physical_markers=as.list(d$markers),fixed_relative_cutoffs=as.list(d$cut),
      marker_set_available=available[[name]],measurement_status=status,
      cell_rows=if(estimable)length(positive) else NA_integer_,
      registered_case_labels=if(estimable)length(unique(case[positive])) else NA_integer_,
      images_with_measurement_support=length(testable[[name]]),
      unmeasured_parent_cell_rows=if(available[[name]]&&parent_resolved)sum(is.na(m)&masks[[d$parent]]) else NA_integer_,
      missing_marker_branch_not_negative=TRUE)
  })
  list(summary=list(status="FIXED_EXTERNAL_RELATIVE_PARENT_STATE_CANDIDATES_NOT_IDENTITY_OR_FUNCTION",
    author_group_reference="fixed provided clustering; ordinal records only, not original labels",
    parent_records=records,resolved_parent_group_counts=lapply(parents,length),state_records=state_records,
    missing_fine_or_unit_rows=sum(!take),parent_selection_excludes_target_state_markers=TRUE,
    images_without_variable_CD3_or_CD45_support=sum(!supported),
    parent_rows_without_variable_CD3_or_CD45_support=sum(take&!parent_take),
    parent_agreement_is_case_balanced_across_images=TRUE,
    state_cutoffs_unchanged_from_previous_relative_construct=TRUE,
    missing_or_unmeasured_states_retained_as_NA=TRUE,
    internal_individual_masks_rank_profiles_and_identifiers_must_not_be_serialized=TRUE),
    internal=list(states=membership,parents=masks,testable_images=testable,
      measurement_available=available,fine=fine,image=image,case=case,eligible=take))
}
