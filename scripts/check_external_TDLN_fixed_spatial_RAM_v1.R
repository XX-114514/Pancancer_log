# Full-size resampling constants on artificial marker/geometry fixtures only.
source("scripts/qualify_LN_unit_geometry_RAM_v1.R")
source("scripts/derive_LN_hierarchical_protein_states_RAM_v2.R")
source("scripts/derive_external_TDLN_parents_states_RAM_v1.R")
source("scripts/measure_LN_state_spatial_family_RAM_v1.R")
source("scripts/measure_LN_regulatory_specificity_RAM_v1.R")
source("scripts/audit_LN_state_spatial_family_RAM_v1.R")
source("scripts/audit_LN_regulatory_specificity_RAM_v1.R")
source("scripts/measure_external_TDLN_fixed_spatial_RAM_v1.R")
source("scripts/audit_external_TDLN_fixed_spatial_RAM_v1.R")

external_TDLN_synthetic_fixture <- function() {
  one<-rbind(data.frame(fine="D4",j=1:320,replicate=0),
    data.frame(fine="D8",j=rep(1:320,each=10),replicate=rep(1:10,320)),
    do.call(rbind,lapply(c("B","MY","TUM","END","FIB"),function(f)
      data.frame(fine=f,j=1:320,replicate=0))))
  meta<-do.call(rbind,lapply(1:7,function(ca)do.call(rbind,lapply(1:2,function(im) {
    z<-one;z$case<-paste0("SYNTHETIC_CASE_CANARY_",ca)
    z$roi<-paste0("SYNTHETIC_ROI_CANARY_",im);z
  }))))
  panel<-c("CD3","CD45","CD4","CD8a","CD14","CD68","CD163","CD20","SOX10","CD31","Coll-I",
    "FoxP3","CTLA-4","Ki-67","PD-1","TIM-3","CD127","TCF1")
  a<-matrix(1,length(panel),nrow(meta),dimnames=list(panel,NULL))
  T<-meta$fine%in%c("D4","D8");a["CD3",T]<-10
  a["CD45",meta$fine%in%c("D4","D8","B","MY")]<-10
  for(item in list(c("CD4","D4"),c("CD8a","D8"),c("CD20","B"),
    c("CD14","MY"),c("CD68","MY"),c("CD163","MY"),c("SOX10","TUM"),c("CD31","END"),c("Coll-I","FIB")))
    a[item[1],meta$fine==item[2]]<-10
  score<-(meta$j*73)%%320;cycle<-(meta$j*13+23)%%320
  d4<-meta$fine=="D4";d8<-meta$fine=="D8"
  for(marker in c("FoxP3","CTLA-4"))a[marker,d4]<-1+score[d4]/5
  a["Ki-67",d4]<-1+cycle[d4]/5
  for(marker in c("PD-1","TIM-3"))a[marker,d8]<-1+score[d8]/5
  a["Ki-67",d8]<-1+(319-score[d8])/5
  a["CD127",d8]<-1+score[d8]/5;a["TCF1",d8]<-1+score[d8]/5
  yy<-c(D4=0,D8=.2,B=50,MY=70,TUM=90,END=110,FIB=130)
  xy<-cbind(meta$j*3+ifelse(d8,.05+meta$replicate*.001,0),yy[meta$fine]+meta$replicate*.001)
  list(a=a,meta=meta,xy=xy,panel=panel,score=score,cycle=cycle)
}

check_external_TDLN_fixed_spatial_RAM <- function() {
  f<-external_TDLN_synthetic_fixture();m<-f$meta
  pc<-list(minimum_cells_per_image_group=20L,minimum_images=5L,minimum_cases=3L,
    minimum_cells=100L,minimum_CD3_relative_rank=.5,minimum_CD45_relative_rank=.5,
    lineage_margin=.1,minimum_case_balanced_agreement=.6,minimum_parent_cells_per_image=20L)
  derive<-function(a=f$a,config=pc,rank_fun=function(z)rank(z,ties.method="average"))
    derive_external_TDLN_parents_states_RAM(a,f$panel,m$fine,m$case,m$roi,config,rank_fun)
  ph<-derive();oracle<-derive(rank_fun=LN_audit_rank_RAM)
  stopifnot(identical(ph$internal$parents,oracle$internal$parents),identical(ph$internal$states,oracle$internal$states),
    identical(ph$internal$parents$CD4,m$fine=="D4"),identical(ph$internal$parents$CD8a,m$fine=="D8"),
    identical(ph$internal$states$R0_FOXP3_CTLA4,m$fine=="D4"&f$score>=240),
    identical(ph$internal$states$E0_PD1_TIM3,m$fine=="D8"&f$score>=224),
    identical(ph$internal$states$E1_CD8_MKI67,m$fine=="D8"&f$score<=79),
    all(is.na(ph$internal$states$E2_CCR7_CD127)))
  target<-c("FoxP3","CTLA-4","Ki-67","PD-1","TIM-3","CD127","TCF1")
  altered<-f$a;for(marker in target)altered[marker,]<-rev(altered[marker,])
  stopifnot(identical(ph$internal$parents,derive(altered)$internal$parents))
  noCD3<-f$a;noCD3["CD3",m$fine%in%c("D4","D8")]<-0
  absent<-derive(noCD3)
  stopifnot(!any(absent$internal$parents$CD4),!any(absent$internal$parents$CD8a),
    all(vapply(absent$summary$state_records,function(z)is.na(z$cell_rows),logical(1))))
  constant<-f$a;constant["Ki-67",m$fine=="D8"]<-3
  cp<-derive(constant);cr<-Filter(function(z)z$state=="E1_CD8_MKI67",cp$summary$state_records)[[1]]
  stopifnot(cr$measurement_status=="NO_IMAGE_WITH_SUFFICIENT_VARIABLE_MARKERS",is.na(cr$cell_rows))
  cfg<-list(k=10L,minimum_parent=20L,minimum_class=10L,minimum_mixed_blocks=2L,
    minimum_movable_class=10L,grids=c(2L,4L,8L),primary_grid=4L,permutations=999L,
    bootstrap_replicates=2000L,minimum_cases=5L,seed=20261023L,
    absolute_deadline_UTC=format(Sys.time()+1200,tz="UTC",usetz=FALSE))
  prod<-measure_external_TDLN_fixed_spatial_RAM(f$xy,m$case,m$roi,ph,cfg)
  audit<-audit_external_TDLN_fixed_spatial_RAM(f$xy,prod,cfg)
  stopifnot(audit$status=="PASS_FULL_EXTERNAL_FIXED_SPATIAL_NUMERIC_AUDIT",length(prod$summary$primary_rows)==15L)
  rows<-prod$summary$primary_rows
  stopifnot(rows[[3]]$median_effect>0,rows[[7]]$median_effect<0,
    rows[[13]]$median_effect>0,rows[[14]]$median_effect<0,
    all(vapply(c(9L,10L,11L,12L,15L),function(j)rows[[j]]$status=="STRUCTURAL_MARKER_ABSENCE"&&
      !is.finite(rows[[j]]$primary_p)&&!is.finite(rows[[j]]$median_effect),logical(1))))
  # Explicit unequal-ROI case weighting oracle, independent of native labels.
  values<-matrix(c(0,10,20,100),ncol=1);md<-data.frame(case=c("A","A","A","B"))
  stopifnot(identical(unname(ETDLN_case_ROI_means(values,md)[,1]),c(10,100)),
    max(abs(ETDLN_case_ROI_means(values,md)-ETDLN_audit_case_ROI_means(values,md)))<1e-12)
  ix<-matrix(rep(1:32,2000),nrow=2000,byrow=TRUE)
  upper<-ETDLN_audit_direction_RAM(rep(1,32),ix)
  stopifnot(upper$exact_sign_p==2^(1-32),upper$bootstrap_lower==1,upper$bootstrap_upper==1)
  zero<-LN_exact_direction_summary_RAM(rep(0,7),7,2000L)
  za<-ETDLN_audit_direction_RAM(rep(0,7),zero$internal$bootstrap_indices)
  stopifnot(zero$summary$exact_sign_p==1,za$exact_sign_p==1)
  deg<-LN_spec_null_summary(rep(0,7),matrix(0,7,999),7)
  stopifnot(deg$summary$status=="DEGENERATE_CONDITIONAL_NULL",!is.finite(deg$summary$permutation_p))
  first<-which(m$case==m$case[1]);small<-ph
  small$internal$states<-lapply(ph$internal$states,function(z)z[first])
  small$internal$parents<-lapply(ph$internal$parents,function(z)z[first]);small$internal$image<-ph$internal$image[first]
  one<-measure_external_TDLN_fixed_spatial_RAM(f$xy[first,],m$case[first],m$roi[first],small,cfg)
  oa<-audit_external_TDLN_fixed_spatial_RAM(f$xy[first,],one,cfg)
  stopifnot(oa$status=="PASS_FULL_EXTERNAL_FIXED_SPATIAL_NUMERIC_AUDIT",
    all(vapply(one$summary$primary_rows,function(z)!is.finite(z$primary_p),logical(1))))
  flat<-f$xy[first,,drop=FALSE];flat[,1]<-0
  noGeometry<-measure_external_TDLN_fixed_spatial_RAM(flat,m$case[first],m$roi[first],small,cfg)
  stopifnot(all(vapply(noGeometry$summary$primary_rows,function(z)!is.finite(z$primary_p),logical(1))),
    noGeometry$summary$coverage_reasons$DEGENERATE_IMAGE_GEOMETRY==4L)
  output<-jsonlite::toJSON(list(phenotype=ph$summary,production=prod$summary,audit=audit),auto_unbox=TRUE,null="null")
  stopifnot(!grepl("SYNTHETIC_CASE_CANARY|SYNTHETIC_ROI_CANARY",output))
  list(status="PASS_FULL_SYNTHETIC_EXTERNAL_PARENT_STATE_FIXED_SPATIAL_AND_AUDIT",
    runtime_R=as.character(getRversion()),full_permutations=999L,bootstrap_replicates=2000L,
    planned_primary_slots=15L,structural_missing_branch_NA_not_negative=TRUE,
    target_marker_perturbation_does_not_change_parents=TRUE,
    independent_rank_and_state_masks_match=TRUE,positive_negative_zero_and_unestimable_controls=TRUE,
    exact_tail_recurrence_checked_without_exponential_enumeration=TRUE,
    full_numeric_audit=audit,no_real_native_or_clinical_or_individual_data_read=TRUE)
}
cat(jsonlite::toJSON(check_external_TDLN_fixed_spatial_RAM(),auto_unbox=TRUE,null="null"),"\n")
