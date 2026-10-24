# A single verified native object, one finite workflow, source aggregates only.
# Generic configuration comes from the environment; native values stay RAM.
external_TDLN_native_workflow_RAM <- function(x,config) {
  av<-function(v,k)attr(v,k,exact=TRUE)
  dl<-function(v){if(is.list(v)&&!isS4(v))v else av(v,"listData")}
  aa<-dl(av(av(x,"assays"),"data"));cd<-dl(av(x,"colData"))
  required<-c("Patient.ID","ROI.ID","Image.name","sample_id","ObjectNumber","pg_clusters_corrected_40")
  stopifnot(all(required%in%names(cd)),all(c("counts","exprs")%in%names(aa)))
  a<-aa$counts;b<-aa$exprs;panel<-unlist(config$physical_panel,use.names=FALSE)
  stopifnot(is.matrix(a),is.numeric(a),is.matrix(b),is.numeric(b),
    identical(dim(a),dim(b)),identical(rownames(a),panel),identical(rownames(b),panel),
    all(is.finite(a)),all(is.finite(b)))
  xy<-av(x,"spatialCoords")
  if(is.null(xy))xy<-dl(av(x,"int_colData"))$spatialCoords
  stopifnot(is.matrix(xy),is.numeric(xy),ncol(xy)==2L,nrow(xy)==ncol(a),
    identical(colnames(xy),c("Pos_X","Pos_Y")),all(lengths(cd[required])==ncol(a)))
  unit<-qualify_external_TDLN_units_RAM(cd$Patient.ID,cd$ROI.ID,cd$Image.name,
    cd$sample_id,cd$ObjectNumber,xy)
  if(!isTRUE(unit$registered_case_ROI_partition_gate_passed))return(list(
    status="EXTERNAL_UNIT_GATE_NOT_PASSED_NO_SPATIAL_EFFECTS",unit=unit,
    clinical_columns_not_selected=TRUE))
  derive<-function(z,rank_fun=function(v)rank(v,ties.method="average"))
    derive_external_TDLN_parents_states_RAM(z,panel,cd$pg_clusters_corrected_40,
      cd$Patient.ID,cd$ROI.ID,config$parent,rank_fun)
  ph<-derive(a);oracle<-derive(a,LN_audit_rank_RAM);alternative<-derive(b)
  stopifnot(identical(ph$internal$parents,oracle$internal$parents),
    identical(ph$internal$states,oracle$internal$states))
  agreement<-function(v,w)list(same_missingness=identical(is.na(v),is.na(w)),
    jointly_measured_rows=sum(!is.na(v)&!is.na(w)),
    discordant_jointly_measured_rows=sum(v!=w,na.rm=TRUE))
  sensitivity<-list(primary_assay="counts",comparison_assay="exprs",
    no_transformation_or_area_adjustment_added=TRUE,
    parent_masks=lapply(names(ph$internal$parents),function(k)
      c(list(parent=k),agreement(ph$internal$parents[[k]],alternative$internal$parents[[k]]))),
    state_masks=lapply(names(ph$internal$states),function(k)
      c(list(state=k),agreement(ph$internal$states[[k]],alternative$internal$states[[k]]))),
    alternative_assay_not_selected_by_effect_or_used_for_extra_hypothesis_search=TRUE)
  out<-list(status="EXTERNAL_PARENT_GATE_NOT_PASSED_NO_SPATIAL_EFFECTS",unit=unit,
    phenotype=ph$summary,independent_rank_masks_exactly_match=TRUE,assay_sensitivity=sensitivity,
    clinical_columns_not_selected=TRUE,author_dyad_or_distance_columns_not_selected=TRUE)
  if(!all(vapply(ph$internal$parents,any,logical(1))))return(out)
  spatial<-measure_external_TDLN_fixed_spatial_RAM(xy,cd$Patient.ID,cd$ROI.ID,ph,config$spatial)
  audit<-audit_external_TDLN_fixed_spatial_RAM(xy,spatial,config$spatial)
  stopifnot(audit$status=="PASS_FULL_EXTERNAL_FIXED_SPATIAL_NUMERIC_AUDIT")
  out$status<-"COMPLETE_EXTERNAL_FIXED_STATE_SPATIAL_EXTENSION_AUDITED"
  out$spatial<-spatial$summary;out$audit<-audit
  out$registered_case_keys_not_independent_patients_or_clinical_validation<-TRUE
  out$relative_construct_not_exact_RNA_CM_projection_or_direct_LR_mechanism<-TRUE
  out
}

if(Sys.getenv("PANCANCER_EXTERNAL_FIXED_MODE")=="NATIVE") {
  phase<-"READ_VERIFIED_INPUT_RAM_ONLY"
  result<-tryCatch({
    cfg<-jsonlite::fromJSON(Sys.getenv("PANCANCER_EXTERNAL_FIXED_CONFIG"),simplifyVector=FALSE)
    cfg$spatial$grids<-as.integer(unlist(cfg$spatial$grids,use.names=FALSE))
    cfg$spatial$absolute_deadline_UTC<-format(Sys.time()+cfg$maximum_analysis_seconds,tz="UTC",usetz=FALSE)
    x<-readRDS(file("stdin","rb"));phase<-"FINITE_UNIT_PARENT_STATE_FIXED_SPATIAL_AND_AUDIT"
    value<-external_TDLN_native_workflow_RAM(x,cfg)
    rm(x);gc()
    list(kind="SOURCE_ONLY_EXTERNAL_FIXED_WORKFLOW_COMPLETE",summary=value,native_RAM_released=TRUE)
  },error=function(e)list(kind="EXTERNAL_FIXED_WORKFLOW_NOT_ACCEPTED",phase=phase,error_class=class(e)[1]))
  cat(jsonlite::toJSON(result,auto_unbox=TRUE,null="null",na="null",digits=16),"\n")
  quit(status=if(result$kind=="SOURCE_ONLY_EXTERNAL_FIXED_WORKFLOW_COMPLETE")0L else 2L)
}
