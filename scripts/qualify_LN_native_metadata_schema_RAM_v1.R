# Schema and physical antibody panel only; no metadata values or native arrays leave RAM.
qualify_LN_native_metadata_schema_RAM <- function(x) {
  slot_value <- function(v, k) attr(v, k, exact=TRUE)
  data_list <- function(v) {
    if(is.null(v))return(list())
    if(is.list(v)&&!isS4(v))return(v)
    z<-slot_value(v,"listData");if(is.list(z))z else list()
  }
  safe_names <- function(v) {
    z<-names(v);if(is.null(z))return(character())
    ok<-!is.na(z)&nchar(z)<=100L&!(grepl(intToUtf8(10L),z,fixed=TRUE)|grepl(intToUtf8(13L),z,fixed=TRUE))
    z[!ok]<-"UNREVIEWED_FIELD_NAME"
    z
  }
  schema <- function(v) lapply(seq_along(v),function(i) {
    cl<-class(v[[i]]);allowed<-c("character","factor","ordered","integer","numeric","logical","matrix","array","data.frame","list","DFrame","SimpleList","Rle")
    list(position=i,field=safe_names(v)[i],known_classes=intersect(cl,allowed),
      unrecognized_class_count=sum(!cl%in%allowed),length=length(v[[i]]),values_and_factor_levels_not_selected=TRUE)
  })
  aa<-data_list(slot_value(slot_value(x,"assays"),"data"))
  stopifnot(length(aa)>0L)
  assay_names<-safe_names(aa)
  assay_schema<-lapply(seq_along(aa),function(i) {
    a<-aa[[i]];d<-slot_value(a,"dim");if(is.null(d))d<-slot_value(a,"Dim")
    cl<-class(a);allowed<-c("matrix","array","dgCMatrix","dgeMatrix")
    list(position=i,name=assay_names[i],known_classes=intersect(cl,allowed),
      unrecognized_class_count=sum(!cl%in%allowed),numeric_in_RAM=is.numeric(a),dimensions=d)
  })
  panels<-lapply(seq_along(aa),function(i) {
    a<-aa[[i]];dn<-slot_value(a,"dimnames");if(is.null(dn))dn<-slot_value(a,"Dimnames")
    z<-if(length(dn)==2L)dn[[1]] else NULL
    if(is.null(z))return(list(assay_position=i,panel_present=FALSE))
    if(!is.character(z)||length(z)>256L)return(list(assay_position=i,panel_present=FALSE,panel_not_accepted=TRUE))
    ok<-!is.na(z)&nchar(z)<=100L&!(grepl(intToUtf8(10L),z,fixed=TRUE)|grepl(intToUtf8(13L),z,fixed=TRUE));z[!ok]<-"UNREVIEWED_PHYSICAL_FEATURE_LABEL"
    list(assay_position=i,panel_present=TRUE,physical_feature_labels=as.list(z),
      panel_labels_not_learned_program_components=TRUE)
  })
  cd<-data_list(slot_value(x,"colData"))
  rd<-data_list(slot_value(slot_value(x,"rowRanges"),"elementMetadata"))
  md<-slot_value(x,"metadata");if(is.null(md))md<-list()
  spatial<-slot_value(x,"spatialCoords")
  if(is.null(spatial))spatial<-data_list(slot_value(x,"int_colData"))[["spatialCoords"]]
  dn<-slot_value(spatial,"dimnames");axis_names<-if(length(dn)==2L)dn[[2]] else NULL
  if(!is.null(axis_names)) {ok<-!is.na(axis_names)&nchar(axis_names)<=100L&!(grepl(intToUtf8(10L),axis_names,fixed=TRUE)|grepl(intToUtf8(13L),axis_names,fixed=TRUE));axis_names[!ok]<-"UNREVIEWED_AXIS_NAME"}
  list(status="NATIVE_SCHEMA_AND_PHYSICAL_PANEL_ONLY_NO_VALUES_OR_EFFECTS",
    assays=assay_schema,physical_panels=panels,colData_schema=schema(cd),
    row_annotation_schema=schema(rd),top_metadata_field_names=as.list(safe_names(md)),
    spatial_axis_names=axis_names,individual_values_factor_levels_and_clinical_outcomes_not_selected=TRUE,
    assay_values_and_coordinates_not_returned=TRUE,external_or_lazy_backing_not_opened=TRUE)
}

if(exists("native_object",inherits=FALSE)) {
  result<-tryCatch(qualify_LN_native_metadata_schema_RAM(native_object),error=function(e)list(status="NATIVE_SCHEMA_STAGE_NOT_ACCEPTED_NO_NATIVE_VALUES_RETURNED",error_class=class(e)[1]))
  result$observed_UTC<-format(Sys.time(),tz="UTC",usetz=TRUE)
  writeLines(jsonlite::toJSON(result,auto_unbox=TRUE,null="null",digits=16),
    file.path(control_root,"private_native_metadata_schema_phase001.json"))
}
