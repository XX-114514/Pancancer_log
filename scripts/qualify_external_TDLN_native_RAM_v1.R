# Structural and physical-marker qualification only. No states or effects.
# Companion structure/schema functions are bound by the execution manifest.
external_marker_aliases <- list(
  CD3=c("cd3","cd3e"), CD4="cd4", CD8=c("cd8","cd8a"),
  FOXP3="foxp3", CTLA4=c("ctla4","cd152"), PD1=c("pd1","pdcd1","cd279"),
  TIM3=c("tim3","havcr2","cd366"), MKI67=c("mki67","ki67"),
  CCR7=c("ccr7","cd197"), IL7R=c("il7r","cd127"),
  PLA2G2D="pla2g2d", CD86="cd86", PDL1=c("pdl1","cd274","cd274pdl1"),
  CD28="cd28", TBX21=c("tbx21","tbet"), IFNG=c("ifng","ifngamma"),
  CD68="cd68", CD163="cd163", CD11c=c("cd11c","itgax"),
  LAMP3=c("lamp3","cd208"), TCF7=c("tcf7","tcf1"), GZMB=c("gzmb","granzymeb"),
  GZMK=c("gzmk","granzymek"), CD45=c("cd45","ptprc"), CD19="cd19", CD20=c("cd20","ms4a1"))

normalize_external_marker <- function(v) {
  elements <- "(?:Nd|Sm|Gd|Dy|Er|Yb|Tb|Ho|Tm|Pr|Eu|Lu|Ce|La|In|Bi|Ir|Pt|Pd|Ag|Y)"
  v<-gsub(paste0("^[0-9]{2,3}",elements,"[._ -]*"),"",v,perl=TRUE,ignore.case=TRUE)
  v<-gsub(paste0("[._ -]*[0-9]{2,3}",elements,"$"),"",v,perl=TRUE,ignore.case=TRUE)
  tolower(gsub("[^A-Za-z0-9]","",v))
}

qualify_external_TDLN <- function(x) {
  structural<-qualify_structure(x)
  schema<-qualify_LN_native_metadata_schema_RAM(x)
  aa<-list_data(getslot(getslot(x,"assays"),"data"))
  axes<-list(top_level=getslot(x,"NAMES"), row_ranges=getslot(getslot(x,"rowRanges"),"NAMES"),
    partition=getslot(getslot(getslot(x,"rowRanges"),"partitioning"),"NAMES"))
  for(prefix in c("rowData","rowRanges")) {
    rd<-list_data(if(prefix=="rowData")getslot(x,"elementMetadata") else getslot(getslot(x,"rowRanges"),"elementMetadata"))
    for(k in intersect(c("marker","marker_name","Marker","antigen","Antigen","symbol","gene_symbol","target","Target","channel_name","channel"),names(rd)))
      axes[[paste0(prefix,"_",k)]]<-rd[[k]]
  }
  profiles<-list()
  for(i in seq_along(aa)) {
    a<-aa[[i]];d<-native_dims(a)
    if(length(d)!=2L||d[1]>256L||d[2]!=structural$source_cells)next
    dn<-getslot(a,"dimnames");if(is.null(dn))dn<-getslot(a,"Dimnames")
    candidates<-axes
    if(length(dn)==2L)candidates$stored_assay<-dn[[1]]
    for(k in names(candidates)) {
      v<-candidates[[k]]
      if(length(v)!=d[1]||!(is.character(v)||is.factor(v)))next
      cleaned<-normalize_external_marker(as.character(v))
      markers<-lapply(names(external_marker_aliases),function(marker) {
        ii<-which(cleaned%in%external_marker_aliases[[marker]])
        z<-list(marker=marker,matching_rows=length(ii),declared_axis_only=TRUE)
        if(length(ii)==1L&&is.matrix(a)&&is.numeric(a)) {
          values<-a[ii,,drop=TRUE];finite<-is.finite(values)
          z$numeric_RAM_assay_checked<-TRUE
          z$finite_cells<-sum(finite);z$nonfinite_cells<-sum(!finite)
          z$nonzero_finite_cells<-sum(values[finite]!=0)
          z$finite_values_vary<-length(unique(values[finite]))>1L
          z$does_not_establish_positive_state_or_antibody_specificity<-TRUE
        }
        z
      })
      profiles[[length(profiles)+1L]]<-list(assay_position=i,axis=k,markers=markers)
    }
  }
  list(status="EXTERNAL_TDLN_STRUCTURE_MARKER_NUMERIC_QUALIFICATION_ONLY",
    structural=structural,schema=schema,predeclared_marker_numeric_profiles=profiles,
    no_states_thresholds_graphs_effects_or_clinical_values_computed=TRUE,
    missing_or_ambiguous_marker_not_filled_as_negative=TRUE,
    patient_independence_region_truth_and_endpoint_semantics_not_accepted=TRUE)
}

external_synthetic <- function() {
  suppressPackageStartupMessages(library(SingleCellExperiment))
  m<-matrix(seq_len(48),4,12,dimnames=list(c("CD4","FoxP3_155Gd","CTLA4","PLA2G2D"),NULL))
  m[3,]<-0
  x<-SingleCellExperiment(assays=list(counts=m),colData=S4Vectors::DataFrame(patient_id=rep(c("INDIVIDUAL_CANARY_A","INDIVIDUAL_CANARY_B"),each=6)))
  out<-qualify_external_TDLN(x)
  stopifnot(out$structural$source_cells==12L,!grepl("INDIVIDUAL_CANARY",jsonlite::toJSON(out,auto_unbox=TRUE)))
  profiles<-out$predeclared_marker_numeric_profiles
  stopifnot(length(profiles)>0L)
  one<-profiles[[length(profiles)]]$markers
  named<-setNames(one,vapply(one,function(z)z$marker,character(1)))
  stopifnot(named$CTLA4$matching_rows==1L,!named$CTLA4$finite_values_vary,
    named$CTLA4$nonzero_finite_cells==0L,named$TIM3$matching_rows==0L,
    named$PLA2G2D$matching_rows==1L)
  rownames(x)<-c("CD4","FOXP3","CTLA4","CD152")
  duplicate<-qualify_external_TDLN(x)$predeclared_marker_numeric_profiles
  stopifnot(any(vapply(duplicate,function(z)any(vapply(z$markers,function(y)y$marker=="CTLA4"&&y$matching_rows==2L&&is.null(y$finite_cells),logical(1))),logical(1))))
  TRUE
}

if(Sys.getenv("PANCANCER_EXTERNAL_TDLN_MODE")=="CHECK") {
  stopifnot(external_synthetic())
  cat('{"status":"PASS_EXTERNAL_TDLN_SYNTHETIC_QUALIFICATION_PRIVACY_AND_MARKER_GATES"}\n')
} else if(Sys.getenv("PANCANCER_EXTERNAL_TDLN_MODE")=="NATIVE") {
  phase<-"READ_VERIFIED_OPEN_DATASET_IN_RAM"
  answer<-tryCatch({
    native_object<-readRDS(file("stdin","rb"))
    phase<-"FINITE_NATIVE_STRUCTURAL_AND_MARKER_QUALIFICATION"
    result<-qualify_external_TDLN(native_object)
    rm(native_object);gc()
    list(kind="SOURCE_ONLY_EXTERNAL_TDLN_QUALIFICATION",summary=result,native_RAM_released=TRUE)
  },error=function(e)list(kind="EXTERNAL_TDLN_NATIVE_NOT_ACCEPTED",phase=phase,error_class=class(e)[1]))
  cat(jsonlite::toJSON(answer,auto_unbox=TRUE,null="null",digits=16),"\n")
  quit(status=if(answer$kind=="SOURCE_ONLY_EXTERNAL_TDLN_QUALIFICATION")0L else 2L)
}
