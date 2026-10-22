# Base attributes only: metadata qualification without missing leaf-class code.
# All native identifiers, protein matrices and coordinates stay in this process.
# No clinical outcome values, author functions or external backing files used.
getslot <- function(x, name) attr(x, name, exact=TRUE)

list_data <- function(x) {
  if (is.null(x)) return(list())
  if (is.list(x) && !isS4(x)) return(x)
  value <- getslot(x, "listData")
  if (is.list(value)) return(value)
  list()
}

native_dims <- function(x) {
  d <- getslot(x, "dim")
  if (is.null(d)) d <- getslot(x, "Dim")
  if (length(d)==2L && is.numeric(d)) return(as.integer(d))
  NULL
}

marker_names <- function(x) {
  aliases <- list(CD3=c("cd3","cd3e"),CD4="cd4",CD8=c("cd8","cd8a"),
    FOXP3="foxp3",PD1=c("pd1","pdcd1"),Ki67=c("ki67","mki67"),
    GZMB=c("gzmb","granzymeb"),GZMK=c("gzmk","granzymek"),
    TCF1=c("tcf1","tcf7"),CCR7="ccr7",CD25=c("cd25","il2ra"),
    LAG3="lag3",TIM3=c("tim3","havcr2"),CTLA4="ctla4",CCR8="ccr8",
    CXCL13="cxcl13",LAMP3=c("lamp3","cd208"),CD1c="cd1c",CD141="cd141",
    CD68="cd68",CD163="cd163",CD20=c("cd20","ms4a1"),
    HLA_DR=c("hladr","hladra"),CD11c=c("cd11c","itgax"),
    CD45=c("cd45","ptprc"),CD45RA="cd45ra",CD45RO="cd45ro",TIA1="tia1")
  elements <- "(?:Nd|Sm|Gd|Dy|Er|Yb|Tb|Ho|Tm|Pr|Eu|Lu|Ce|La|In|Bi|Ir|Pt|Pd|Ag|Y)"
  clean <- gsub(paste0("^[0-9]{2,3}",elements,"[._ -]*"),"",x,perl=TRUE,ignore.case=TRUE)
  clean <- gsub(paste0("[._ -]*[0-9]{2,3}",elements,"$"),"",clean,perl=TRUE,ignore.case=TRUE)
  clean <- tolower(gsub("[^A-Za-z0-9]","",clean))
  lapply(names(aliases),function(k)list(marker=k,matching_feature_rows=sum(clean %in% aliases[[k]])))
}

qualify_structure <- function(x) {
  accepted <- c("SpatialExperiment","SingleCellExperiment","RangedSummarizedExperiment","SummarizedExperiment")
  classes <- as.character(getslot(x,"class"))
  cd <- getslot(x,"colData"); columns <- list_data(cd)
  n <- getslot(cd,"nrows")
  if (is.null(n) && is.data.frame(cd)) n <- .row_names_info(cd,2L)
  stopifnot(length(n)==1L, is.finite(n), n>0)
  a <- getslot(x,"assays"); aa <- list_data(getslot(a,"data"))
  if (!length(aa) && is.list(a)) aa <- a
  known <- c("counts","exprs","intensity","normcounts","normalized","scaled","arcsinh")
  assays <- lapply(seq_along(aa),function(i) {
    nm <- if(is.null(names(aa)))NA_character_ else names(aa)[i]; d <- native_dims(aa[[i]])
    list(position=i,recognized_assay_name=if(!is.na(nm)&&nm %in% known)nm else "UNREVIEWED_NAME",
         native_dimensions=d,external_or_lazy_assay_not_opened=is.null(d))
  })
  compatible <- which(vapply(aa,function(a){d<-native_dims(a);length(d)==2L&&d[2]==n&&d[1]<=256L},logical(1)))
  features <- NULL;feature_profiles<-list()
  if (length(compatible)) {
    a <- aa[[compatible[1]]]; d <- native_dims(a)
    dn <- getslot(a,"dimnames");if(is.null(dn))dn<-getslot(a,"Dimnames")
    if(length(dn)==2L)features<-dn[[1]]
    if(is.null(features))features<-getslot(x,"NAMES")
    if(is.null(features))features<-getslot(getslot(x,"rowRanges"),"NAMES")
    if(!is.null(features))stopifnot(length(features)==d[1])
    axes<-list(stored_assay_rownames=if(length(dn)==2L)dn[[1]] else NULL,
      top_level_names=getslot(x,"NAMES"),row_ranges_names=getslot(getslot(x,"rowRanges"),"NAMES"),
      row_ranges_partition_names=getslot(getslot(getslot(x,"rowRanges"),"partitioning"),"NAMES"))
    for(prefix in c("rowData","rowRanges")) {
      rd<-list_data(if(prefix=="rowData")getslot(x,"elementMetadata") else getslot(getslot(x,"rowRanges"),"elementMetadata"))
      for(k in intersect(c("marker","marker_name","Marker","antigen","Antigen","symbol","gene_symbol","target","Target"),names(rd)))
        axes[[paste0(prefix,"_",k)]]<-rd[[k]]
    }
    for(k in names(axes))if(!is.null(axes[[k]])&&length(axes[[k]])==d[1]&&(is.character(axes[[k]])||is.factor(axes[[k]])))
      feature_profiles[[k]]<-marker_names(as.character(axes[[k]]))
  }
  spatial <- getslot(x,"spatialCoords")
  if(is.null(spatial))spatial<-list_data(getslot(x,"int_colData"))[["spatialCoords"]]
  dims <- native_dims(spatial)
  coords <- list(native_dimensions=dims,coordinates_not_exported=TRUE,
                 physical_units_and_ROI_identity_not_accepted=TRUE)
  if(length(dims)==2L&&dims[1]==n&&dims[2]>=2L&&is.numeric(spatial)) {
    coords$first_two_axes_finite_rows<-sum(rowSums(is.finite(spatial[,1:2,drop=FALSE]))==2L)
  }
  idfields<-c("patient_id","patient","Patient","donor_id","donor","subject_id",
    "sample_id","sample","sampleID","image_id","image","ROI","roi_id","slide_id","slide","ImageNumber")
  semantic_key_names<-names(columns)[grepl("^(patient|donor|subject|case|sample|image|slide|roi|fov)([_. ]?(id|number|key|name))?$",names(columns),ignore.case=TRUE)]
  units<-lapply(union(intersect(idfields,names(columns)),semantic_key_names),function(k) {
    v<-columns[[k]];stopifnot(length(v)==n)
    v<-as.character(v);ok<-!is.na(v)&!tolower(trimws(v))%in%c("","na","nan","none","unknown")
    list(field=k,nonmissing_cells=sum(ok),distinct_literal_keys=length(unique(v[ok])),
         individual_values_not_exported=TRUE,patient_semantics_not_inferred=TRUE)
  })
  schemafields<-c("cell_type","celltype","CellType","annotation","phenotype","cluster",
    "stage","Stage","tumor_stage","clinical_stage","LN_station","station","location","tissue","timepoint","treatment")
  list(status="COMPLETE_RAM_ONLY_STRUCTURE_AND_PREDECLARED_MARKER_QUALIFICATION_NO_EFFECTS",
    known_object_classes=intersect(accepted,classes),unrecognized_class_count=sum(!classes%in%accepted),
    source_cells=as.integer(n),assays=assays,native_protein_axis_candidates=length(compatible),
    predeclared_marker_matches=if(!is.null(features))marker_names(features) else NULL,
    all_predeclared_feature_annotation_axis_profiles=feature_profiles,
    semantic_marker_axis_not_selected_or_adopted=TRUE,
    feature_identity_unaccepted=is.null(features),spatial=coords,unit_key_schema_and_counts=units,
    known_metadata_fields_present=intersect(schemafields,names(columns)),
    other_colData_fields_not_selected=sum(!names(columns)%in%c(idfields,schemafields,semantic_key_names)),
    clinical_outcome_values_not_selected=TRUE,assay_intensity_values_not_exported_or_interpreted=TRUE,
    original_cell_or_patient_identifiers_not_returned=TRUE,identity_patient_spatial_or_clinical_effect_not_accepted=TRUE)
}

synthetic <- function() {
  suppressPackageStartupMessages(library(SingleCellExperiment))
  methods::setClass("SyntheticSpatialRAM",contains="SingleCellExperiment",slots=c(spatialCoords="matrix"))
  m<-matrix(seq_len(60),6,10,dimnames=list(c("CD3","CD4","CD8a","FoxP3_155Gd","PD-1","Ki67"),NULL))
  sce<-SingleCellExperiment(assays=list(counts=m),colData=S4Vectors::DataFrame(patient_id=rep(c("FAKE_A","FAKE_B"),each=5)))
  unnamed<-sce;a<-getslot(unnamed,"assays");ad<-getslot(a,"data");ld<-getslot(ad,"listData")
  names(ld)<-NULL;attr(ad,"listData")<-ld;attr(a,"data")<-ad;attr(unnamed,"assays")<-a
  stopifnot(qualify_structure(unnamed)$assays[[1]]$recognized_assay_name=="UNREVIEWED_NAME")
  isotope<-sce;SummarizedExperiment::rowData(isotope)<-S4Vectors::DataFrame(marker=rownames(m))
  rownames(isotope)<-c("170Er","154Sm","146Nd","155Gd","168Er","162Dy")
  qi<-qualify_structure(isotope);profiles<-qi$all_predeclared_feature_annotation_axis_profiles
  stopifnot(length(profiles)>=2L,any(vapply(profiles,function(z)sum(vapply(z,function(v)v$matching_feature_rows,numeric(1)))==0,logical(1))),
            any(vapply(profiles,function(z)sum(vapply(z,function(v)v$matching_feature_rows,numeric(1)))==6,logical(1))))
  x<-methods::new("SyntheticSpatialRAM",sce,spatialCoords=matrix(seq_len(20),10,2))
  raw<-serialize(x,NULL);methods::removeClass("SyntheticSpatialRAM");x<-unserialize(raw)
  q<-qualify_structure(x)
  stopifnot(q$source_cells==10L,q$spatial$first_two_axes_finite_rows==10L,
            q$unit_key_schema_and_counts[[1]]$distinct_literal_keys==2L)
  marker<-setNames(vapply(q$predeclared_marker_matches,function(v)v$matching_feature_rows,numeric(1)),
                   vapply(q$predeclared_marker_matches,function(v)v$marker,character(1)))
  stopifnot(all(marker[c("CD3","CD4","CD8","FOXP3","PD1","Ki67")]==1L))
  stopifnot(!grepl("FAKE_A",jsonlite::toJSON(q,auto_unbox=TRUE)))
  TRUE
}
