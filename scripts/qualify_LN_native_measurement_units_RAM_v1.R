# Requires the separately SHA-bound marker and geometry helpers. RAM only.
qualify_LN_native_measurement_units_RAM <- function(x, marker_rows, expected_panel, fields) {
  av<-function(x,k)attr(x,k,exact=TRUE)
  dl<-function(x){if(is.list(x)&&!isS4(x))return(x);v<-av(x,"listData");if(is.list(v))v else list()}
  assays<-dl(av(av(x,"assays"),"data"));cd<-dl(av(x,"colData"))
  stopifnot(all(c("counts","asin")%in%names(assays)),all(unlist(fields)%in%names(cd)))
  a<-assays[["counts"]];b<-assays[["asin"]]
  stopifnot(is.matrix(a),is.numeric(a),is.matrix(b),is.numeric(b),identical(dim(a),dim(b)),
    identical(rownames(a),expected_panel),identical(rownames(b),expected_panel))
  n<-ncol(a);v<-lapply(fields,function(k){z<-cd[[k]];stopifnot((is.character(z)||is.factor(z)),length(z)==n);as.character(z)})
  valid<-function(z)!is.na(z)&!tolower(trimws(z))%in%c("","na","nan","none","null","unknown","n/a")
  distinct<-function(z)length(unique(z[valid(z)]))
  nesting<-function(child,parent) {
    ok<-valid(child)&valid(parent);edges<-unique(data.frame(child=child[ok],parent=parent[ok]))
    t<-table(edges$child)
    list(complete_rows=sum(ok),child_labels=length(t),child_labels_with_multiple_parents=sum(t>1L),
      maximum_parents_per_child=if(length(t))as.integer(max(t)) else NULL)
  }
  original_to_corrected<-nesting(v$original_patient,v$corrected_patient)
  corrected_to_original<-nesting(v$corrected_patient,v$original_patient)
  equal_complete<-valid(v$original_patient)&valid(v$corrected_patient)
  patient_partitions<-list(original_distinct_labels=distinct(v$original_patient),corrected_distinct_labels=distinct(v$corrected_patient),
    jointly_complete_rows=sum(equal_complete),different_literal_label_rows=sum(v$original_patient[equal_complete]!=v$corrected_patient[equal_complete]),
    original_to_corrected=original_to_corrected,corrected_to_original=corrected_to_original,
    same_partition_on_jointly_complete_rows=sum(equal_complete)>0L&&original_to_corrected$child_labels_with_multiple_parents==0L&&corrected_to_original$child_labels_with_multiple_parents==0L,
    literal_partition_not_independent_patient_review=TRUE)
  hierarchy<-list(node_to_original_patient=nesting(v$node,v$original_patient),node_to_corrected_patient=nesting(v$node,v$corrected_patient),
    image_to_node=nesting(v$image,v$node),image_to_original_patient=nesting(v$image,v$original_patient),
    image_to_corrected_patient=nesting(v$image,v$corrected_patient),image_to_raw_ROI=nesting(v$image,v$ROI),
    sample_to_image=nesting(v$sample,v$image),image_to_sample=nesting(v$image,v$sample))
  complete<-valid(v$corrected_patient)&valid(v$node)&valid(v$image)
  roi_key<-unit_key_RAM(list(v$corrected_patient,v$node,v$image));roi_key[!complete]<-NA_character_
  numeric_profiles<-list(counts=profile_LN_marker_assay_RAM(a,marker_rows,roi_key),asin=profile_LN_marker_assay_RAM(b,marker_rows,roi_key))
  transform<-lapply(c(1,5,10,20,100),function(cofactor) {
    total<-nonfinite<-violations<-0;maxerr<-NA_real_
    for(j in seq_len(nrow(a))) {
      raw<-as.numeric(a[j,]);target<-as.numeric(b[j,]);ok<-is.finite(raw)&is.finite(target)
      total<-total+sum(ok);nonfinite<-nonfinite+sum(!ok)
      if(any(ok)) {
        predicted<-asinh(raw[ok]/cofactor);err<-abs(predicted-target[ok]);violations<-violations+sum(err>1e-8*pmax(1,abs(target[ok])))
        maxerr<-max(c(maxerr,err),na.rm=TRUE)
      }
    }
    list(candidate_function="asinh(counts/cofactor)",cofactor=cofactor,finite_pairs=total,nonfinite_pairs=nonfinite,
      tolerance_violations=violations,max_absolute_error=maxerr,empirically_equal_at_predeclared_tolerance=total>0&&violations==0)
  })
  type_records<-function(z,parent=NULL) {
    labels<-sort(unique(z[valid(z)]));if(length(labels)>128L)return(list(status="TOO_MANY_LABELS_NOT_RELEASED"))
    records<-lapply(labels,function(label) {
      if(nchar(label)>120L||grepl(intToUtf8(10L),label,fixed=TRUE)||grepl(intToUtf8(13L),label,fixed=TRUE))return(list(label="UNREVIEWED_TYPE_LABEL"))
      i<-which(!is.na(z)&z==label)
      list(author_type=label,cell_rows=length(i),candidate_corrected_patient_labels=distinct(v$corrected_patient[i]),
        source_image_labels=distinct(v$image[i]),coarse_parent_labels=if(is.null(parent))NULL else as.list(sort(unique(parent[i][valid(parent[i])]))))
    })
    list(status="AUTHOR_ANNOTATION_SOURCE_AGGREGATES_ONLY",missing_rows=sum(!valid(z)),records=records)
  }
  category_summary<-function(z,allowed,patient) {
    z<-toupper(trimws(z));z<-sub("^STAGE[ ]+","",z)
    labels<-c(allowed,"MISSING","UNREVIEWED_NONMISSING");code<-ifelse(!valid(z),"MISSING",ifelse(z%in%allowed,z,"UNREVIEWED_NONMISSING"))
    lapply(labels,function(k){i<-which(code==k);list(category=k,cell_rows=length(i),candidate_patient_labels=distinct(patient[i]),
      patient_node_tuples={ok<-valid(patient[i])&valid(v$node[i]);length(unique(unit_key_RAM(list(patient[i][ok],v$node[i][ok]))))},source_image_labels=distinct(v$image[i]))})
  }
  paired_count<-function(patient) {
    ok<-valid(patient)&toupper(trimws(v$station))%in%c("N1","N2")
    groups<-split(toupper(trimws(v$station[ok])),patient[ok])
    sum(vapply(groups,function(z)all(c("N1","N2")%in%z),logical(1)))
  }
  coordinates<-av(x,"spatialCoords");if(is.null(coordinates))coordinates<-dl(av(x,"int_colData"))[["spatialCoords"]]
  node_key<-unit_key_RAM(list(v$corrected_patient,v$node));node_key[!valid(v$corrected_patient)|!valid(v$node)]<-NA_character_
  geometry<-qualify_LN_unit_geometry_RAM(coordinates[,1:2,drop=FALSE],data.frame(patient=v$corrected_patient,node=v$node,ROI=v$image,stringsAsFactors=FALSE))
  list(status="NATIVE_NUMERIC_AUTHOR_TYPE_AND_UNIT_QUALIFICATION_NO_STATE_OR_EFFECT",cells=n,
    selected_schema_fields=fields,patient_partition_checks=patient_partitions,unit_hierarchy=hierarchy,
    candidate_corrected_patient_node_tuples=length(unique(unit_key_RAM(list(v$corrected_patient[complete],v$node[complete])))),
    numeric_profiles=numeric_profiles,predeclared_transform_checks=transform,
    author_coarse_types=type_records(v$coarse_type),author_fine_types=type_records(v$fine_type,v$coarse_type),
    reported_stage_categories=category_summary(v$stage,c("I","IA","IB","II","IIA","IIB","III","IIIA","IIIB","IIIC","IV","IVA","IVB"),v$corrected_patient),
    reported_station_categories=category_summary(v$station,c("N1","N2"),v$corrected_patient),
    original_patient_labels_with_both_stations=paired_count(v$original_patient),corrected_patient_labels_with_both_stations=paired_count(v$corrected_patient),
    stage_to_patient_consistency=nesting(v$corrected_patient,v$stage),station_to_node_consistency=nesting(node_key,v$station),
    geometry_under_candidate_corrected_patient_node_image_mapping=geometry,
    geometry_grouping="literal corrected-patient x node x image; raw ROI label is not assumed globally unique",
    physical_units_and_patient_semantics_not_inferred_from_names=TRUE,
    composition_derived_neighbor_matrix_and_niches_not_selected=TRUE,
    no_new_state_thresholds_or_spatial_or_clinical_effects=TRUE,
    no_individual_identifiers_records_coordinates_or_omics_arrays_returned=TRUE)
}
