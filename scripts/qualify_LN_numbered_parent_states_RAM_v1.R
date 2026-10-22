# Protein-relative candidate annotation; no spatial association, no expression matrix export.
qualify_LN_numbered_parent_states_RAM <- function(x, expected_panel, config) {
  av<-function(x,k)attr(x,k,exact=TRUE)
  dl<-function(x){if(is.list(x)&&!isS4(x))return(x);v<-av(x,"listData");if(is.list(v))v else list()}
  aa<-dl(av(av(x,"assays"),"data"));cd<-dl(av(x,"colData"));rd<-dl(av(av(x,"rowRanges"),"elementMetadata"))
  assay<-aa[["asin"]];stopifnot(is.matrix(assay),identical(rownames(assay),expected_panel))
  n<-ncol(assay);fine<-as.character(cd[[config$fine_field]]);image<-as.character(cd[[config$image_field]]);patient<-as.character(cd[[config$patient_field]])
  valid<-function(z)!is.na(z)&!tolower(trimws(z))%in%c("","na","unknown","none","null")
  good<-valid(fine)&valid(image)&valid(patient);flevels<-sort(unique(fine[good]));ilevels<-unique(image[good]);plevels<-unique(patient[good])
  fi<-match(fine,flevels);ii<-match(image,ilevels);pi<-match(patient,plevels)
  gid<-(ii-1L)*length(flevels)+fi;groups<-split(which(good),gid[good]);group_index<-as.integer(names(groups));group_fine<-(group_index-1L)%%length(flevels)+1L
  group_image<-(group_index-1L)%/%length(flevels)+1L;group_n<-lengths(groups);eligible_group<-group_n>=config$minimum_cells_per_image_group
  image_groups<-split(which(good),ii[good]);image_patient<-pi[match(ilevels,image)]
  marker_positions<-match(config$qualification_markers,expected_panel)
  stopifnot(!anyNA(marker_positions));ranks<-matrix(NA_real_,length(marker_positions),length(groups),dimnames=list(config$qualification_markers,NULL))
  marker_flags<-lapply(seq_along(marker_positions),function(j) {
    row<-marker_positions[j];z<-as.numeric(assay[row,]);rz<-rep(NA_real_,n)
    for(i in image_groups){ok<-is.finite(z[i]);if(any(ok))rz[i[ok]]<-rank(z[i[ok]],ties.method="average")/(sum(ok)+1)}
    for(g in which(eligible_group)){a<-rz[groups[[g]]];a<-a[is.finite(a)];if(length(a)>=config$minimum_cells_per_image_group)ranks[j,g]<<-median(a)}
    flag<-function(k){v<-rd[[k]];if(is.null(v)||length(v)!=length(expected_panel))return("MISSING");w<-v[row];if(length(w)!=1||!is.finite(w)||!w%in%c(0,1))"UNREVIEWED" else if(w==1)"ONE" else "ZERO"}
    list(marker=config$qualification_markers[j],finite_rows=sum(is.finite(z)),nonconstant=length(unique(z[is.finite(z)]))>1L,
      binary_row_field_keep=flag("keep"),binary_row_field_ilastik=flag("ilastik"),binary_row_field_deepcell=flag("deepcell"))
  })
  candidate<-function(r) {
    take<-function(k){v<-r[k];if(length(v)!=1L||!is.finite(v))return(NA_real_);v}
    if(any(!is.finite(r)))return(character())
    myeloid<-mean(r[c("CD14","CD68","CD163")]);dc<-mean(r[c("CD11c","CD1c","CD141","HLADRB")]);opponent<-max(myeloid,r[c("CD19","CD56","EpCAM","SMA")])
    anchor<-max(r[c("CD45RA","CD45RO")])>=.55
    flags<-c(CD4_Tlike=r["CD4"]>=.70&&r["CD8a"]<=.60&&opponent<=.65&&r["CD4"]-opponent>=.10&&anchor,
      CD8a_Tlike=r["CD8a"]>=.70&&r["CD4"]<=.60&&opponent<=.65&&r["CD8a"]-opponent>=.10&&anchor,
      B_like=r["CD19"]>=.70&&r["CD19"]-max(myeloid,r[c("CD4","CD8a","CD56","EpCAM","SMA")])>=.10,
      Myeloid_like=myeloid>=.65&&myeloid-max(r[c("CD19","CD4","CD8a","EpCAM","SMA")])>=.05,
      DC_like=dc>=.65&&dc-max(r[c("CD19","CD4","CD8a","EpCAM","SMA")])>=.05,
      NK_like=r["CD56"]>=.70&&r["CD56"]-max(myeloid,r[c("CD19","CD4","EpCAM","SMA")])>=.10,
      Epithelial_like=r["EpCAM"]>=.70&&r["EpCAM"]-max(myeloid,r[c("CD19","CD4","CD8a","CD56","SMA")])>=.10,
      Stromal_like=r["SMA"]>=.70&&r["SMA"]-max(myeloid,r[c("CD19","CD4","CD8a","CD56","EpCAM")])>=.10)
    # Scalar logical operations drop input names; class keys are fixed.
    names(flags)<-c("CD4_Tlike","CD8a_Tlike","B_like","Myeloid_like","DC_like","NK_like","Epithelial_like","Stromal_like")
    names(flags)[flags]
  }
  profiles<-matrix(NA_real_,nrow(ranks),length(flevels),dimnames=list(rownames(ranks),flevels));records<-vector("list",length(flevels));group_classes<-lapply(seq_len(ncol(ranks)),function(g)candidate(ranks[,g]))
  for(k in seq_along(flevels)) {
    g<-which(group_fine==k&eligible_group);cases<-unique(image_patient[group_image[g]])
    for(j in seq_len(nrow(ranks))) {
      per_case<-vapply(cases,function(p){z<-ranks[j,g[image_patient[group_image[g]]==p]];z<-z[is.finite(z)];if(length(z))median(z) else NA_real_},numeric(1));per_case<-per_case[is.finite(per_case)]
      if(length(per_case))profiles[j,k]<-median(per_case)
    }
    classes<-candidate(profiles[,k]);support<-vapply(classes,function(c)if(length(g))mean(vapply(group_classes[g],function(z)c%in%z,logical(1))) else 0,numeric(1))
    eligible<-length(g)>=config$minimum_images&&length(cases)>=config$minimum_patients&&sum(good&fine==flevels[k],na.rm=TRUE)>=config$minimum_cells
    accepted<-classes[support>=config$minimum_image_agreement];if(!eligible)accepted<-character()
    r<-profiles[,k];high<-function(markers,cutoffs)all(is.finite(r[markers]))&&all(r[markers]>=cutoffs)
    states<-list(CD4_FOXP3_CTLA4_enriched="CD4_Tlike"%in%accepted&&high(c("FoxP3","CTLA4"),c(.70,.65)),
      MKI67_enriched=length(accepted)>0&&high("MKI67",.70),
      CD8_PD1_TIM3_enriched="CD8a_Tlike"%in%accepted&&high(c("PD1","TIM3"),c(.65,.65)),
      CD8_CCR7_CD127_enriched="CD8a_Tlike"%in%accepted&&high(c("CCR7","CD127"),c(.65,.65)))
    records[[k]]<-list(author_fine_group=flevels[k],candidate_lineages=as.list(accepted),ambiguous_lineage=length(accepted)>1L,
      cell_rows=sum(good&fine==flevels[k],na.rm=TRUE),eligible_image_groups=length(g),candidate_patient_labels=length(cases),
      minimum_image_agreement_satisfied=length(accepted)>0,source_support_fraction=if(length(accepted))as.list(support[accepted]) else NULL,
      predeclared_state_support_flags=states,relative_marker_candidate_not_formal_identity=TRUE)
  }
  stage<-as.character(cd[[config$stage_field]]);values<-sort(unique(stage[valid(stage)]))
  allowed_stage<-function(s){z<-tolower(trimws(s));z<-sub("^stage[ _-]*","",z);grepl("^[1-4iv]+[abc]?([ _/-]+[1-4iv]+[abc]?)?$",z)||z%in%c("early","late","advanced","intermediate","low","high","normal","healthy","benign","early-stage","late-stage","early_stage","late_stage","ais","mia")}
  stage_categories<-lapply(values,function(s){ok<-nchar(s)<=40L&&allowed_stage(s);list(source_stage_code=if(ok)s else "UNREVIEWED_CODE_NOT_DISCLOSED",cell_rows=sum(stage==s,na.rm=TRUE),not_interpreted_as_AJCC_or_pathologic_N=TRUE)})
  list(summary=list(status="NUMBERED_PARENT_RELATIVE_PROTEIN_QUALIFICATION_NO_SPATIAL_EFFECTS",source_fine_group_count=length(flevels),
    candidate_records=records,marker_measurement_flags=marker_flags,source_stage_category_codes=stage_categories,
    rank_profiles_and_individual_data_not_exported=TRUE,coordinates_neighbor_composition_and_clinical_outcomes_not_selected=TRUE,
    classification_thresholds_are_relative_not_absolute_protein_positivity=TRUE),
    internal=list(fine_labels=flevels,records=records,relative_rank_profiles=profiles,marker_flags=marker_flags))
}
