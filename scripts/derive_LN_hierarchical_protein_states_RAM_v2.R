# One prospective hierarchical alternative. All individual vectors stay RAM only.
LN_rank_profiles_RAM <- function(a, markers, fine, image, patient, take,
  minimum_cells=20L,rank_fun=function(z)rank(z,ties.method="average")) {
  f<-sort(unique(fine[take]));im<-unique(image[take]);fi<-match(fine,f);ii<-match(image,im)
  gids<-(ii-1L)*length(f)+fi;groups<-split(which(take),gids[take]);keys<-as.integer(names(groups));gf<-(keys-1L)%%length(f)+1L;gi<-(keys-1L)%/%length(f)+1L
  eligible<-lengths(groups)>=minimum_cells;imgs<-split(which(take),ii[take]);case<-patient[match(im,image)]
  gprofile<-matrix(NA_real_,length(markers),length(groups),dimnames=list(markers,NULL))
  for(j in seq_along(markers)) {
    z<-as.numeric(a[markers[j],]);rz<-rep(NA_real_,length(z))
    for(i in imgs){ok<-is.finite(z[i]);if(any(ok))rz[i[ok]]<-rank_fun(z[i[ok]])/(sum(ok)+1)}
    gprofile[j,]<-vapply(seq_along(groups),function(k){v<-rz[groups[[k]]];v<-v[is.finite(v)];if(eligible[k]&&length(v)>=minimum_cells)median(v) else NA_real_},numeric(1))
  }
  profiles<-matrix(NA_real_,length(markers),length(f),dimnames=list(markers,f))
  for(k in seq_along(f))for(j in seq_along(markers)) {
    g<-which(gf==k&eligible);bycase<-split(g,case[gi[g]])
    cp<-vapply(bycase,function(i){z<-gprofile[j,i];z<-z[is.finite(z)];if(length(z))median(z) else NA_real_},numeric(1));cp<-cp[is.finite(cp)]
    if(length(cp))profiles[j,k]<-median(cp)
  }
  list(fine=f,image=im,groups=groups,group_fine=gf,group_image=gi,eligible=eligible,
    image_patient=case,group_profiles=gprofile,profiles=profiles)
}

derive_LN_hierarchical_protein_states_RAM <- function(x,config,assay_name="asin",rank_fun=function(z)rank(z,ties.method="average"),area_adjust=FALSE,residual_fun=function(y,x){fit<-lm.fit(cbind(1,x),y,tol=1e-10);if(fit$rank<2L)rep(NA_real_,length(y)) else fit$residuals}) {
  av<-function(v,k)attr(v,k,exact=TRUE);dl<-function(v){if(is.list(v)&&!isS4(v))return(v);attr(v,"listData",exact=TRUE)}
  a<-dl(av(av(x,"assays"),"data"))[[assay_name]];cd<-dl(av(x,"colData"));stopifnot(identical(rownames(a),unlist(config$physical_panel)))
  if(area_adjust&&assay_name=="counts")a<-asinh(a)
  fine<-as.character(cd[[config$fine_field]]);image<-as.character(cd[[config$image_field]]);patient<-as.character(cd[[config$patient_field]])
  valid<-function(z)!is.na(z)&!tolower(trimws(z))%in%c("","na","unknown","none","null")
  take<-valid(fine)&valid(image)&valid(patient);n<-ncol(a)
  lm<-c("CD4","CD8a","CD14","CD68","CD163","CD19","CD56","SMA","CD45RA","CD45RO")
  p<-LN_rank_profiles_RAM(a,lm,fine,image,patient,take,config$minimum_cells_per_image_group,rank_fun)
  rule<-function(r,k) {
    if(any(!is.finite(r)))return(FALSE)
    other<-if(k=="CD4")"CD8a" else "CD4"
    opponent<-max(mean(r[c("CD14","CD68","CD163")]),r[c("CD19","CD56","SMA",other)])
    r[k]-opponent>=config$lineage_margin&&max(r[c("CD45RA","CD45RO")])>=config$minimum_CD45_relative_rank
  }
  coverage<-function(p,k){g<-which(p$group_fine==k&p$eligible);list(g=g,cases=length(unique(p$image_patient[p$group_image[g]])),cells=sum(take&fine==p$fine[k]))}
  parents<-list(CD4=character(),CD8a=character());parent_records<-list()
  for(k in seq_along(p$fine)) {
    cover<-coverage(p,k);eligible<-length(cover$g)>=config$minimum_images&&cover$cases>=config$minimum_patients&&cover$cells>=config$minimum_cells
    passed<-character();agreements<-list()
    for(label in names(parents)) {
      agreement<-if(length(cover$g))mean(vapply(cover$g,function(g)rule(p$group_profiles[,g],label),logical(1))) else 0
      ok<-eligible&&rule(p$profiles[,k],label)&&agreement>=config$minimum_image_agreement
      if(ok){parents[[label]]<-c(parents[[label]],p$fine[k]);passed<-c(passed,label)}
      agreements[[label]]<-agreement
    }
    parent_records[[k]]<-list(author_fine_group=p$fine[k],candidate_lineages=as.list(passed),eligible_images=length(cover$g),candidate_patient_labels=cover$cases,cell_rows=cover$cells,
      image_agreement=agreements,relative_marker_candidate_not_formal_identity=TRUE)
  }
  stopifnot(!length(intersect(parents$CD4,parents$CD8a)))
  parent_masks<-list(CD4=take&fine%in%parents$CD4,CD8a=take&fine%in%parents$CD8a)
  area<-if(area_adjust)as.numeric(cd[[config$area_field]]) else rep(1,n)
  if(area_adjust)parent_masks<-lapply(parent_masks,function(v)v&is.finite(area)&area>0)
  definitions<-list(
    R0_FOXP3_CTLA4=list(parent="CD4",markers=c("FoxP3","CTLA4"),cut=c(.75,.60)),
    R1_FOXP3_CTLA4_MKI67=list(parent="CD4",markers=c("FoxP3","CTLA4","MKI67"),cut=c(.75,.60,.75)),
    R2_general_CD4_MKI67=list(parent="CD4",markers="MKI67",cut=.75),
    E0_PD1_TIM3=list(parent="CD8a",markers=c("PD1","TIM3"),cut=c(.70,.70)),
    E1_CD8_MKI67=list(parent="CD8a",markers="MKI67",cut=.75),
    E2_CCR7_CD127=list(parent="CD8a",markers=c("CCR7","CD127"),cut=c(.70,.70)))
  membership<-setNames(rep(list(rep(FALSE,n)),length(definitions)),names(definitions));testable_images<-setNames(rep(list(character()),length(definitions)),names(definitions))
  for(lineage in names(parents)) {
    inside<-parent_masks[[lineage]];if(!any(inside))next
    defs<-definitions[vapply(definitions,function(d)d$parent==lineage,logical(1))]
    markers<-unique(unlist(lapply(defs,function(d)d$markers),use.names=FALSE));imgs<-split(which(inside),image[inside])
    rz<-setNames(rep(list(rep(NA_real_,n)),length(markers)),markers)
    varying<-matrix(FALSE,length(imgs),length(markers),dimnames=list(names(imgs),markers))
    for(marker in markers) {
      z<-as.numeric(a[marker,])
      for(im in names(imgs)) {
        i<-imgs[[im]];i<-i[is.finite(z[i])]
        if(length(i)<config$minimum_parent_cells_per_image||length(unique(z[i]))<2L)next
        y<-z[i]
        if(area_adjust){lx<-log(area[i]);if(var(lx)<=1e-12)next;y<-residual_fun(y,lx);if(any(!is.finite(y))||diff(range(y))<=1e-10*max(1,diff(range(z[i]))))next}
        rz[[marker]][i]<-rank_fun(y)/(length(i)+1);varying[im,marker]<-TRUE
      }
    }
    for(name in names(defs)) {
      def<-defs[[name]];good_images<-rownames(varying)[rowSums(varying[,def$markers,drop=FALSE])==length(def$markers)]
      mask<-inside&image%in%good_images
      for(j in seq_along(def$markers)){v<-rz[[def$markers[j]]];mask<-mask&is.finite(v)&v>=def$cut[j]}
      mask[is.na(mask)]<-FALSE;membership[[name]]<-mask;testable_images[[name]]<-good_images
    }
  }
  state_records<-lapply(names(definitions),function(name) {
    m<-membership[[name]];def<-definitions[[name]]
    list(state=name,parent=def$parent,fixed_markers=as.list(def$markers),fixed_relative_cutoffs=as.list(def$cut),
      cell_rows=sum(m),candidate_patient_labels=length(unique(patient[m])),image_groups_with_state=length(unique(image[m])),
      images_with_measurement_support=length(testable_images[[name]]),
      within_parent_relative_state_not_absolute_positivity=TRUE)
  })
  list(summary=list(status="FIXED_HIERARCHICAL_RELATIVE_MARKER_QUALIFICATION_ONLY",parent_records=parent_records,state_records=state_records,
      resolved_parent_groups=lapply(parents,as.list),state_splitting_independent_of_author_fine_group_resolution=TRUE,
      cell_area_adjusted_state_scoring=area_adjust,relative_markers_not_absolute_positivity_or_TCR_truth=TRUE,rank_defined_state_fractions_not_absolute_clinical_abundance=TRUE,
      all_individual_membership_and_rank_profiles_RAM_only=TRUE),
    internal=list(states=membership,parents=parent_masks,fine=fine,image=image,patient=patient,eligible=take,testable_images=testable_images))
}
