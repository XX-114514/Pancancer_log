# Metadata and coordinate qualification only. Inputs and intermediate keys stay
# in this R process; the result contains source aggregates and no unit records.
# This never selects fields from a native object or establishes clinical units.

unit_key_RAM <- function(columns) {
  stopifnot(length(columns) > 0L)
  stopifnot(length(unique(lengths(columns)))==1L)
  if (!length(columns[[1L]])) return(character())
  encoded <- lapply(columns, function(v) {
    z <- enc2utf8(as.character(v))
    paste0(nchar(z, type="bytes"), ":", z)
  })
  Reduce(paste0, encoded)
}

qualify_LN_unit_geometry_RAM <- function(coordinates, units) {
  stopifnot(is.matrix(coordinates), is.numeric(coordinates),
    ncol(coordinates)==2L, is.data.frame(units),
    identical(names(units), c("patient", "node", "ROI")),
    nrow(coordinates)==nrow(units), nrow(units)>0L,
    all(vapply(units, function(v) is.character(v)||is.factor(v), logical(1))))
  units <- lapply(units, as.character)
  missing <- lapply(units, function(v) is.na(v)|tolower(trimws(v)) %in%
    c("", "na", "nan", "none", "null", "unknown", "n/a"))
  complete <- !Reduce(`|`, missing)
  finite <- is.finite(coordinates[,1]) & is.finite(coordinates[,2])
  take <- which(complete & finite)
  # Length-prefixed literal tuples, not delimiter/prefix guesses or patient
  # inference. Repeated station/ROI labels can occur in distinct patients.
  eligible <- lapply(units, function(v) v[take])
  patient_key <- eligible$patient
  node_key <- unit_key_RAM(eligible[c("patient", "node")])
  roi_key <- unit_key_RAM(eligible)
  groups <- split(seq_along(take), roi_key)
  nearest <- numeric(length(take))
  nearest[] <- NA_real_
  duplicates <- singleton <- 0L
  for (g in groups) {
    xy <- coordinates[take[g],,drop=FALSE]
    duplicates <- duplicates + sum(duplicated(as.data.frame(xy)))
    if (length(g)==1L) {
      singleton <- singleton+1L
      next
    }
    nn <- RANN::nn2(xy, query=xy, k=2L, eps=0)
    for (i in seq_along(g)) {
      # With coincident centroids the query itself need not be first, or even
      # among the first two ties. Pick an actual nonself index explicitly.
      j <- which(nn$nn.idx[i,] != i)[1L]
      stopifnot(!is.na(j))
      nearest[g[i]] <- nn$nn.dists[i,j]
    }
  }
  nesting <- function(child, parent) {
    edges <- unique(data.frame(child=child, parent=parent))
    sizes <- table(edges$child)
    list(distinct_child_labels=length(sizes),
      child_labels_with_multiple_parents=sum(sizes>1L),
      maximum_parents_per_child=if(length(sizes))as.integer(max(sizes)) else NULL)
  }
  summary <- list(status="COMPLETE_UNIT_AND_NATIVE_COORDINATE_QUALIFICATION_NO_EFFECTS",
    input_cell_rows=nrow(coordinates), complete_unit_rows=sum(complete),
    finite_coordinate_rows=sum(finite), jointly_eligible_rows=length(take),
    missing_unit_rows=lapply(missing, sum),
    eligible_distinct_literal_patients=length(unique(patient_key)),
    eligible_distinct_patient_node_tuples=length(unique(node_key)),
    eligible_distinct_patient_node_ROI_tuples=length(groups),
    raw_node_label_to_patient=nesting(eligible$node, patient_key),
    raw_ROI_label_to_patient_node=nesting(eligible$ROI, node_key),
    patient_node_ROI_to_patient_node=nesting(roi_key, node_key),
    ROI_singletons=singleton, repeated_coordinate_rows_within_ROI=duplicates,
    nonself_nearest_distance_rows=sum(is.finite(nearest)),
    nonself_zero_distance_rows=sum(nearest==0, na.rm=TRUE),
    nonself_nearest_distance_quantiles_native_units=if(any(is.finite(nearest)))
      as.list(setNames(as.numeric(quantile(nearest[is.finite(nearest)],
        probs=c(0,.25,.5,.75,1), names=FALSE, type=7)),c("min","q25","median","q75","max"))) else NULL,
    all_neighbors_restricted_to_literal_patient_node_ROI=TRUE,
    duplicate_coordinates_not_deleted_or_counted_as_new_patients=TRUE,
    physical_coordinate_units_not_inferred=TRUE,
    literal_unit_counts_not_accepted_independent_patient_semantics=TRUE,
    no_state_thresholds_graph_effects_or_clinical_values=TRUE,
    no_identifiers_individual_distances_or_records_in_result=TRUE)
  # Keys and per-cell distances remain local; return the aggregate alone.
  summary
}

check_LN_unit_geometry_synthetic_RAM <- function() {
  a <- data.frame(patient=c("FAKE_A","FAKE_A","FAKE_B","FAKE_B","FAKE_B"),
    node=rep("shared-node",5), ROI=c("shared-ROI","shared-ROI","shared-ROI","shared-ROI","single"))
  xy <- matrix(c(0,0, 3,4, 0,0, 0,0, 9,9), ncol=2, byrow=TRUE)
  q <- qualify_LN_unit_geometry_RAM(xy,a)
  stopifnot(q$eligible_distinct_literal_patients==2L,
    q$eligible_distinct_patient_node_tuples==2L,
    q$eligible_distinct_patient_node_ROI_tuples==3L,
    q$raw_node_label_to_patient$child_labels_with_multiple_parents==1L,
    q$raw_ROI_label_to_patient_node$child_labels_with_multiple_parents==1L,
    q$patient_node_ROI_to_patient_node$child_labels_with_multiple_parents==0L,
    q$ROI_singletons==1L, q$repeated_coordinate_rows_within_ROI==1L,
    q$nonself_nearest_distance_rows==4L, q$nonself_zero_distance_rows==2L)
  # Independent full pairwise distances inside explicit unit tuples. This
  # checks nearest-distance values, not just the production implementation.
  expected <- c(5,5,0,0)
  oracle <- as.list(setNames(as.numeric(quantile(expected,c(0,.25,.5,.75,1),
    names=FALSE,type=7)),c("min","q25","median","q75","max")))
  stopifnot(identical(q$nonself_nearest_distance_quantiles_native_units,oracle))
  perm <- c(5,3,1,4,2)
  stopifnot(identical(q,qualify_LN_unit_geometry_RAM(xy[perm,,drop=FALSE],a[perm,,drop=FALSE])))
  bad <- a;bad$patient[1] <- NA_character_;bad$ROI[5] <- "unknown"
  zz <- xy;zz[2,1] <- Inf
  qb <- qualify_LN_unit_geometry_RAM(zz,bad)
  stopifnot(qb$jointly_eligible_rows==2L,qb$nonself_zero_distance_rows==2L)
  empty <- a;empty$patient[] <- NA_character_
  qe <- qualify_LN_unit_geometry_RAM(xy,empty)
  stopifnot(qe$jointly_eligible_rows==0L,is.null(qe$nonself_nearest_distance_quantiles_native_units))
  # Plain delimiter joining would merge these distinct tuples.
  k <- unit_key_RAM(list(c("a:b","a"),c("c","b:c")))
  stopifnot(length(unique(k))==2L)
  stopifnot(!grepl("FAKE_|shared-node|shared-ROI",jsonlite::toJSON(q,auto_unbox=TRUE)))
  # Larger deterministic geometry with ties, independently checked by dist().
  x <- cbind(rep(0:8,each=7),rep(0:6,9))
  x <- rbind(x,x[1:3,,drop=FALSE])
  u <- data.frame(patient=rep("SYNTHETIC",nrow(x)),node="N",ROI="R")
  d <- as.matrix(dist(x));diag(d) <- Inf
  truth <- apply(d,1,min)
  qd <- qualify_LN_unit_geometry_RAM(x,u)
  observed <- unlist(qd$nonself_nearest_distance_quantiles_native_units,use.names=FALSE)
  stopifnot(max(abs(observed-as.numeric(quantile(truth,c(0,.25,.5,.75,1),names=FALSE))))<1e-12)
  list(status="PASS_SYNTHETIC_NESTING_CROSS_ROI_TIES_AND_PAIRWISE_DISTANCE_ORACLE",
    no_real_source_object_or_individual_records_read=TRUE,
    original_live_worker_untouched=TRUE, actual_native_fields_not_selected=TRUE,
    native_measurement_or_spatial_effect_not_accepted=TRUE)
}
