# Usage: bash scripts/with_r.sh scripts/run_one_section.R input.rds output_dir
# Input: one independently annotated section, raw counts, micron coordinates,
# cell_type + section_id columns. Caller must audit mapping to frozen V5 states.
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
source("scripts/check_environment.R")
spe <- readRDS(args[1])
stopifnot(inherits(spe,"SpatialExperiment"))
cd <- SummarizedExperiment::colData(spe)
stopifnot(all(c("cell_type","section_id") %in% names(cd)),
          !anyNA(cd$cell_type), !anyNA(cd$section_id),
          length(unique(cd$section_id))==1L)
xy <- SpatialExperiment::spatialCoords(spe)
stopifnot(ncol(xy)==2L,all(is.finite(xy)),!is.null(rownames(spe)),
          !anyDuplicated(rownames(spe)),!anyDuplicated(colnames(spe)))
x <- SummarizedExperiment::assay(spe,"counts")
v <- if(inherits(x,"sparseMatrix")) x@x else as.vector(x)
stopifnot(all(is.finite(v)),all(v>=0),all(abs(v-round(v))<1e-8))
if (dir.exists(args[2])) stop("Output directory exists; use a new run directory")
dir.create(args[2],recursive=TRUE)
a <- read.delim("tables/focus_axes_definition.tsv",check.names=FALSE)
# Complex receptor handling is deliberately gated pending package semantics audit.
a$status <- vapply(seq_len(nrow(a)),function(i) {
  if (grepl("_",a$ligand[i]) || grepl("_",a$receptor[i])) return("complex_pending_audit")
  if (!all(c(a$ligand[i],a$receptor[i]) %in% rownames(spe))) return("gene_not_measured")
  if (!all(c(a$source[i],a$target[i]) %in% cd$cell_type)) return("state_not_identified")
  "eligible"
},character(1))
write.table(a,file.path(args[2],"axis_eligibility.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
ok <- a[a$status=="eligible",]
if (!nrow(ok)) stop("No eligible frozen axes; do not substitute broader identities")
lr <- unique(data.frame(ligand.symbol=ok$ligand,receptor.symbol=ok$receptor))
# Nearby-mode sensitivity is explicit; this pilot is not a final mechanism-specific model.
set.seed(1337)
res <- blisa::blisa(spe,LR_df=lr,group="cell_type",bin_size=50,dmax=150,
                    nsim=999,cpu_threads=2,default_mode="nearby",platform="xenium")
saveRDS(res,file.path(args[2],"blisa_raw.rds"))
# Preserve all raw tests; apply BH over the tested LR x bin family per section.
p <- res$LR_results$all_pval
if (is.null(p)) stop("Expected all_pval absent; inspect installed API")
flat <- unlist(p,use.names=FALSE)
res$LR_results$all_q_BH_section <- split(p.adjust(flat,"BH"),rep(seq_along(p),lengths(p)))
saveRDS(res,file.path(args[2],"blisa_with_section_BH.rds"))
writeLines(capture.output(sessionInfo()),file.path(args[2],"sessionInfo.txt"))
writeLines("Exploratory section run: require HH direction and BH q<=0.05; CCI_scores are descriptive raw scores, not BH-filtered or causal.",file.path(args[2],"STATUS.txt"))
