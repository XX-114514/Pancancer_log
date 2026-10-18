# Read fixed official outputs, reconstruct every supported score, then project.
# No fitting of genes/PCA/DIALOGUE on heldout material and no clinical fields.
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]; models <- args[2]; out <- args[3]
stopifnot(file.exists(file.path(models,'private_terminal.json')), !dir.exists(out))
dir.create(out,recursive=TRUE,mode='0700')
suppressPackageStartupMessages(library(Matrix))
views <- c('R_FOXP3_capture','E_effector_capture','Myeloid')
official <- setNames(gsub('_','',views,fixed=TRUE),views)
fits <- data.frame(k=c(2,2,2,3,2),seed=c(20261018,20261019,20261020,20261018,20261018),
 mode=c(rep('source_adjusted',4),'source_permissive'))
checks <- list(); scope <- list(); weights <- list(); cases <- list(); errors <- list()
cell.file <- gzfile(file.path(out,'private_all_cell_score_reconstruction.tsv.gz'),'wt')
first <- TRUE
read.part <- function(view,part) {
 p <- file.path(input,'private_inputs',paste0(view,'_',part))
 md <- read.delim(paste0(p,'_metadata.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
 genes <- readLines(paste0(p,'_genes.tsv'))
 m <- readMM(gzfile(paste0(p,'_tpm.mtx.gz')))
 stopifnot(nrow(m)==length(genes),ncol(m)==nrow(md),!anyDuplicated(md$cell),
           !anyDuplicated(genes))
 rownames(m) <- genes
 list(md=md,m=m)
}
for(i in seq_len(nrow(fits))) {
 fit <- fits[i,]; id <- paste0('K',fit$k,'_seed',fit$seed,'_',fit$mode)
 path <- file.path(models,id)
 result.path <- list.files(path,pattern='^DLG[.]output_.*[.]rds$',full.names=TRUE)
 if(length(result.path)!=1L) {
  for(view in views) for(program in paste0('MCP',seq_len(fit$k))) {
   key <- paste(id,view,program,sep='|')
   scope[[key]] <- data.frame(fit_id=id,view=view,official_view=official[[view]],
       program=program,status='OFFICIAL_FINAL_OBJECT_ABSENT',active_genes=0)
  }
  next
 }
 r <- readRDS(result.path)
 stopifnot(setequal(r$cell.types,unname(official)),
           setequal(names(r$scores),unname(official)),
           setequal(names(r$gene.pval),unname(official)))
 controls <- c('cellQ','geneQ',if(fit$mode=='source_adjusted') 'sourceCC','fraction','logN')
 stopifnot(setequal(as.character(r$param$conf),controls),
           setequal(as.character(r$param$covar),controls),is.null(r$param$pheno))
 for(view in views) {
  tr <- read.part(view,'discovery'); te <- read.part(view,'TNBC_holdout')
  stopifnot(identical(rownames(tr$m),rownames(te$m)),
      length(unique(tr$md$sample))==38L,length(unique(te$md$sample))==28L,
      all(tr$md$partition=='reviewed_case_discovery'),
      all(te$md$partition=='TNBC_source_holdout_sample_identity_pending'),
      all(te$md$dataset=='GSE305078'),all(is.na(te$md$sourceCC)),
      !any(te$md$cell %in% tr$md$cell),!any(te$md$sample %in% tr$md$sample))
  saved <- r$scores[[official[[view]]]]
  ix <- match(tr$md$cell,saved$cells)
  stopifnot(nrow(saved)==nrow(tr$md),!anyNA(ix),!anyDuplicated(saved$cells),
            identical(as.character(saved$samples[ix]),as.character(tr$md$sample)))
  gp <- r$gene.pval[[official[[view]]]]
  for(program in paste0('MCP',seq_len(fit$k))) {
   key <- paste(id,view,program,sep='|')
   active <- if(is.null(gp)) NULL else gp[gp$program==program &
       is.finite(gp$coef) & gp$coef>0,,drop=FALSE]
   n <- if(is.null(active)) 0L else nrow(active)
   scope[[key]] <- data.frame(fit_id=id,view=view,official_view=official[[view]],
      program=program,status='PENDING',active_genes=n)
   if(n==0L) {
    scope[[key]]$status <- 'NO_ACTIVE_NNLS_GENE_PROJECTION_UNSUPPORTED'
    next
   }
   one <- tryCatch({
    stopifnot(program %in% colnames(saved),!anyDuplicated(active$genes),
       all(active$genes %in% rownames(tr$m)),all(active$genes %in% rownames(te$m)))
    a <- as.matrix(tr$m[active$genes,,drop=FALSE])
    b <- as.matrix(te$m[active$genes,,drop=FALSE])
    mu <- rowMeans(a); sd <- apply(a,1,stats::sd)
    stopifnot(all(is.finite(mu)),all(is.finite(sd)),all(sd>0))
    w <- active$coef * ifelse(active$up,1,-1)
    za <- sweep(sweep(a,1,mu,'-'),1,sd,'/')
    zb <- sweep(sweep(b,1,mu,'-'),1,sd,'/')
    raw <- as.numeric(crossprod(w,za)); raw.h <- as.numeric(crossprod(w,zb))
    D <- cbind('(Intercept)'=1,as.matrix(tr$md[,controls,drop=FALSE]))
    stopifnot(all(is.finite(D)),qr(D)$rank==ncol(D))
    nuisance <- lm.fit(D,raw); score <- as.numeric(nuisance$residuals)
    err <- max(abs(score-as.numeric(saved[ix,program])))
    stopifnot(is.finite(err),err<1e-6)
    # A new source has no fitted source coefficient. Leave metadata NA.
    # Within-heldout-source associations are invariant to this unknown offset;
    # these scores cannot compare absolute levels across sources.
    common.controls <- setdiff(controls,'sourceCC')
    Dh <- cbind('(Intercept)'=1,as.matrix(te$md[,common.controls,drop=FALSE]))
    coef <- setNames(nuisance$coefficients,colnames(D))
    held <- raw.h-as.numeric(Dh %*% coef[colnames(Dh)])
    stopifnot(all(is.finite(score)),all(is.finite(held)))
    local.cases <- list(); local.cells <- list()
    for(part in c('discovery','TNBC_holdout')) {
     md <- if(part=='discovery') tr$md else te$md
     sc <- if(part=='discovery') score else held
     original <- if(part=='discovery') as.numeric(saved[ix,program]) else rep(NA_real_,length(sc))
     row <- data.frame(fit_id=id,view=view,program=program,partition=part,
       cell=md$cell,sample=md$sample,dataset=md$dataset,score=sc,
       official_saved_score=original,source_offset_unknown=part=='TNBC_holdout' && fit$mode=='source_adjusted')
     local.cells[[part]] <- row
     agg <- aggregate(sc,list(sample=md$sample,dataset=md$dataset),mean)
     names(agg)[3] <- 'score'; agg$fit_id <- id; agg$view <- view
     agg$program <- program; agg$partition <- part
     agg$source_offset_unknown <- part=='TNBC_holdout' && fit$mode=='source_adjusted'
     local.cases[[part]] <- agg
    }
    list(check=data.frame(fit_id=id,view=view,program=program,active_genes=n,
       discovery_cells=length(score),heldout_cells=length(held),max_saved_score_error=err,
       discovery_design_rank=qr(D)$rank,training_only_transform=TRUE,
       unknown_source_offset_not_estimated=fit$mode=='source_adjusted'),
       cells=do.call(rbind,local.cells),cases=do.call(rbind,local.cases),
       weights=data.frame(fit_id=id,view=view,program=program,gene=active$genes,
           weight=w,discovery_mean=mu,discovery_sd=sd))
   },error=function(e) {
    errors[[key]] <<- data.frame(fit_id=id,view=view,program=program,message=conditionMessage(e))
    NULL
   })
   if(is.null(one)) {
    scope[[key]]$status <- 'TECHNICAL_RECONSTRUCTION_FAILURE_REVIEW_REQUIRED'
    next
   }
   checks[[key]] <- one$check; weights[[key]] <- one$weights; cases[[key]] <- one$cases
   write.table(one$cells,cell.file,sep='\t',quote=FALSE,row.names=FALSE,
       col.names=first); first <- FALSE
   scope[[key]]$status <- 'OFFICIAL_PARITY_AND_FROZEN_SOURCE_PROJECTION_PASS'
  }
 }
}
close(cell.file)
write.rows <- function(rows,name) {
 if(length(rows)) write.table(do.call(rbind,rows),file.path(out,name),sep='\t',
     quote=FALSE,row.names=FALSE)
}
write.rows(scope,'private_all_planned_projection_slots.tsv')
write.rows(checks,'private_all_score_reconstruction_checks.tsv')
write.rows(errors,'private_all_projection_errors.tsv')
write.rows(weights,'private_frozen_gene_weights.tsv')
write.rows(cases,'private_all_case_mean_scores.tsv')
stopifnot(length(scope)==33L)
writeLines('All planned slots retained. Independent numeric audit, seed stability and source association inference are pending.',
 file.path(out,'private_interpretation_boundary.txt'))
cat('FIXED_OUTPUT_READER_FINISHED_INDEPENDENT_ACCEPTANCE_PENDING\n')
