source("scripts/check_environment.R")
library(SpatialExperiment)
set.seed(1337)
xy <- as.matrix(expand.grid(x=seq(0,990,10), y=seq(0,990,10)))
hot <- xy[,1] > 350 & xy[,1] < 650 & xy[,2] > 350 & xy[,2] < 650
counts <- rbind(LIG=rpois(nrow(xy), ifelse(hot,12,1)),
                REC=rpois(nrow(xy), ifelse(hot,12,1)),
                BG=rpois(nrow(xy), 4))
colnames(counts) <- paste0("synthetic_",seq_len(nrow(xy)))
spe <- SpatialExperiment(assays=list(counts=Matrix::Matrix(counts,sparse=TRUE)),
                        spatialCoords=xy,
                        colData=S4Vectors::DataFrame(cell_type=rep(c("A","B"),length.out=nrow(xy))))
lr <- data.frame(ligand.symbol="LIG",receptor.symbol="REC",annotation="Secreted Signaling")
res <- blisa::blisa(spe, LR_df=lr, group="cell_type", bin_size=50, dmax=150,
                    nsim=99, cpu_threads=2, min_ligand=1, min_receptor=1,
                    platform="xenium")
stopifnot(inherits(res,"blisa"),nrow(res$LR_results)>0,!is.null(res$CCI_scores))
saveRDS(res,"outputs/synthetic_smoke.rds")
writeLines("PASS: synthetic API/CCI execution only; not real-data validation", "logs/smoke_status.txt")
