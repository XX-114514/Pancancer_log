# Run from the task directory via scripts/with_r.sh.
# Existing project libraries are read-only dependencies; installs go to env/R-library.
lib <- normalizePath("env/R-library", mustWork=TRUE)
options(timeout=180, Ncpus=2, repos=c(
  ChenLaboratory="https://chenlaboratory.r-universe.dev",
  CRAN="https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", lib=lib)
for (p in c("SpatialExperiment", "SummarizedExperiment", "ComplexHeatmap")) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, lib=lib, ask=FALSE, update=FALSE)
}
for (p in c("sf", "spdep", "fastLISA", "Matrix", "ggplot2", "viridisLite")) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p, lib=lib)
}
dest <- "env/blisa_1.0.3.tar.gz"
download.file("https://chenlaboratory.r-universe.dev/src/contrib/blisa_1.0.3.tar.gz", dest, mode="wb")
install.packages(dest, repos=NULL, type="source", lib=lib)
source("scripts/check_environment.R")
write.table(installed.packages()[,c("Package","Version","LibPath")],
            "logs/installed_packages.tsv", sep="\t", row.names=FALSE, quote=FALSE)
