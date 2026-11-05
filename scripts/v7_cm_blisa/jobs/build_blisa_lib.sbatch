#!/bin/bash
#SBATCH -J blisa-lib
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; set +a
module load ImageMagick/7.1.1
ENV=${PC_ROOT}/envs/conda/scrna_r_full
LIB=${PC_ROOT}/envs/rlib/blisa_r43
export PATH=$ENV/bin:$PATH PKG_CONFIG_PATH=$(dirname $(dirname $(command -v magick)))/lib/pkgconfig R_LIBS_USER=$LIB R_LIBS=$LIB MAKEFLAGS=-j8
BLISA_LIB=$LIB Rscript -e '
lib <- Sys.getenv("BLISA_LIB"); .libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
install.packages("magick", lib=lib)
BiocManager::install(c("SpatialExperiment"), lib=lib, update=FALSE, ask=FALSE, version="3.18")
install.packages("fastLISA", lib=lib)
install.packages("blisa", lib=lib, repos=c("https://chenlaboratory.r-universe.dev","https://cloud.r-project.org"))
suppressPackageStartupMessages(library(blisa)); cat("blisa", as.character(packageVersion("blisa")), "fastLISA", as.character(packageVersion("fastLISA")), "SpatialExperiment", as.character(packageVersion("SpatialExperiment")), "\n")
print(sessionInfo())'
