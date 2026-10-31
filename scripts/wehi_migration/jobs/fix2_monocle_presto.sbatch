#!/bin/bash -l
#SBATCH --job-name=env-Rfix2
#SBATCH --partition=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%j.out
set -eo pipefail
module load micromamba; export MAMBA_ROOT_PREFIX=/vast/scratch/users/$USER/micromamba
P=/vast/scratch/users/xiong.x/pancancer/envs/conda/scrna_r_full
eval "$(micromamba shell hook -s bash)"; micromamba activate "$P"
export CPATH="$P/include/freetype2:$P/include"
Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); if (!requireNamespace("terra", quietly=TRUE)) install.packages("terra")
for (x in list(c("monocle3","cole-trapnell-lab/monocle3","1.3.1"), c("monocle3","cole-trapnell-lab/monocle3","1.3.1"))) {
  ok <- FALSE
  for (ref in c(x[3], paste0("v", x[3]))) {
    r <- try(remotes::install_github(x[2], ref = ref, dependencies = FALSE, upgrade = "never", force = TRUE), silent = TRUE)
    if (!inherits(r, "try-error") && requireNamespace(x[1], quietly = TRUE)) { ok <- TRUE; break }
  }
  cat(x[1], if (ok) as.character(packageVersion(x[1])) else "FAILED", "\n")
}'
for p in monocle3 presto; do Rscript -e "suppressMessages(library($p)); cat('$p', as.character(packageVersion('$p')), 'LOAD_OK\n')" 2>&1 | tail -1 || true; done
