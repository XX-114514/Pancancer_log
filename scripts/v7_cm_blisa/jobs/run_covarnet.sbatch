#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; set +a
ENV=${PC_ROOT}/envs/conda/scrna_r_covarnet2
export PATH=$ENV/bin:$PATH OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS= R_LIBS_USER=
Rscript $RUN/scripts/run_TME_CoVarNet_scheme_v1.R $RUN/config/covarnet_${BRANCH}.json
