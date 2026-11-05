#!/bin/bash
#SBATCH -J cm_loso
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o logs/loso/%x-%A_%a.out
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; set +a
ENV=${PC_ROOT}/envs/conda/scrna_r_covarnet2
export PATH=$ENV/bin:$PATH OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS= R_LIBS_USER=
s=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $RUN/config/loso/sources.txt)
[ -f $RUN/outputs/loso/$s/CoVarNet_execution_report_v1.json ] && exit 0
Rscript $RUN/scripts/run_TME_CoVarNet_scheme_v1.R $RUN/config/loso/$s.json
