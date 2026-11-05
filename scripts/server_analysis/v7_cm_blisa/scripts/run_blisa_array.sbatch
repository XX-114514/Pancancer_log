#!/bin/bash
#SBATCH -J v7_blisa
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH -o logs/blisa/%x-%A_%a.out
set -euo pipefail
module load ImageMagick/7.1.1
export PATH=${PC_ROOT}/envs/conda/scrna_r_full/bin:$PATH R_LIBS=${PC_ROOT}/envs/rlib/blisa_r43
mapfile -t SECS < <(sed 's|.*/||; s|\.counts\.h5ad$||' $RUN/config/xenium_tumour_sections.txt)
BINS=(50 25 100)
i=$SLURM_ARRAY_TASK_ID; sec=${SECS[$((i % ${#SECS[@]}))]}; bin=${BINS[$((i / ${#SECS[@]}))]}
[ -f $RUN/xenium/$sec/blisa_bin$bin/summary.tsv ] && { echo "skip $sec $bin"; exit 0; }
Rscript $RUN/scripts/run_blisa_section.R $RUN/xenium/$sec $bin 8
