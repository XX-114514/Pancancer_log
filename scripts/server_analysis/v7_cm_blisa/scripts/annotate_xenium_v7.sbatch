#!/bin/bash
#SBATCH -J v7_xen_annot
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 12:00:00
#SBATCH -o %x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=16
${PC_ROOT}/envs/conda/scanpy/bin/python -I $RUN/scripts/annotate_xenium_v7.py $RUN $(cat $RUN/config/xenium_tumour_sections.txt)
