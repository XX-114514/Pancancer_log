#!/bin/bash
#SBATCH -J v7_cm_blisa_agg
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 06:00:00
#SBATCH -o logs/%x-%j.out
set -euo pipefail
${PC_ROOT}/envs/conda/scanpy/bin/python -I $RUN/scripts/aggregate_cm_blisa.py $RUN $BRANCH $SEED $BIN $MEMB
