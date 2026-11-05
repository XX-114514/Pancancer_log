#!/bin/bash -l
#SBATCH --job-name=nmf-smoke
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>/scripts/nmf_benchmark,PC_ROOT=<root> smoke_nmf.sbatch
set -eo pipefail
module load micromamba
export MAMBA_ROOT_PREFIX=${PC_ROOT}/envs/micromamba
eval "$(micromamba shell hook -s bash)"
O=${NMF_SMOKE_OUT:?Set NMF_SMOKE_OUT to a new directory}; mkdir "$O"
rc=0
micromamba activate "$PC_ROOT/envs/conda/nmf_py"
python -I "$CODE/smoke_py.py" "$O" || rc=1
micromamba deactivate
micromamba activate "$PC_ROOT/envs/conda/scrna_r_full"
Rscript "$CODE/smoke_r.R" "$PC_ROOT/envs/rlib/nmf_r43" || rc=1
echo "smoke rc=$rc"; exit $rc
