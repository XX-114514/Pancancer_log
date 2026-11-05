#!/bin/bash -l
#SBATCH --job-name=cm-bench
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%A_%a.out
# Usage: sbatch --array=0-N --export=ALL,CODE=<repo>/scripts/nmf_benchmark/cm,PC_ROOT=<root>,RUN=<run_dir>,TASKS=<tasks.tsv> cm_bench.sbatch
# TASKS: one "method<TAB>seed<TAB>K" line per array index.
set -eo pipefail  # no -u: conda activate.d scripts reference unset vars
read -r METHOD SEED K < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TASKS")
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
module load micromamba
export MAMBA_ROOT_PREFIX=${PC_ROOT}/envs/micromamba
eval "$(micromamba shell hook -s bash)"
case "$METHOD" in
  NMF_nsNMF)  # frozen CoVarNet env (NMF 0.28), as the reference run
    micromamba activate "$PC_ROOT/envs/conda/scrna_r_covarnet2"; export R_LIBS= R_LIBS_USER=
    Rscript "$CODE/run_r_method.R" "$RUN" "$METHOD" "$SEED" "$K" ;;
  RcppML_*|CoGAPS|LIGER_*|log1pNMF_*)
    micromamba activate "$PC_ROOT/envs/conda/scrna_r_full"
    Rscript "$CODE/run_r_method.R" "$RUN" "$METHOD" "$SEED" "$K" "$PC_ROOT/envs/rlib/nmf_r43" ;;
  *)
    micromamba activate "$PC_ROOT/envs/conda/nmf_py"
    python -I "$CODE/run_py_method.py" "$RUN" "$METHOD" "$SEED" "$K" ;;
esac
