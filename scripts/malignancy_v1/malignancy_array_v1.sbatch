#!/bin/bash -l
#SBATCH --job-name=malig
#SBATCH --partition=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=2-00:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%A_%a.out
# Usage: sbatch --array=0-N --export=ALL,CODE=<repo>,PC_ROOT=<root>,LIST=<shard list>,OUT=<dir> malignancy_array_v1.sbatch
set -euo pipefail
SHARD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${LIST:?}"); SRC=$(basename "$SHARD" .h5ad)
[ -s "$OUT/$SRC.tsv.gz" ] && { echo "done already $SRC"; exit 0; }
echo "repo_commit $(git -C "${CODE:?}" rev-parse HEAD) source $SRC start $(date -u +%FT%TZ)"
W=$OUT/work/$SRC; mkdir -p "$W"
# scATOMIC uses MAGIC imputation (its default) via Rmagic/reticulate: point it at a plain venv with magic-impute
# (a conda-prefix Python makes reticulate look for a conda binary, which WEHI does not have)
export RETICULATE_PYTHON=$PC_ROOT/envs/py_magic_venv/bin/python
"$PC_ROOT/envs/conda/scanpy/bin/python" -I -u -W ignore "$CODE/scripts/malignancy_v1/malignancy_source_v1.py" "$SHARD" "$W" "$OUT/$SRC.tsv.gz.part" \
  "$PC_ROOT/envs/conda/scrna_r_full/bin/Rscript" "$CODE/scripts/malignancy_v1/malignancy_unit_worker.R" "$PC_ROOT/resources/infercnv/hg38_gencode_v27.txt" 8 2026110201
mv "$OUT/$SRC.tsv.gz.part" "$OUT/$SRC.tsv.gz"; mv "$OUT/$SRC.tsv.gz.part.units.jsonl" "$OUT/$SRC.units.jsonl"
echo "end $(date -u +%FT%TZ)"
