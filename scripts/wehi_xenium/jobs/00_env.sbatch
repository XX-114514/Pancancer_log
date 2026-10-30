#!/bin/bash -l
#SBATCH --job-name=pc-spatial-env
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
# Submit from a run/log directory with XENIUM_CODE_ROOT exported (sbatch passes the environment).
set -euo pipefail
source "${XENIUM_CODE_ROOT:?export XENIUM_CODE_ROOT}/config.sh"
module load micromamba/latest
YML="${XENIUM_CODE_ROOT}/envs/pc-spatial.yml"
if [[ -x "$ENV_PREFIX/bin/python" ]]; then
  micromamba install -y --prefix "$ENV_PREFIX" -f "$YML"
else
  micromamba create -y --prefix "$ENV_PREFIX" -f "$YML"
fi
"$ENV_PREFIX/bin/python" -c 'import scanpy, anndata, squidpy, spatialdata, spatialdata_io, pyarrow, h5py; print("Imports OK")'
# Lock exports go to the data root; copy reviewed versions into envs/ when the env changes.
mkdir -p "${META_ROOT}/env_lock"
micromamba list -p "$ENV_PREFIX" --explicit > "${META_ROOT}/env_lock/pc-spatial.explicit.txt"
"$ENV_PREFIX/bin/python" -m pip list --format=freeze > "${META_ROOT}/env_lock/pc-spatial.versions.txt"
