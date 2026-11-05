#!/bin/bash -l
#SBATCH --job-name=env-nmf
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,CODE=<repo>/scripts/nmf_benchmark,PC_ROOT=<root>[,PART=py|r|all][,REINSTALL=pkgA:pkgB] build_nmf_envs.sbatch
# py: conda env nmf_py (cNMF from source + latest scikit-learn/scanpy)
# r : overlay lib envs/rlib/nmf_r43 on scrna_r_full (RcppML 1.0, CoGAPS, log1pNMF, rliger from GitHub HEAD)
set -eo pipefail  # no -u: conda activate.d scripts reference unset vars
R0=$PC_ROOT/resources/nmf_benchmark; SRC=$R0/src
module load micromamba
export MAMBA_ROOT_PREFIX=${PC_ROOT}/envs/micromamba
eval "$(micromamba shell hook -s bash)"
PART=${PART:-all}
if [ "$PART" != r ]; then
  P=$PC_ROOT/envs/conda/nmf_py
  [ -e "$P" ] || micromamba create -y -p "$P" -c conda-forge --override-channels \
      python=3.11 "scikit-learn>=1.5" scanpy anndata pandas numpy scipy matplotlib palettable pyyaml \
      pytest ipykernel squidpy h5py
  micromamba activate "$P"
  pip install --no-deps --no-build-isolation "$SRC/cNMF"
  python -I -c "import cnmf, sklearn, scanpy; print('cnmf', cnmf.__version__ if hasattr(cnmf,'__version__') else '?', 'sklearn', sklearn.__version__, 'scanpy', scanpy.__version__)"
  micromamba env export -p "$P" > "$R0/nmf_py_env_export.yml"
  micromamba deactivate
fi
if [ "$PART" != py ]; then
  # RcppPlanc needs hwloc headers; the conda compiler's sysroot ignores /usr/include, so take
  # hwloc from a small side prefix rather than re-solving scrna_r_full (that swaps r-base/ICU).
  H=$PC_ROOT/envs/conda/nmf_sysdeps
  [ -e "$H" ] || micromamba create -y -p "$H" -c conda-forge --override-channels "libhwloc>=2,<3" "openblas=0.3.34"
  micromamba activate "$PC_ROOT/envs/conda/scrna_r_full"
  export CPATH=$H/include LIBRARY_PATH=$H/lib PKG_CONFIG_PATH=$H/lib/pkgconfig CMAKE_PREFIX_PATH=$H
  export LDFLAGS="${LDFLAGS} -L$H/lib -Wl,-rpath,$H/lib"
  # RcppPlanc's cmake otherwise picks the system FlexiBLAS/hwloc in /lib64 -> segfault on load
  E=$PC_ROOT/envs/conda/scrna_r_full
  # (cblas.h comes from openblas in nmf_sysdeps; the library linked is R's own libblas -> openblasp 0.3.34)
  export CMAKE_ARGS="-DCMAKE_DISABLE_FIND_PACKAGE_MKL=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_OpenBLAS=TRUE -DBLAS_LIBRARIES=$E/lib/libblas.so -DLAPACK_LIBRARIES=$E/lib/liblapack.so -DCBLAS_H_DIR=$H/include -DCMAKE_PREFIX_PATH=$H;$E -DCMAKE_IGNORE_PREFIX_PATH=/usr;/"
  # R links the final .so with Makeconf LDFLAGS, not the environment's: add the hwloc rpath there
  export R_MAKEVARS_USER=$PC_ROOT/resources/nmf_benchmark/Makevars.nmf
  printf 'LDFLAGS += -L%s/lib -Wl,-rpath,%s/lib\n' "$H" "$H" > "$R_MAKEVARS_USER"
  Rscript "$CODE/install_nmf_r.R" "$PC_ROOT/envs/rlib/nmf_r43" "$SRC"
fi
for n in cNMF RcppML CoGAPS log1pNMF liger; do
  printf '%s\t%s\t%s\n' "$n" "$(git -C "$SRC/$n" rev-parse HEAD)" "$(git -C "$SRC/$n" remote get-url origin)"
done > "$R0/source_commits.tsv"
echo done
