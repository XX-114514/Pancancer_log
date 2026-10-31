#!/bin/bash -l
#SBATCH --job-name=env-Rsubset
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,MIG_CODE_ROOT=<repo>/scripts/wehi_migration,PC_ROOT=<root>,NAME=<env>,PKGS=a:b:c[,RESUME=1] jobs/build_R_subset_env.sbatch
# PKGS is colon-separated: sbatch --export splits its argument on commas.
# RESUME=1 reuses an existing env and only (re)installs packages that are not yet exact.
# conda-forge r-base 4.3.3 (same build as the old scrna_r env) + toolchain, then the exact
# old-server versions of PKGS and their dependency closure from CRAN archive/Bioconductor 3.18.
set -euo pipefail
echo "repo_commit $(git -C "${MIG_CODE_ROOT:?}" rev-parse HEAD) dirty $(git -C "$MIG_CODE_ROOT" status --porcelain | wc -l)"
module load micromamba
export MAMBA_ROOT_PREFIX=/vast/scratch/users/$USER/micromamba
P=$PC_ROOT/envs/conda/${NAME:?}
if [ "${RESUME:-0}" = 1 ]; then
  [ -d "$P" ] || { echo "missing $P" >&2; exit 2; }
else
  [ -e "$P" ] && { echo "exists $P" >&2; exit 2; }
  # Everything in one solve with the r-base build pinned: adding packages later let the solver
  # swap the r-base build and ICU underneath already compiled packages.
  micromamba create -y -p "$P" -c conda-forge --override-channels "r-base=4.3.3=h5074ccb_16" compilers make pkg-config glpk gmp libxml2 libxml2-devel zlib liblzma-devel r-remotes r-biocmanager
fi
mkdir -p "$P/rebuild_report"
# Activate (not `micromamba run`): compiler activation scripts must set CC/CXX, and configure
# scripts need an open fd 5, which `micromamba run` does not provide.
eval "$(micromamba shell hook -s bash)"; micromamba activate "$P"
Rscript "$MIG_CODE_ROOT/wehi_install_R_exact.R" "$PC_ROOT/envs/specs/scrna_r/R_packages.tsv" "${PKGS:?}" "$P/rebuild_report/R_exact_install.tsv"
