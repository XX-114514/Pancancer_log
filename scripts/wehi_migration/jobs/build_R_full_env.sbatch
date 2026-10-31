#!/bin/bash -l
#SBATCH --job-name=env-Rfull
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,MIG_CODE_ROOT=<repo>/scripts/wehi_migration,PC_ROOT=<root>,NAME=<env>[,RESUME=1] jobs/build_R_full_env.sbatch
# Full old scrna_r R library on WEHI: conda-forge r-base 4.3.3 (same build as the old env) and the
# system libraries the old env carried, in one solve; then every CRAN/Bioconductor package at its
# old version (wehi_install_R_exact.R), then GitHub/other packages (wehi_install_R_github.R).
set -euo pipefail
echo "repo_commit $(git -C "${MIG_CODE_ROOT:?}" rev-parse HEAD) dirty $(git -C "$MIG_CODE_ROOT" status --porcelain | wc -l)"
module load micromamba
export MAMBA_ROOT_PREFIX=/vast/scratch/users/$USER/micromamba
SPEC=$PC_ROOT/envs/specs/scrna_r
P=$PC_ROOT/envs/conda/${NAME:?}
if [ "${RESUME:-0}" = 1 ]; then
  [ -d "$P" ] || { echo "missing $P" >&2; exit 2; }
else
  [ -e "$P" ] && { echo "exists $P" >&2; exit 2; }
  micromamba create -y -p "$P" -c conda-forge --override-channels "r-base=4.3.3=h5074ccb_16" \
    compilers make pkg-config cmake glpk gmp mpfr gsl libxml2 libxml2-devel zlib liblzma-devel \
    hdf5 jags gdal geos proj udunits2 fftw cairo harfbuzz fribidi freetype libpng libtiff libjpeg-turbo \
    libgit2 libcurl openssl imagemagick r-remotes r-biocmanager
fi
mkdir -p "$P/rebuild_report"
eval "$(micromamba shell hook -s bash)"; micromamba activate "$P"
PKGS=$(awk -F'\t' 'NR>1 && $11!="base" && $11!="github" && $11!="other" {print $1}' "$SPEC/R_packages.tsv" | paste -sd:)
Rscript "$MIG_CODE_ROOT/wehi_install_R_exact.R" "$SPEC/R_packages.tsv" "$PKGS" "$P/rebuild_report/R_exact_install.tsv"
Rscript "$MIG_CODE_ROOT/wehi_install_R_github.R" "$SPEC/R_packages.tsv" "$SPEC/R_nonstandard_library" "$P/rebuild_report/R_github_install.tsv"
Rscript -e 'ip <- installed.packages()[, c("Package","Version")]; write.table(ip, "'"$P"'/rebuild_report/R_installed_wehi.tsv", sep="\t", quote=FALSE, row.names=FALSE)'
for p in Seurat CellChat copykat scATOMIC infercnv NMF SingleR monocle3 GSVA; do
  Rscript -e "suppressMessages(library($p)); cat('$p', as.character(packageVersion('$p')), 'LOAD_OK\n')" 2>&1 | tail -1
done
