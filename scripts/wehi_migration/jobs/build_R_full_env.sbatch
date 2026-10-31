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
# system libraries the old env carried, in one solve. Pins (2026-10-31 dry runs): libxml2 2.13 as in the
# old env (libxml2-devel exists only for >=2.14 and then clashes with libgdal-core); hdf5 1.14 and JAGS 4.3
# for hdf5r 1.3.12 / rjags 4-17; gcc 11 as in the old env (gcc 15 breaks readxl and other older C++);
# udunits2 and X11 headers for units/sf and Cairo; imagemagick left out (needs ICU >=78, the r-base build needs 75.1); then every CRAN/Bioconductor package at its
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
    "gcc_linux-64=11" "gxx_linux-64=11" "gfortran_linux-64=11" make pkg-config cmake glpk gmp mpfr gsl \
    "libxml2=2.13" zlib liblzma-devel "hdf5=1.14" "jags=4.3" libgdal-core geos proj udunits2 \
    cairo harfbuzz pango fribidi freetype libpng libtiff libjpeg-turbo libgit2 libcurl openssl \
    xorg-libx11 xorg-xorgproto xorg-libxt xorg-libxrender xorg-libxext r-remotes r-biocmanager
fi
mkdir -p "$P/rebuild_report"
eval "$(micromamba shell hook -s bash)"; micromamba activate "$P"
# systemfonts includes <ft2build.h> without pkg-config flags; units looks for udunits2.h
export CPATH="$P/include/freetype2:$P/include${CPATH:+:$CPATH}"
export UDUNITS2_INCLUDE=$P/include UDUNITS2_LIBS=$P/lib
PKGS=$(awk -F'\t' 'NR>1 && $11!="base" && $11!="github" && $11!="other" {print $1}' "$SPEC/R_packages.tsv" | paste -sd:)
Rscript "$MIG_CODE_ROOT/wehi_install_R_exact.R" "$SPEC/R_packages.tsv" "$PKGS" "$P/rebuild_report/R_exact_install.tsv"
Rscript "$MIG_CODE_ROOT/wehi_install_R_github.R" "$SPEC/R_packages.tsv" "$SPEC/R_nonstandard_library" "$P/rebuild_report/R_github_install.tsv"
Rscript -e 'ip <- installed.packages()[, c("Package","Version")]; write.table(ip, "'"$P"'/rebuild_report/R_installed_wehi.tsv", sep="\t", quote=FALSE, row.names=FALSE)'
for p in Seurat CellChat copykat scATOMIC infercnv NMF SingleR monocle3 GSVA; do
  Rscript -e "suppressMessages(library($p)); cat('$p', as.character(packageVersion('$p')), 'LOAD_OK\n')" 2>&1 | tail -1
done
