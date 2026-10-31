#!/bin/bash -l
#SBATCH --job-name=env-Rfix
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --constraint=AVX512
#SBATCH --output=%x-%j.out
# Usage: sbatch --export=ALL,MIG_CODE_ROOT=<repo>/scripts/wehi_migration,PC_ROOT=<root>,NAME=<env> jobs/fix_R_full_env.sbatch
# Second stage for build_R_full_env.sbatch (2026-11-01 findings):
# - ragg 1.2.6 needs -fpermissive with the current freetype headers (unsigned char* vs char*);
# - lme4 1.1-34 does not compile against Matrix 1.6-5's CHOLMOD headers, so lme4 1.1-35.5 is used
#   and reported as a version deviation;
# then the exact installer is rerun for every package still missing or mismatched.
set -euo pipefail
echo "repo_commit $(git -C "${MIG_CODE_ROOT:?}" rev-parse HEAD) dirty $(git -C "$MIG_CODE_ROOT" status --porcelain | wc -l)"
module load micromamba
export MAMBA_ROOT_PREFIX=/vast/scratch/users/$USER/micromamba
SPEC=$PC_ROOT/envs/specs/scrna_r
P=$PC_ROOT/envs/conda/${NAME:?}
REP=$P/rebuild_report
set +u; eval "$(micromamba shell hook -s bash)"; micromamba activate "$P"; set -u
export CPATH="$P/include/freetype2:$P/include${CPATH:+:$CPATH}"
export UDUNITS2_INCLUDE=$P/include UDUNITS2_LIBS=$P/lib
MV=$(mktemp); echo 'CXXFLAGS += -fpermissive' > "$MV"
R_MAKEVARS_USER=$MV Rscript -e 'remotes::install_version("ragg", version = "1.2.6", dependencies = FALSE, upgrade = "never")'
Rscript -e 'remotes::install_version("lme4", version = "1.1-35.5", dependencies = FALSE, upgrade = "never")'
printf 'package\told\twehi\treason\nlme4\t1.1-34\t%s\t1.1-34 does not compile against Matrix 1.6-5 CHOLMOD headers\n' \
  "$(Rscript -e 'cat(as.character(packageVersion("lme4")))')" > "$REP/R_version_deviations.tsv"
PKGS=$(awk -F'\t' 'NR>1 && $4!="installed_exact" && $4!="already_exact" && $4!="not_in_old_env" && $1!="lme4" {print $1}' "$REP/R_exact_install.tsv" | paste -sd:)
echo "retry: $PKGS"
Rscript "$MIG_CODE_ROOT/wehi_install_R_exact.R" "$SPEC/R_packages.tsv" "$PKGS" "$REP/R_exact_install_fix.tsv"
Rscript "$MIG_CODE_ROOT/wehi_install_R_github.R" "$SPEC/R_packages.tsv" "$SPEC/R_nonstandard_library" "$REP/R_github_install_fix.tsv"
Rscript -e 'ip <- installed.packages()[, c("Package","Version")]; write.table(ip, "'"$REP"'/R_installed_wehi.tsv", sep="\t", quote=FALSE, row.names=FALSE)'
for p in Seurat CellChat copykat scATOMIC infercnv NMF SingleR monocle3 GSVA liana; do
  Rscript -e "suppressMessages(library($p)); cat('$p', as.character(packageVersion('$p')), 'LOAD_OK\n')" 2>&1 | tail -1 || true
done
