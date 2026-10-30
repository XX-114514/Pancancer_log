#!/usr/bin/env bash
# Old server: record conda environments exactly enough to rebuild them elsewhere.
# Usage: old_export_env.sh CONDA_BIN OUT_DIR ENV_PREFIX [ENV_PREFIX ...]
# Per environment: conda explicit list and export, pip freeze, and for R the installed package
# table with install origin (Repository/RemoteType/RemoteSha). R packages not from CRAN or
# Bioconductor (e.g. GitHub installs) are copied as built package directories, since a spec
# alone cannot reinstall them reliably.
set -euo pipefail
conda="$1"; out="$2"; shift 2
mkdir -p "$out"
for prefix in "$@"; do
  name="$(basename "$prefix")"; d="$out/$name"; mkdir -p "$d"
  "$conda" list -p "$prefix" --explicit --md5 > "$d/conda_explicit.txt"
  "$conda" env export -p "$prefix" --no-builds > "$d/conda_env.yml"
  if [[ -x "$prefix/bin/python" ]]; then
    "$prefix/bin/python" -m pip freeze --all > "$d/pip_freeze.txt" 2>/dev/null || true
  fi
  if [[ -x "$prefix/bin/Rscript" ]]; then
    "$prefix/bin/Rscript" --vanilla -e '
      ip <- installed.packages(fields = c("Repository", "RemoteType", "RemoteSha", "RemoteUrl", "RemoteRepo", "RemoteUsername", "biocViews"))
      keep <- c("Package", "LibPath", "Version", "Priority", "Repository", "RemoteType", "RemoteSha", "RemoteUrl", "RemoteRepo", "RemoteUsername", "biocViews")
      t <- as.data.frame(ip[, keep], stringsAsFactors = FALSE)
      t$origin <- ifelse(!is.na(t$RemoteType) & !(t$RemoteType %in% c("standard", "cran", "bioc")), t$RemoteType,
                  ifelse(!is.na(t$Priority), "base",
                  ifelse(!is.na(t$Repository) & t$Repository == "CRAN", "CRAN",
                  ifelse(!is.na(t$biocViews), "Bioconductor", "other"))))
      t$biocViews <- NULL
      write.table(t, file = commandArgs(TRUE)[1], sep = "\t", quote = FALSE, row.names = FALSE)
      writeLines(capture.output(sessionInfo()), commandArgs(TRUE)[2])
    ' "$d/R_packages.tsv" "$d/R_sessionInfo.txt"
    awk -F'\t' 'NR > 1 && $NF != "CRAN" && $NF != "Bioconductor" && $NF != "base" {print $2 "/" $1}' "$d/R_packages.tsv" > "$d/R_nonstandard_packages.txt"
    if [[ -s "$d/R_nonstandard_packages.txt" ]]; then
      mkdir -p "$d/R_nonstandard_library"
      while read -r p; do cp -a "$p" "$d/R_nonstandard_library/"; done < "$d/R_nonstandard_packages.txt"
    fi
  fi
done
