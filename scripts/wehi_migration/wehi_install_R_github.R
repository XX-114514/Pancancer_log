# Install the old env's GitHub/other R packages at the recorded commit (R_packages.tsv RemoteSha);
# without a recorded commit, or if the build fails, copy the old compiled package from the
# exported R_nonstandard_library and record that it is a copy.
# Usage: Rscript wehi_install_R_github.R R_packages.tsv R_nonstandard_library REPORT.tsv
args <- commandArgs(trailingOnly = TRUE)
old <- read.delim(args[[1]], stringsAsFactors = FALSE)
libdir <- args[[2]]; report <- args[[3]]
gh <- old[old$origin %in% c("github", "other"), ]
lib <- .libPaths()[1]
rows <- list()
for (i in seq_len(nrow(gh))) {
  p <- gh$Package[i]; v <- gh$Version[i]; sha <- gh$RemoteSha[i]
  repo <- if (!is.na(gh$RemoteUsername[i])) paste0(gh$RemoteUsername[i], "/", gh$RemoteRepo[i]) else NA
  how <- "failed"
  if (!is.na(repo) && !is.na(sha)) {
    ok <- tryCatch({ remotes::install_github(repo, ref = sha, dependencies = FALSE, upgrade = "never", force = TRUE); TRUE },
                   error = function(e) { message(p, ": ", conditionMessage(e)); FALSE })
    if (ok && requireNamespace(p, quietly = TRUE)) how <- "github_exact_commit"
  }
  if (how == "failed" && dir.exists(file.path(libdir, p))) {
    unlink(file.path(lib, p), recursive = TRUE)
    file.copy(file.path(libdir, p), lib, recursive = TRUE)
    how <- if (requireNamespace(p, quietly = TRUE)) "copied_old_build" else "copied_old_build_LOAD_FAILED"
  }
  have <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  rows[[p]] <- c(p, v, have, ifelse(is.na(sha), "", sha), how)
}
out <- do.call(rbind, rows); colnames(out) <- c("package", "old", "wehi", "old_sha", "result")
write.table(out, report, sep = "\t", quote = FALSE, row.names = FALSE)
print(out)
