# Install R packages at the exact versions recorded on the old server (R_packages.tsv), from
# the CRAN archive / Bioconductor, resolving the recursive Depends/Imports/LinkingTo closure.
# Usage: Rscript wehi_install_R_exact.R R_packages.tsv pkg1:pkg2:... REPORT.tsv  (":" or ",")
# Base/recommended packages that ship with r-base are left alone; a version mismatch there is
# reported, not changed. Nothing is upgraded beyond the recorded version.
args <- commandArgs(trailingOnly = TRUE)
old <- read.delim(args[[1]], stringsAsFactors = FALSE)
want <- strsplit(args[[2]], "[,:]")[[1]]
report <- args[[3]]
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
ver <- setNames(old$Version, old$Package)
origin <- setNames(old$origin, old$Package)
bioc <- BiocManager::repositories(version = "3.18")
db <- available.packages(repos = c(bioc, CRAN = "https://cloud.r-project.org"))
deps <- tools::package_dependencies(want, db = db, recursive = TRUE, which = c("Depends", "Imports", "LinkingTo"))
need <- unique(c(unlist(deps), want))
basepk <- rownames(installed.packages(priority = c("base", "recommended")))
need <- setdiff(need, c(basepk, "R"))
# topological order over the closure
direct <- tools::package_dependencies(need, db = db, which = c("Depends", "Imports", "LinkingTo"))
order <- character(0)
visit <- function(p) {
  if (p %in% order || !(p %in% need)) return(invisible())
  for (d in direct[[p]]) visit(d)
  order <<- c(order, p)
}
for (p in need) visit(p)
same <- function(a, b) !is.na(a) && !is.na(b) && package_version(a) == package_version(b)
# Up to three passes: when the R 4.3 package index lacks a package's current entry (e.g. a newer
# Matrix needs a newer R), its dependencies cannot be ordered first; retrying resolves that.
rows <- list()
todo <- order
for (pass in 1:3) {
  failed <- character(0)
  for (p in todo) {
    v <- ver[p]
    have <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
    if (is.na(v)) { rows[[p]] <- c(p, "", have, "not_in_old_env"); next }
    if (same(have, v)) { if (is.null(rows[[p]])) rows[[p]] <- c(p, v, have, "already_exact"); next }
    ok <- tryCatch({
      if (grepl("Bioc", origin[p])) {
        # BiocManager only offers the latest 3.18 patch; try the exact recorded tarball first.
        base <- "https://bioconductor.org/packages/3.18"
        urls <- c(sprintf("%s/%s/src/contrib/%s_%s.tar.gz", base, c("bioc", "data/annotation", "data/experiment"), p, v),
                  sprintf("%s/%s/src/contrib/Archive/%s/%s_%s.tar.gz", base, c("bioc", "data/annotation", "data/experiment"), p, p, v))
        done <- FALSE
        for (u in urls) {
          f <- tempfile(fileext = ".tar.gz")
          if (!inherits(try(suppressWarnings(download.file(u, f, quiet = TRUE)), silent = TRUE), "try-error")) {
            install.packages(f, repos = NULL, type = "source"); done <- TRUE; break
          }
        }
        if (!done) BiocManager::install(p, version = "3.18", update = FALSE, ask = FALSE, dependencies = FALSE)
      }
      else remotes::install_version(p, version = v, dependencies = FALSE, upgrade = "never", quiet = FALSE)
      TRUE
    }, error = function(e) { message(p, ": ", conditionMessage(e)); FALSE })
    have <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
    exact <- isTRUE(ok) && same(have, unname(v))
    rows[[p]] <- c(p, v, have, if (exact) "installed_exact" else "MISMATCH")
    if (!exact) failed <- c(failed, p)
  }
  if (!length(failed)) break
  todo <- failed
}
for (p in intersect(basepk, old$Package)) {
  have <- as.character(packageVersion(p))
  if (!same(have, ver[p])) rows[[paste0("base:", p)]] <- c(p, ver[p], have, "base_version_differs")
}
out <- do.call(rbind, rows); colnames(out) <- c("package", "old", "wehi", "result")
write.table(out, report, sep = "\t", quote = FALSE, row.names = FALSE)
print(table(out[, "result"]))
