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
rows <- list()
for (p in order) {
  v <- ver[p]
  have <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  if (is.na(v)) { rows[[p]] <- c(p, "", have, "not_in_old_env"); next }
  if (same(have, v)) { rows[[p]] <- c(p, v, have, "already_exact"); next }
  ok <- tryCatch({
    if (grepl("Bioc", origin[p])) BiocManager::install(p, version = "3.18", update = FALSE, ask = FALSE, dependencies = FALSE)
    else remotes::install_version(p, version = v, dependencies = FALSE, upgrade = "never", quiet = FALSE)
    TRUE
  }, error = function(e) { message(p, ": ", conditionMessage(e)); FALSE })
  have <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  rows[[p]] <- c(p, v, have, if (isTRUE(ok) && same(have, unname(v))) "installed_exact" else "MISMATCH")
}
for (p in intersect(basepk, old$Package)) {
  have <- as.character(packageVersion(p))
  if (!same(have, ver[p])) rows[[paste0("base:", p)]] <- c(p, ver[p], have, "base_version_differs")
}
out <- do.call(rbind, rows); colnames(out) <- c("package", "old", "wehi", "result")
write.table(out, report, sep = "\t", quote = FALSE, row.names = FALSE)
print(table(out[, "result"]))
