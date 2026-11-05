# Install the R NMF methods into an overlay library on top of scrna_r_full.
# Usage: Rscript install_nmf_r.R <overlay_lib> <src_dir>
a <- commandArgs(TRUE); lib <- a[1]; src <- a[2]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = BiocManager::repositories(), Ncpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
need <- c("fgsea", "pracma", "distr", "leidenAlg", "RcppPlanc")
# installed.packages(), not requireNamespace(): a broken build can segfault on load.
# REINSTALL=a:b forces those packages to be rebuilt.
redo <- strsplit(Sys.getenv("REINSTALL"), ":")[[1]]
need <- union(setdiff(need, rownames(installed.packages())), intersect(redo, need))
# fgsea 1.28 (Bioc 3.18) does not compile against BH 1.87 (boost::math::pow constexpr clash):
# build it against BH 1.84 kept in a throwaway lib; BH is header-only, so nothing links to it at runtime.
if ("fgsea" %in% need) {
  bh <- file.path(tempdir(), "bh184"); dir.create(bh)
  remotes::install_version("BH", "1.84.0-0", lib = bh, upgrade = "never")
  withr::with_libpaths(c(bh, .libPaths()), BiocManager::install("fgsea", lib = lib, update = FALSE, ask = FALSE))
  need <- setdiff(need, "fgsea")
}
if (length(need)) BiocManager::install(need, lib = lib, update = FALSE, ask = FALSE, force = TRUE)
# GitHub HEADs cloned under src/ (commits recorded in versions.tsv)
gh <- c(RcppML = "RcppML", CoGAPS = "CoGAPS", log1pNMF = "log1pNMF", rliger = "liger")
gh <- gh[!names(gh) %in% rownames(installed.packages(lib)) | names(gh) %in% redo]
for (p in gh)
  remotes::install_local(file.path(src, p), lib = lib, dependencies = TRUE, upgrade = "never",
                         force = TRUE, build_vignettes = FALSE, INSTALL_opts = "--no-test-load")
pk <- c("RcppML", "CoGAPS", "log1pNMF", "rliger", "RcppPlanc", "fgsea")
# load each package in its own R process so one crash cannot hide the others
rs <- file.path(R.home("bin"), "Rscript")
ok <- vapply(pk, function(p) system2(rs, c("-e", shQuote(sprintf(
  ".libPaths(c('%s', .libPaths())); suppressPackageStartupMessages(library(%s))", lib, p)))) == 0, TRUE)
v <- data.frame(package = pk, version = vapply(pk, function(p) as.character(packageVersion(p)), ""),
                libpath = vapply(pk, function(p) dirname(find.package(p)), ""), loads = ok)
print(v, row.names = FALSE)
write.table(v, file.path(lib, "R_versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(ok)) stop("failed to load: ", paste(pk[!ok], collapse = ", "))
