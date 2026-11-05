# One R NMF method, one seed, K CMs on the shared CM-benchmark inputs.
# Usage: Rscript run_r_method.R <run_dir> <method> <seed> <K> [overlay_lib]
# Writes <run_dir>/methods/<method>/seed_<seed>/{W.tsv (nodes x K), H.tsv (K x samples), info.json}
a <- commandArgs(TRUE); run <- a[1]; method <- a[2]; seed <- as.integer(a[3]); K <- as.integer(a[4])
if (length(a) >= 5) .libPaths(c(a[5], .libPaths()))
X <- as.matrix(read.delim(file.path(run, "inputs/X_minmax.tsv"), row.names = 1, check.names = FALSE))
C <- as.matrix(read.delim(file.path(run, "inputs/counts.tsv"), row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(run, "inputs/metadata.tsv"), check.names = FALSE)
stopifnot(identical(meta$sampleID, colnames(X)), identical(colnames(C), colnames(X)))
out <- file.path(run, "methods", method, sprintf("seed_%d", seed)); if (dir.exists(out)) stop("Output already exists; use a new run directory"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
nrun <- 30L  # best-of-30 restarts where a method has no restarts of its own (as CoVarNet nsNMF nrun=30)
t0 <- Sys.time(); info <- list(method = method, seed = seed, K = K)
if (method == "NMF_nsNMF") {
  # exact CoVarNet call (unit weights -> weighted == X)
  suppressPackageStartupMessages(library(NMF))
  m <- nmf(X, K, method = "nsNMF", nrun = nrun, seed = seed, .options = "-p", maxIter = 2000L)
  W <- basis(m); H <- coef(m); info$objective <- deviance(m); info$input <- "X_minmax"
  info$package <- paste("NMF", packageVersion("NMF"))
} else if (method == "RcppML_mse") {
  suppressPackageStartupMessages(library(RcppML))
  best <- NULL
  for (r in seq_len(nrun)) {
    m <- RcppML::nmf(X, K, seed = seed * 100L + r, tol = 1e-6, maxit = 1000, verbose = FALSE)
    if (is.null(best) || m@misc$loss < best@misc$loss) best <- m
  }
  W <- best@w %*% diag(best@d); H <- best@h; info$objective <- best@misc$loss; info$input <- "X_minmax"
  info$package <- paste("RcppML", packageVersion("RcppML"))
} else if (method == "CoGAPS") {
  suppressPackageStartupMessages(library(CoGAPS))
  m <- CoGAPS(X, nPatterns = K, nIterations = 20000L, seed = seed, messages = FALSE, nThreads = 1L,
              checkpointOutFile = file.path(out, "gaps_checkpoint.out"))
  W <- m@featureLoadings; H <- t(m@sampleFactors); info$objective <- min(m@metadata$chisq); info$input <- "X_minmax"
  info$nIterations <- 20000L; info$package <- paste("CoGAPS", packageVersion("CoGAPS"))
} else if (method == "LIGER_iNMF") {
  # iNMF core of rliger (RcppPlanc::inmf), one dataset per source with >=30 samples, the rest pooled
  # (at >=10 most starts abort on singular solves). X is already min-max scaled, so rliger's
  # normalize/selectGenes/scaleNotCenter are not applied. Best of 5 successful starts.
  min_n <- 30L; fits <- list(); tried <- 0L
  while (length(fits) < 5L && tried < 20L) {
    tried <- tried + 1L; f <- file.path(out, sprintf("start_%02d.rds", tried))
    rc <- system2(file.path(R.home("bin"), "Rscript"), c(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "liger_start.R"),
                  run, K, seed * 100L + tried, min_n, f, a[5]), stdout = FALSE, stderr = FALSE)
    if (rc == 0 && file.exists(f)) { fits[[length(fits) + 1L]] <- readRDS(f); unlink(f) }
  }
  stopifnot(length(fits) > 0)
  best <- fits[[which.min(vapply(fits, function(x) x$m$objErr, 0))]]; groups <- best$groups; best <- best$m
  W <- best$W; H <- matrix(0, K, ncol(X))
  for (g in names(groups)) H[, groups[[g]]] <- t(best$H[[match(g, names(groups))]])
  info$objective <- best$objErr; info$input <- "X_minmax split by source"; info$lambda <- 5; info$min_samples_per_dataset <- min_n
  info$starts_tried <- tried; info$starts_ok <- length(fits)
  info$datasets <- as.list(lengths(groups)); info$package <- paste("RcppPlanc", packageVersion("RcppPlanc"))
} else if (method == "log1pNMF_counts") {
  # count-native model: raw node cell counts, size factor = sample total (package default)
  suppressPackageStartupMessages({library(log1pNMF); library(Matrix)})
  set.seed(seed)
  m <- fit_poisson_log1p_nmf(as(t(C), "CsparseMatrix"), K = K, cc = 1, init_method = "random",
                             control = list(maxiter = 500, verbose = FALSE, threads = 1))
  W <- m$FF; H <- t(m$LL); info$input <- "counts"
  info$package <- paste("log1pNMF", packageVersion("log1pNMF"))
} else stop("unknown method ", method)
info$seconds <- as.numeric(Sys.time() - t0, units = "secs")
W <- as.matrix(W); H <- as.matrix(H)
dimnames(W) <- list(rownames(X), sprintf("F%02d", seq_len(K))); dimnames(H) <- list(colnames(W), colnames(X))
write.table(W, file.path(out, "W.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(H, file.path(out, "H.tsv"), sep = "\t", quote = FALSE, col.names = NA)
jsonlite::write_json(info, file.path(out, "info.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(method, seed, "done", round(info$seconds, 1), "s\n")
