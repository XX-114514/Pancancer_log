# Smoke test: each R NMF method on one simulated Poisson count matrix (genes x cells, K=4).
lib <- commandArgs(TRUE)[1]; .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({library(Matrix); library(RcppML); library(CoGAPS); library(log1pNMF); library(rliger)})
set.seed(1); G <- 500; N <- 400; K <- 4
W <- matrix(rgamma(G * K, 0.3, 1), G); H <- matrix(rgamma(K * N, 0.5, 1), K)
Y <- matrix(rpois(G * N, W %*% H * 3), G, dimnames = list(paste0("g", 1:G), paste0("c", 1:N)))
Y <- Y[rowSums(Y) > 0, ]; Ys <- as(Y, "CsparseMatrix"); X <- log1p(t(t(Y) / colSums(Y)) * 1e4)
res <- list()
tm <- function(name, expr) { t0 <- Sys.time(); r <- try(expr); ok <- !inherits(r, "try-error")
  res[[name]] <<- data.frame(method = name, ok = ok, sec = round(as.numeric(Sys.time() - t0, units = "secs"), 1)); r }
m <- tm("RcppML_mse", RcppML::nmf(X, K, seed = 1, verbose = FALSE)); print(dim(m@w))
m <- tm("RcppML_gp", RcppML::nmf(Ys, K, loss = "gp", seed = 1, verbose = FALSE)); print(dim(m@w))
m <- tm("CoGAPS", CoGAPS(as.matrix(X), nPatterns = K, nIterations = 2000, seed = 1, messages = FALSE)); print(dim(m@featureLoadings))
m <- tm("log1pNMF", fit_poisson_log1p_nmf(t(Ys), K = K, cc = 1, control = list(maxiter = 20, verbose = FALSE))); print(dim(m$LL))
m <- tm("rliger_iNMF", { a <- 1:200; lig <- createLiger(list(b1 = Ys[, a], b2 = Ys[, -a]), verbose = FALSE)
  lig <- normalize(lig, verbose = FALSE); lig <- selectGenes(lig, thresh = 0, verbose = FALSE)
  lig <- scaleNotCenter(lig, verbose = FALSE); runINMF(lig, k = K, nIteration = 10, verbose = FALSE) }); print(dim(getMatrix(m, "W")))
out <- do.call(rbind, res); print(out, row.names = FALSE)
if (!all(out$ok)) quit(status = 1)
