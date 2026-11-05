required <- c("blisa", "sf", "spdep", "fastLISA", "Matrix", "SpatialExperiment",
              "SummarizedExperiment", "ComplexHeatmap", "ggplot2", "viridisLite")
cat(R.version.string, "\n")
status <- vapply(required, function(p) requireNamespace(p, quietly=TRUE), logical(1))
for (p in required) cat(p, if (status[p]) as.character(packageVersion(p)) else "MISSING", "\n")
if (!all(status)) quit(status=2)
if (as.character(packageVersion("blisa")) != "1.0.3") stop("Runner targets blisa 1.0.3; audit API before changing version")
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")
