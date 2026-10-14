# Load the verified private runtime without changing the shared R library.
# PANCANCER_DIALOGUE_RUNTIME_ROOT points to the method pilot directory.
root <- Sys.getenv('PANCANCER_DIALOGUE_RUNTIME_ROOT')
if (root=='')stop('PANCANCER_DIALOGUE_RUNTIME_ROOT is required')
private <- file.path(root,c('library_v7','library_v5','library_v2','library'))
if (!all(dir.exists(private)))stop('Verified private runtime is incomplete')
.libPaths(c(private,.libPaths()))
if (!startsWith(normalizePath(find.package('lme4')),normalizePath(private[2])))stop('Compatible private lme4 did not resolve')
required <- c('DIALOGUE','lmerTest','PMA','nnls','assertthat')
if (!all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))stop('DIALOGUE runtime dependency missing')
if (as.character(packageVersion('lme4'))!='1.1.35.5')stop('lme4 runtime version differs from acceptance')
suppressPackageStartupMessages(library(DIALOGUE))
