# CM07 胶原×DDR1 在 TCGA: 纯度校正 + 表达匹配基因阴性对照 + ccRCC 专项。探索性。
library(survival); library(parallel)
a <- commandArgs(TRUE); RUN <- a[1]; E2 <- a[2]; set.seed(20261106)
E <- readRDS(file.path(E2, "inputs/expr.rds")); cl <- read.delim(file.path(E2, "inputs/clin.tsv"), stringsAsFactors = FALSE)
cv <- read.delim(file.path(E2, "inputs/covariates.tsv"), stringsAsFactors = FALSE); par <- read.delim(file.path(E2, "inputs/parents.tsv"), check.names = FALSE)
stopifnot(identical(cl$sample, colnames(E)), identical(cv$sample, cl$sample), identical(par$sample, cl$sample))
OUT <- file.path(RUN, "outputs")
zv <- function(v) { s <- sd(v, na.rm = TRUE); if (!is.finite(s) || s == 0) return(v * 0); (v - mean(v, na.rm = TRUE)) / s }
sc <- function(gs, i) rowMeans(apply(E[gs, i, drop = FALSE], 1, zv))  # apply 结果为 样本×基因
COL <- c("COL1A1","COL1A2","COL3A1")
canc <- names(which(table(cl$cancer) >= 50))
fit <- function(y, L, R, cov = NULL) {
  X <- cbind(L = L, R = R, LxR = L * R, cov); ok <- complete.cases(X, y)
  f <- lm.fit(cbind(1, X[ok, , drop = FALSE]), y[ok]); rdf <- sum(ok) - f$rank; s2 <- sum(f$residuals^2) / rdf
  V <- chol2inv(qr.R(f$qr)) * s2; b <- f$coefficients; se <- sqrt(diag(V))
  data.frame(term = c("L","R","LxR"), beta = b[2:4], se = se[2:4], t = b[2:4] / se[2:4], p = 2 * pt(-abs(b[2:4] / se[2:4]), rdf), n = sum(ok))
}
res <- list(); nulls <- list(); assoc <- list()
scor <- function(a, b) { o <- is.finite(a) & is.finite(b); if (sum(o) < 30) return(NA_real_); cor(a[o], b[o], method = "s") }
for (ca in canc) {
  i <- which(cl$cancer == ca); L <- zv(sc(COL, i)); R <- zv(E["DDR1", i])
  ys <- list(cyt = zv(sc(c("GZMA","PRF1"), i)), cd8 = zv(sc(c("CD8A","CD8B","GZMK"), i)))
  pur <- cv$purity_absolute[i]; str <- zv(cv$estimate_stromal[i]); leu <- cv$leukocyte_fraction[i]
  assoc[[ca]] <- data.frame(cancer = ca, n = length(i), r_DDR1_purity = scor(R, pur), r_DDR1_stromal = scor(R, str), r_DDR1_leuk = scor(R, leu),
    r_DDR1_col = cor(R, L, method = "s"), r_DDR1_EPCAM = cor(R, E["EPCAM", i], method = "s"), r_DDR1_LRP2 = cor(R, E["LRP2", i], method = "s"), r_cyt_purity = scor(ys$cyt, pur), n_purity = sum(!is.na(pur)))
  has_pur <- sum(is.finite(pur)) >= 50
  for (yn in names(ys)) for (mn in if (has_pur) c("raw","purity","purity_stromal") else "raw") {
    cov <- switch(mn, raw = NULL, purity = cbind(pur = zv(pur)), purity_stromal = cbind(pur = zv(pur), str = str))
    r <- fit(ys[[yn]], L, R, cov); r$cancer <- ca; r$y <- yn; r$model <- mn; res[[length(res) + 1]] <- r
  }
  # 阴性对照: DDR1 换成该癌种内表达均值和 SD 最近的 500 个基因中随机 300 个
  ex <- E[, i, drop = FALSE]; m <- rowMeans(ex); s <- apply(ex, 1, sd); ok <- s > 0 & !(rownames(ex) %in% c("DDR1", COL))
  fv <- scale(cbind(m[ok], s[ok])); tgt <- (c(m["DDR1"], s["DDR1"]) - attr(fv, "scaled:center")) / attr(fv, "scaled:scale")
  dd <- sqrt((fv[, 1] - tgt[1])^2 + (fv[, 2] - tgt[2])^2); pool <- names(sort(dd))[1:500]; pick <- sample(pool, 300)
  covp <- if (has_pur) cbind(pur = zv(pur)) else NULL
  for (yn in names(ys)) {
    z <- t(sapply(pick, function(g) { r <- fit(ys[[yn]], L, zv(ex[g, ]), covp); c(tR = r$t[2], tLxR = r$t[3]) }))
    ob <- fit(ys[[yn]], L, R, covp)
    nulls[[length(nulls) + 1]] <- data.frame(cancer = ca, y = yn, obs_tR = ob$t[2], obs_tLxR = ob$t[3], null_med_tR = median(z[, 1]), null_med_tLxR = median(z[, 2]),
      p_R_neg = (1 + sum(z[, 1] <= ob$t[2])) / 301, p_LxR_neg = (1 + sum(z[, 2] <= ob$t[3])) / 301)
  }
}
res <- do.call(rbind, res); res$q <- ave(res$p, res$term, res$y, res$model, FUN = function(p) p.adjust(p, "BH"))
write.table(res, file.path(OUT, "TCGA_receiver_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
nl <- do.call(rbind, nulls); nl$q_LxR <- ave(nl$p_LxR_neg, nl$y, FUN = function(p) p.adjust(p, "BH")); nl$q_R <- ave(nl$p_R_neg, nl$y, FUN = function(p) p.adjust(p, "BH"))
write.table(nl, file.path(OUT, "TCGA_receiver_gene_null.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
as <- do.call(rbind, assoc); write.table(as, file.path(OUT, "TCGA_DDR1_associations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# 生存: 逐癌种 OS ~ 胶原 + DDR1 + 交互 + 年龄; 再加分期和纯度
sv <- list(); stg <- ifelse(is.na(cv$stage), "NA", cv$stage)
for (ca in canc) { i <- which(cl$cancer == ca); L <- zv(sc(COL, i)); R <- zv(E["DDR1", i]); d <- data.frame(t = cl$OS.time[i], e = cl$OS[i], L, R, age = cl$age[i], stage = stg[i], pur = cv$purity_absolute[i], desmo = par$desmoplasia_myCAF[i], cm07 = par$CM07[i]); d <- d[d$t > 0 & is.finite(d$t) & is.finite(d$e), ]
  if (sum(d$e) < 20) next
  for (mn in c("base","stage","stage_purity","desmo")) { fm <- switch(mn, base = Surv(t,e) ~ L*R + age, stage = Surv(t,e) ~ L*R + age + stage, stage_purity = if (sum(is.finite(d$pur)) >= 50) Surv(t,e) ~ L*R + age + stage + pur else Surv(t,e) ~ L*R + age + stage, desmo = Surv(t,e) ~ L*R + age + desmo)
    f <- tryCatch(suppressWarnings(coxph(fm, d)), error = function(e) NULL); if (is.null(f)) next; ss <- summary(f)$coefficients
    for (tm in c("L","R","L:R")) sv[[length(sv) + 1]] <- data.frame(cancer = ca, model = mn, term = tm, HR = ss[tm, 2], lo = exp(ss[tm, 1] - 1.96 * ss[tm, 3]), hi = exp(ss[tm, 1] + 1.96 * ss[tm, 3]), p = ss[tm, 5], n = f$n, events = f$nevent) } }
sv <- do.call(rbind, sv); sv$q <- ave(sv$p, sv$term, sv$model, FUN = function(p) p.adjust(p, "BH")); write.table(sv, file.path(OUT, "TCGA_survival_by_cancer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("done\n")
