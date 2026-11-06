#!/usr/bin/env Rscript
# E2 | 轴增量的零分布与混杂校正
# 用法: Rscript e2_axis_increment.R <config.R> <step>
# step: repro | marginal | family | null | robust | decide | all
# 依赖: R 基础包 + survival + parallel
# 规则: repro 未通过时其余步骤拒绝运行; 输出目录内已有结果不覆盖。

suppressPackageStartupMessages({library(survival); library(parallel)})
args <- commandArgs(TRUE)
if (length(args) < 2) stop("usage: e2_axis_increment.R <config.R> <step>")
cfg <- new.env(); sys.source(args[1], envir = cfg); cfg <- as.list(cfg)
step <- args[2]
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(cfg$out_dir, f)
log_msg <- function(...) cat(format(Sys.time(), "%FT%T%z"), ..., "\n")

write_new <- function(x, f) {
  p <- out(f)
  if (file.exists(p)) stop("refuse to overwrite: ", p)
  write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(p)
}
sha256 <- function(p) {
  o <- tryCatch(system2("sha256sum", shQuote(p), stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (!length(o)) o <- system2("shasum", c("-a", "256", shQuote(p)), stdout = TRUE)
  sub(" .*", "", o[1])
}

## ---------- 载入 ----------
load_all <- function() {
  log_msg("loading inputs")
  expr <- if (grepl("\\.rds$", cfg$expr)) readRDS(cfg$expr) else
    as.matrix(read.delim(cfg$expr, row.names = 1, check.names = FALSE))
  clin <- read.delim(cfg$clin, stringsAsFactors = FALSE, check.names = FALSE)
  parents <- read.delim(cfg$parents, check.names = FALSE)
  axes <- read.delim(cfg$axes, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(clin$sample), !anyDuplicated(clin$patient),
            !anyDuplicated(parents$sample), !anyDuplicated(axes$axis_id),
            !anyDuplicated(rownames(expr)))
  ids <- Reduce(intersect, list(clin$sample, parents$sample, colnames(expr)))
  clin <- clin[match(ids, clin$sample), ]
  parents <- parents[match(ids, parents$sample), ]
  expr <- expr[, ids, drop = FALSE]
  clin$age_z <- zby(clin$age, rep(1, nrow(clin)))
  # 基因在癌种内 z 标准化; 癌种内零方差的基因该癌种记 0
  log_msg("within-cancer z:", nrow(expr), "genes x", ncol(expr), "samples")
  Z <- expr
  for (ca in unique(clin$cancer)) {
    i <- which(clin$cancer == ca)
    m <- rowMeans(expr[, i, drop = FALSE], na.rm = TRUE)
    s <- sqrt(rowSums((expr[, i, drop = FALSE] - m)^2, na.rm = TRUE) / (length(i) - 1))
    s[!is.finite(s) | s == 0] <- Inf
    Z[, i] <- (expr[, i, drop = FALSE] - m) / s
  }
  list(expr = expr, Z = Z, clin = clin, parents = parents, axes = axes)
}
zby <- function(x, g) {
  x <- as.numeric(x)
  ave(x, g, FUN = function(v) {
    s <- sd(v, na.rm = TRUE); if (!is.finite(s) || s == 0) return(v * 0)
    (v - mean(v, na.rm = TRUE)) / s
  })
}
genes_of <- function(s) trimws(strsplit(s, "[,;]")[[1]])
score_set <- function(Z, g) {
  g <- intersect(g, rownames(Z))
  if (!length(g)) return(rep(NA_real_, ncol(Z)))
  colMeans(Z[g, , drop = FALSE])
}
rescale_fun <- function(v, cancer, mode) switch(mode,
  none = v,
  global_sd = v / sd(v, na.rm = TRUE),
  within_cancer_z = zby(v, cancer))
proxy_of <- function(L, R, form, cancer, mode) {
  v <- switch(form, min = pmin(L, R), mean = (L + R) / 2, prod = L * R)
  rescale_fun(v, cancer, mode)
}

## ---------- Cox ----------
make_y <- function(clin, ep) {
  tt <- as.numeric(clin[[paste0(ep, ".time")]]); ee <- as.numeric(clin[[ep]])
  list(t = tt, e = ee, ok = is.finite(tt) & is.finite(ee) & tt > 0)
}
cox_fit <- function(X, y, strat) {
  # X: 已无缺失的数值矩阵; strat: 整数层
  fit <- suppressWarnings(survival::coxph.fit(
    x = X, y = Surv(y$t, y$e), strata = as.integer(factor(strat)),
    offset = NULL, init = NULL, control = survival::coxph.control(),
    weights = NULL, method = "efron", rownames = NULL, resid = FALSE))
  list(beta = fit$coefficients, se = sqrt(diag(as.matrix(fit$var))),
       ll = fit$loglik[2], n = nrow(X), events = sum(y$e))
}
fit_cols <- function(cols, y, strat) {
  X <- do.call(cbind, cols)
  keep <- y$ok & stats::complete.cases(X)
  r <- cox_fit(X[keep, , drop = FALSE],
               list(t = y$t[keep], e = y$e[keep]), strat[keep])
  names(r$beta) <- names(r$se) <- colnames(X); r$keep <- keep; r
}
hr <- function(f, nm) c(HR = exp(f$beta[[nm]]), lo = exp(f$beta[[nm]] - 1.96 * f$se[[nm]]),
                        hi = exp(f$beta[[nm]] + 1.96 * f$se[[nm]]),
                        p = 2 * pnorm(-abs(f$beta[[nm]] / f$se[[nm]])))
# 模型 A/C 必须在同一行集合上比较
model_ac <- function(par, prox, age, y, strat) {
  ok <- y$ok & is.finite(par) & is.finite(prox) & is.finite(age)
  yy <- list(t = y$t[ok], e = y$e[ok]); st <- strat[ok]
  A <- cox_fit(cbind(parent = par[ok], age = age[ok]), yy, st)
  C <- cox_fit(cbind(parent = par[ok], proxy = prox[ok], age = age[ok]), yy, st)
  names(A$beta) <- names(A$se) <- c("parent", "age")
  names(C$beta) <- names(C$se) <- c("parent", "proxy", "age")
  list(A = A, C = C, lrt = 2 * (C$ll - A$ll),
       lrt_p = pchisq(pmax(2 * (C$ll - A$ll), 0), 1, lower.tail = FALSE))
}

gate <- function() {
  f <- out("repro_gate.tsv")
  if (!file.exists(f)) stop("repro has not been run; no new effects may be computed")
  g <- read.delim(f)
  if (!all(g$pass)) stop("repro gate FAILED; no new effects may be computed. See ", f)
  readLines(out("repro_mode.txt"))[1]
}

## ---------- step: repro ----------
step_repro <- function(d) {
  pub <- read.delim(cfg$published, stringsAsFactors = FALSE)
  res <- list()
  for (mode in c("none", "global_sd", "within_cancer_z")) {
    for (i in seq_len(nrow(pub))) {
      a <- d$axes[d$axes$axis_id == pub$axis[i], ]
      if (!nrow(a)) stop("published axis missing from axis table: ", pub$axis[i])
      L <- score_set(d$Z, genes_of(a$ligand_genes)); R <- score_set(d$Z, genes_of(a$receptor_genes))
      form <- sub("_.*", "", pub$proxy[i])
      px <- proxy_of(L, R, form, d$clin$cancer, mode)
      m <- model_ac(d$parents[[a$parent]], px, d$clin$age_z, make_y(d$clin, pub$endpoint[i]), d$clin$cancer)
      res[[length(res) + 1]] <- data.frame(mode = mode, axis = pub$axis[i], proxy = pub$proxy[i],
        endpoint = pub$endpoint[i], HR_proxy = exp(m$C$beta[["proxy"]]), HR_pub = pub$HR_proxy_in_C[i],
        HR_parent = exp(m$C$beta[["parent"]]), HR_parent_pub = pub$HR_parent_in_C[i], n = m$C$n,
        n_pub = pub$n[i], d_proxy = abs(m$C$beta[["proxy"]] - log(pub$HR_proxy_in_C[i])),
        d_parent = abs(m$C$beta[["parent"]] - log(pub$HR_parent_in_C[i])))
    }
  }
  res <- do.call(rbind, res)
  res$pass <- res$d_proxy <= cfg$repro_tol & res$d_parent <= cfg$repro_tol & res$n == res$n_pub
  by_mode <- tapply(res$pass, res$mode, all)
  best <- if (any(by_mode)) names(which(by_mode))[1] else NA
  log_msg("repro pass rate by proxy scaling:", paste(names(by_mode), tapply(res$pass, res$mode, mean), collapse = "; "))
  write_new(res, "repro_all_modes.tsv")
  # 诊断: 模块单独 HR 在两张公开表间不一致, 这里并列是否校正年龄的结果
  ep1 <- unique(pub$endpoint)
  diag <- do.call(rbind, lapply(unique(d$axes$parent[d$axes$axis_id %in% pub$axis]), function(pn) do.call(rbind, lapply(ep1, function(ep) {
    y <- make_y(d$clin, ep); st <- d$clin$cancer
    f1 <- fit_cols(list(parent = cbind(parent = d$parents[[pn]]), age = cbind(age = d$clin$age_z)), y, st)
    ok <- is.finite(d$clin$age_z)
    y0 <- y; y0$ok <- y$ok & ok
    f0 <- fit_cols(list(parent = cbind(parent = d$parents[[pn]])), y0, st)
    y1 <- y; f2 <- fit_cols(list(parent = cbind(parent = d$parents[[pn]])), y1, st)
    data.frame(parent = pn, endpoint = ep, HR_age_adjusted = exp(f1$beta[["parent"]]), HR_no_age_same_rows = exp(f0$beta[["parent"]]),
               HR_no_age_all_rows = exp(f2$beta[["parent"]]), n_age = f1$n, n_all = f2$n)
  }))))
  write_new(diag, "repro_parent_alone_diagnostic.tsv")
  g <- if (is.na(best)) res[res$mode == "none", ] else res[res$mode == best, ]
  if (is.na(best)) g$pass <- FALSE
  write_new(g, "repro_gate.tsv")
  writeLines(ifelse(is.na(best), "NONE", best), out("repro_mode.txt"))
  if (is.na(best)) stop("repro FAILED under all scalings: stop and reconcile with the original script")
  log_msg("repro PASS with scaling =", best)
}

## ---------- step: marginal ----------
# 每个候选基因的表达均值、癌种内 SD、单基因分层 Cox z (OS), 供匹配用
step_marginal <- function(d) {
  gate()
  expr <- d$expr; cl <- d$clin
  keep <- rowMeans(expr > cfg$expr_floor, na.rm = TRUE) >= cfg$min_detect_frac
  genes <- rownames(expr)[keep]
  y <- make_y(cl, cfg$primary_endpoint); st <- cl$cancer
  log_msg("marginal Cox for", length(genes), "genes")
  zs <- unlist(mclapply(genes, function(g) {
    f <- tryCatch(fit_cols(list(g = cbind(g = d$Z[g, ]), age = cbind(age = cl$age_z)), y, st), error = function(e) NULL)
    if (is.null(f)) NA_real_ else f$beta[[1]] / f$se[[1]]
  }, mc.cores = cfg$cores))
  m <- rowMeans(expr[genes, , drop = FALSE])
  s <- sapply(split(seq_along(cl$cancer), cl$cancer), function(i) apply(expr[genes, i, drop = FALSE], 1, sd))
  write_new(data.frame(gene = genes, mean_expr = m, within_sd = rowMeans(s), cox_z = zs), "gene_marginals.tsv")
}

## ---------- step: family ----------
family_row <- function(a, form, ep, d, mode) {
  L <- score_set(d$Z, genes_of(a$ligand_genes)); R <- score_set(d$Z, genes_of(a$receptor_genes))
  if (anyNA(L) || anyNA(R) || !a$parent %in% names(d$parents)) return(NULL)
  px <- proxy_of(L, R, form, d$clin$cancer, mode)
  y <- make_y(d$clin, ep)
  m <- model_ac(d$parents[[a$parent]], px, d$clin$age_z, y, d$clin$cancer)
  fB <- fit_cols(list(proxy = cbind(proxy = px), age = cbind(age = d$clin$age_z)), y, d$clin$cancer)
  data.frame(axis = a$axis_id, tier = a$family, parent = a$parent, form = form, endpoint = ep,
    n = m$C$n, events = m$C$events, z_proxy = m$C$beta[["proxy"]] / m$C$se[["proxy"]],
    HR_proxy = exp(m$C$beta[["proxy"]]), lo = exp(m$C$beta[["proxy"]] - 1.96 * m$C$se[["proxy"]]),
    hi = exp(m$C$beta[["proxy"]] + 1.96 * m$C$se[["proxy"]]),
    HR_parent_A = exp(m$A$beta[["parent"]]), HR_parent_C = exp(m$C$beta[["parent"]]),
    lrt = m$lrt, lrt_p = m$lrt_p, HR_proxy_B = exp(fB$beta[["proxy"]]),
    p_proxy_B = 2 * pnorm(-abs(fB$beta[["proxy"]] / fB$se[["proxy"]])))
}
step_family <- function(d) {
  mode <- gate()
  grid <- expand.grid(i = seq_len(nrow(d$axes)), form = cfg$proxy_forms, ep = cfg$endpoints, stringsAsFactors = FALSE)
  res <- do.call(rbind, mclapply(seq_len(nrow(grid)), function(k)
    family_row(d$axes[grid$i[k], ], grid$form[k], grid$ep[k], d, mode), mc.cores = cfg$cores))
  # BH 家族: 每个 终点×形式 内的全部轴
  res$lrt_q <- ave(res$lrt_p, res$form, res$endpoint, FUN = function(p) p.adjust(p, "BH"))
  write_new(res, "family_all_axes.tsv")
  # N1: 配体—受体资源中的全部真实配对, 分别以 2 个主轴的父程序为条件
  lr <- read.delim(cfg$lr_resource, stringsAsFactors = FALSE)
  prim <- d$axes[d$axes$family == "primary", ]
  n1 <- do.call(rbind, lapply(seq_len(nrow(prim)), function(p) {
    do.call(rbind, mclapply(seq_len(nrow(lr)), function(k) {
      a <- data.frame(axis_id = paste(lr$ligand_genes[k], lr$receptor_genes[k], sep = "|"), family = "N1",
                      parent = prim$parent[p], ligand_genes = lr$ligand_genes[k],
                      receptor_genes = lr$receptor_genes[k], stringsAsFactors = FALSE)
      r <- family_row(a, "min", cfg$primary_endpoint, d, mode)
      if (!is.null(r)) r$target <- prim$axis_id[p]
      r
    }, mc.cores = cfg$cores))
  }))
  write_new(n1, "null_N1_real_pairs.tsv")
}

## ---------- step: null ----------
# 匹配抽样: 每个基因从特征空间 K 近邻中抽; 排除轴基因、父程序基因、资源内已知配对基因对
knn_pool <- function(g, M, feats, K, exclude) {
  if (!g %in% M$gene) stop("gene not in marginal table: ", g)
  X <- scale(as.matrix(M[, feats, drop = FALSE]))
  rownames(X) <- M$gene
  d <- sqrt(colSums((t(X) - X[g, ])^2))
  d[names(d) %in% exclude] <- Inf
  names(sort(d))[seq_len(K)]
}
step_null <- function(d) {
  mode <- gate()
  M <- read.delim(out("gene_marginals.tsv"), stringsAsFactors = FALSE)
  M <- M[stats::complete.cases(M), ]
  pg <- read.delim(cfg$parent_genes, stringsAsFactors = FALSE)
  lr <- read.delim(cfg$lr_resource, stringsAsFactors = FALSE)
  lr_keys <- unique(unlist(lapply(seq_len(nrow(lr)), function(k)
    as.vector(outer(genes_of(lr$ligand_genes[k]), genes_of(lr$receptor_genes[k]), paste)))))
  y <- make_y(d$clin, cfg$primary_endpoint); st <- d$clin$cancer; age <- d$clin$age_z
  prim <- d$axes[d$axes$family == "primary", ]
  families <- list(
    N2a = list(fix = "none", feats = c("mean_expr", "within_sd"), K = cfg$K_expr),
    N2b = list(fix = "none", feats = c("mean_expr", "within_sd", "cox_z"), K = cfg$K_prog),
    N3L = list(fix = "ligand", feats = c("mean_expr", "within_sd"), K = cfg$K_expr),
    N3R = list(fix = "receptor", feats = c("mean_expr", "within_sd"), K = cfg$K_expr))
  summ <- list()
  for (p in seq_len(nrow(prim))) {
    a <- prim[p, ]; Lg <- genes_of(a$ligand_genes); Rg <- genes_of(a$receptor_genes)
    Lg <- intersect(Lg, rownames(d$Z)); Rg <- intersect(Rg, rownames(d$Z))
    par <- d$parents[[a$parent]]
    excl <- unique(c(Lg, Rg, pg$gene[pg$signature == a$parent], cfg$null_exclude))
    obs_L <- score_set(d$Z, Lg); obs_R <- score_set(d$Z, Rg)
    obs <- model_ac(par, proxy_of(obs_L, obs_R, "min", st, mode), age, y, st)
    obs_z <- obs$C$beta[["proxy"]] / obs$C$se[["proxy"]]
    obs_dpar <- obs$C$beta[["parent"]] - obs$A$beta[["parent"]]
    sgn <- if (a$direction == "adverse") 1 else -1
    for (nm in names(families)) {
      fm <- families[[nm]]
      Lp <- lapply(Lg, knn_pool, M, fm$feats, fm$K, excl)
      Rp <- lapply(Rg, knn_pool, M, fm$feats, fm$K, excl)
      RNGkind("L'Ecuyer-CMRG"); set.seed(cfg$seed + 1000 * p + match(nm, names(families)))
      draws <- mclapply(seq_len(cfg$B), function(b) {
        for (tries in 1:50) {
          l <- if (fm$fix == "ligand") Lg else vapply(Lp, sample, "", 1)
          r <- if (fm$fix == "receptor") Rg else vapply(Rp, sample, "", 1)
          if (anyDuplicated(c(l, r))) next
          if (any(as.vector(outer(l, r, paste)) %in% lr_keys)) next
          px <- proxy_of(score_set(d$Z, l), score_set(d$Z, r), "min", st, mode)
          m <- tryCatch(model_ac(par, px, age, y, st), error = function(e) NULL)
          if (is.null(m)) return(c(z = NA, lrt = NA, dpar = NA))
          return(c(z = m$C$beta[["proxy"]] / m$C$se[["proxy"]], lrt = m$lrt,
                   dpar = m$C$beta[["parent"]] - m$A$beta[["parent"]]))
        }
        c(z = NA, lrt = NA, dpar = NA)
      }, mc.cores = cfg$cores)
      draws <- do.call(rbind, draws); fails <- sum(is.na(draws[, "z"]))
      if (fails > cfg$B * 0.01) stop("more than 1% null draws failed: ", a$axis_id, " ", nm)
      z <- draws[is.finite(draws[, "z"]), "z"]
      summ[[length(summ) + 1]] <- data.frame(axis = a$axis_id, family = nm, B = length(z), failed = fails,
        obs_z = obs_z, null_median_z = median(z), null_q95_z = quantile(z, .95 ),
        p_emp = (1 + sum(sgn * z >= sgn * obs_z)) / (length(z) + 1),
        obs_dpar = obs_dpar, frac_dpar_ge_obs = mean(abs(draws[, "dpar"]) >= abs(obs_dpar), na.rm = TRUE))
      write_new(as.data.frame(draws), sprintf("null_draws_%s_%s.tsv", a$axis_id, nm))
      log_msg(a$axis_id, nm, "p_emp =", tail(summ, 1)[[1]]$p_emp)
    }
  }
  write_new(do.call(rbind, summ), "null_summary.tsv")
}

## ---------- step: robust ----------
dl_meta <- function(b, s) {
  w <- 1 / s^2; mu <- sum(w * b) / sum(w); Q <- sum(w * (b - mu)^2); k <- length(b)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  w2 <- 1 / (s^2 + tau2); m <- sum(w2 * b) / sum(w2); se <- sqrt(1 / sum(w2))
  c(HR = exp(m), lo = exp(m - 1.96 * se), hi = exp(m + 1.96 * se),
    I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0, tau2 = tau2, k = k)
}
step_robust <- function(d) {
  mode <- gate()
  cl <- d$clin; st <- cl$cancer; age <- cl$age_z
  cv <- read.delim(cfg$covariates, check.names = FALSE, stringsAsFactors = FALSE)
  cv <- cv[match(cl$sample, cv$sample), ]
  num <- function(nms) { x <- lapply(nms, function(n) zby(cv[[n]], st)); names(x) <- nms; x }
  stage <- cv[[cfg$stage_col]]; stage[is.na(stage) | stage == ""] <- "missing"
  S <- model.matrix(~ factor(stage) - 1)[, -1, drop = FALSE]
  colnames(S) <- paste0("stage_", seq_len(ncol(S)))
  sex <- as.numeric(factor(cv[[cfg$sex_col]])) - 1
  base <- c(list(sex = sex), split(S, col(S), drop = TRUE) |> setNames(colnames(S)))
  sets <- list(
    S0_age = list(),
    S1_stage_sex = base,
    S2_plus_ESTIMATE = c(base, num(cfg$estimate_cols)),
    S3_plus_fractions = c(base, num(cfg$fraction_cols)),
    S4_full = c(base, num(c(cfg$estimate_cols, cfg$fraction_cols))),
    S5_full_plus_desmoplasia = c(base, num(c(cfg$estimate_cols, cfg$fraction_cols)),
                                 list(desmo = d$parents[[cfg$desmoplasia_col]])),
    S6_ABSOLUTE_subset = c(base, num(cfg$purity_col)))
  prim <- d$axes[d$axes$family %in% c("primary", "control"), ]
  rob <- loco <- het <- ends <- boot <- ov <- list()
  for (p in seq_len(nrow(prim))) {
    a <- prim[p, ]
    Lg <- genes_of(a$ligand_genes); Rg <- genes_of(a$receptor_genes)
    L <- score_set(d$Z, Lg); R <- score_set(d$Z, Rg)
    px <- proxy_of(L, R, "min", st, mode); par <- d$parents[[a$parent]]
    for (ep in cfg$endpoints) {
      y <- make_y(cl, ep)
      for (sn in names(sets)) {
        f <- fit_cols(c(list(proxy = cbind(proxy = px), parent = cbind(parent = par), age = cbind(age = age)),
                        lapply(sets[[sn]], function(v) as.matrix(v))), y, st)
        rob[[length(rob) + 1]] <- data.frame(axis = a$axis_id, endpoint = ep, model = sn, n = f$n, events = f$events,
          t(hr(f, "proxy")), HR_parent = exp(f$beta[["parent"]]))
      }
      # 替换父程序版本(V7)
      if (!is.null(cfg$alt_parent_map) && a$parent %in% names(cfg$alt_parent_map)) {
        alt <- d$parents[[cfg$alt_parent_map[[a$parent]]]]
        m <- model_ac(alt, px, age, y, st)
        rob[[length(rob) + 1]] <- data.frame(axis = a$axis_id, endpoint = ep, model = "S7_alt_parent_V7",
          n = m$C$n, events = m$C$events, t(hr(m$C, "proxy")), HR_parent = exp(m$C$beta[["parent"]]))
      }
      # 留一癌种
      for (ca in unique(st)) {
        y2 <- y; y2$ok <- y$ok & st != ca
        m <- model_ac(par, px, age, y2, st)
        loco[[length(loco) + 1]] <- data.frame(axis = a$axis_id, endpoint = ep, left_out = ca, t(hr(m$C, "proxy")), lrt_p = m$lrt_p)
      }
      # 逐癌种与随机效应合并
      pc <- do.call(rbind, lapply(unique(st), function(ca) {
        i <- y$ok & st == ca
        if (sum(i) < cfg$min_cancer_n || sum(y$e[i]) < cfg$min_cancer_events) return(NULL)
        f <- tryCatch(cox_fit(cbind(parent = par[i], proxy = px[i], age = age[i]),
                              list(t = y$t[i], e = y$e[i]), rep(1L, sum(i))), error = function(e) NULL)
        if (is.null(f)) return(NULL)
        data.frame(axis = a$axis_id, endpoint = ep, cancer = ca, n = f$n, events = f$events, b = f$beta[2], se = f$se[2])
      }))
      het[[length(het) + 1]] <- pc
      mm <- dl_meta(pc$b, pc$se)
      ends[[length(ends) + 1]] <- data.frame(axis = a$axis_id, endpoint = ep, term = "RE_per_cancer", t(mm),
        frac_same_dir = mean(sign(pc$b) == sign(log(mm[["HR"]]))))
      # 嵌套: 代理是否超出单端
      nest <- function(extra, lab) {
        ok <- y$ok & is.finite(px) & is.finite(par)
        base_x <- cbind(parent = par[ok], age = age[ok], extra[ok, , drop = FALSE])
        yy <- list(t = y$t[ok], e = y$e[ok])
        f0 <- cox_fit(base_x, yy, st[ok]); f1 <- cox_fit(cbind(base_x, proxy = px[ok]), yy, st[ok])
        ends[[length(ends) + 1]] <<- data.frame(axis = a$axis_id, endpoint = ep, term = lab,
          HR = exp(tail(f1$beta, 1)), lo = exp(tail(f1$beta, 1) - 1.96 * tail(f1$se, 1)),
          hi = exp(tail(f1$beta, 1) + 1.96 * tail(f1$se, 1)), I2 = NA, tau2 = NA, k = NA,
          frac_same_dir = NA, lrt_p = pchisq(max(2 * (f1$ll - f0$ll), 0), 1, lower.tail = FALSE))
      }
      nest(cbind(L = L), "proxy_given_ligand"); nest(cbind(R = R), "proxy_given_receptor")
      nest(cbind(L = L, R = R), "proxy_given_both_ends")
    }
    # 共线性与基因重叠(E7)
    pg <- read.delim(cfg$parent_genes, stringsAsFactors = FALSE)
    shared <- intersect(pg$gene[pg$signature == a$parent], c(Lg, Rg))
    r <- cor(par, px, use = "complete.obs")
    ov[[length(ov) + 1]] <- data.frame(axis = a$axis_id, parent = a$parent, r_parent_proxy = r, VIF = 1 / (1 - r^2),
      n_shared_genes = length(shared), shared = paste(shared, collapse = ","))
    # 患者层 bootstrap: 父程序系数变化
    y <- make_y(cl, cfg$primary_endpoint)
    RNGkind("L'Ecuyer-CMRG"); set.seed(cfg$seed + 7 * p)
    bs <- unlist(mclapply(seq_len(cfg$B_boot), function(b) {
      i <- unlist(lapply(split(seq_along(st), st), function(ix) ix[sample.int(length(ix), replace = TRUE)]))
      m <- tryCatch(model_ac(par[i], px[i], age[i], list(t = y$t[i], e = y$e[i], ok = y$ok[i]), st[i]), error = function(e) NULL)
      if (is.null(m)) NA else m$C$beta[["parent"]] - m$A$beta[["parent"]]
    }, mc.cores = cfg$cores))
    obs <- model_ac(par, px, age, y, st)
    boot[[length(boot) + 1]] <- data.frame(axis = a$axis_id, endpoint = cfg$primary_endpoint,
      d_parent = obs$C$beta[["parent"]] - obs$A$beta[["parent"]],
      lo = quantile(bs, .025, na.rm = TRUE), hi = quantile(bs, .975, na.rm = TRUE),
      HR_parent_A = exp(obs$A$beta[["parent"]]), HR_parent_C = exp(obs$C$beta[["parent"]]), t(setNames(hr(obs$C, "parent"), c("C_HR", "C_lo", "C_hi", "C_p"))),
      t(setNames(hr(obs$C, "proxy"), c("px_HR", "px_lo", "px_hi", "px_p"))))
  }
  rb <- function(l) { n <- unique(unlist(lapply(l, names))); do.call(rbind, lapply(l, function(x) { x[setdiff(n, names(x))] <- NA; x[n] })) }
  write_new(do.call(rbind, rob), "robust_covariates.tsv")
  write_new(do.call(rbind, loco), "robust_leave_one_cancer.tsv")
  write_new(do.call(rbind, het), "robust_per_cancer.tsv")
  write_new(rb(ends), "robust_meta_and_ends.tsv")
  write_new(do.call(rbind, ov), "robust_overlap_vif.tsv")
  write_new(do.call(rbind, boot), "robust_parent_shift_bootstrap.tsv")
}

## ---------- step: decide ----------
step_decide <- function(d) {
  gate()
  ns <- read.delim(out("null_summary.tsv")); rob <- read.delim(out("robust_covariates.tsv"))
  lo <- read.delim(out("robust_leave_one_cancer.tsv")); me <- read.delim(out("robust_meta_and_ends.tsv"))
  bo <- read.delim(out("robust_parent_shift_bootstrap.tsv")); fam <- read.delim(out("family_all_axes.tsv"))
  ep <- cfg$primary_endpoint; al <- cfg$alpha
  prim <- d$axes[d$axes$family == "primary", ]
  rows <- lapply(seq_len(nrow(prim)), function(p) {
    ax <- prim$axis_id[p]; adv <- prim$direction[p] == "adverse"
    g <- function(f) ns$p_emp[ns$axis == ax & ns$family == f]
    agree <- function(lo_, hi_) if (adv) lo_ > 1 else hi_ < 1
    s4 <- rob[rob$axis == ax & rob$endpoint == ep & rob$model == "S4_full", ]
    s0 <- rob[rob$axis == ax & rob$endpoint == ep & rob$model == "S0_age", ]
    retained <- log(s4$HR) / log(s0$HR)
    ng <- me[me$axis == ax & me$endpoint == ep & me$term %in% c("proxy_given_ligand", "proxy_given_receptor"), ]
    l <- lo[lo$axis == ax & lo$endpoint == ep, ]
    m <- me[me$axis == ax & me$endpoint == ep & me$term == "RE_per_cancer", ]
    data.frame(axis = ax, p_N2a = g("N2a"), p_N2b = g("N2b"), p_N3L = g("N3L"), p_N3R = g("N3R"),
      p_gene = g("N2a"), p_axis_IUT = max(g("N2a"), g("N3L"), g("N3R")),
      S4_HR = s4$HR, S4_lo = s4$lo, S4_hi = s4$hi, S4_retained = retained,
      S4_pass = agree(s4$lo, s4$hi) && retained >= cfg$min_retained,
      nested_pass = all(mapply(agree, ng$lo, ng$hi)),
      LOCO_all_same = all(if (adv) l$HR > 1 else l$HR < 1), LOCO_max_p = max(l$lrt_p),
      RE_HR = m$HR, RE_lo = m$lo, RE_hi = m$hi, I2 = m$I2, RE_pass = agree(m$lo, m$hi))
  })
  r <- do.call(rbind, rows)
  r$p_gene_holm <- p.adjust(r$p_gene, "holm"); r$p_axis_holm <- p.adjust(r$p_axis_IUT, "holm")
  r$robust <- r$S4_pass & r$LOCO_all_same & r$RE_pass
  r$tier <- with(r, ifelse(p_axis_holm < al & nested_pass & robust, "A_axis_robust",
              ifelse(p_axis_holm < al & nested_pass, "A_minus_axis_confounding_sensitive",
              ifelse(p_gene_holm < al & robust, "B_gene_level_robust",
              ifelse(p_gene_holm < al, "B_minus_gene_level_confounding_sensitive", "C_not_above_null")))))
  # 关系分类(解耦/吸收/归入)
  fr <- fam[fam$form == "min" & fam$endpoint == ep, ]
  cls <- do.call(rbind, lapply(d$axes$axis_id[d$axes$family %in% c("primary", "control")], function(ax) {
    b <- bo[bo$axis == ax, ]; f <- fr[fr$axis == ax, ]
    shift <- b$lo > 0 | b$hi < 0
    marginal_signal <- f$p_proxy_B < al && f$HR_proxy_B > 1 && b$HR_parent_A > 1
    label <- if (!(f$p_proxy_B < al) && !(f$lrt_p < al)) "no_signal"
      else if (marginal_signal && f$lo >= cfg$equiv_lo && f$hi <= cfg$equiv_hi) "subsumed_equivalence_met"
      else if (!(f$lrt_p < al)) "indeterminate_no_increment_equivalence_not_met"
      else if (shift && b$C_hi < 1 && f$HR_proxy > 1) "decoupled"
      else if (shift && b$C_lo <= 1 && b$C_hi >= 1 && f$HR_proxy > 1) "absorbed"
      else "increment_other"
    data.frame(axis = ax, HR_proxy_B = f$HR_proxy_B, p_proxy_B = f$p_proxy_B, HR_proxy = f$HR_proxy, lo = f$lo, hi = f$hi, lrt_p = f$lrt_p, HR_parent_A = b$HR_parent_A,
               HR_parent_C = b$HR_parent_C, parent_C_lo = b$C_lo, parent_C_hi = b$C_hi,
               parent_shift_lo = b$lo, parent_shift_hi = b$hi, class = label)
  }))
  # 特异性: 全家族中有多少轴通过
  pr <- fr[fr$tier %in% c("priority", "focus", "background"), ]
  spec <- data.frame(n_axes = nrow(pr), frac_q05_adverse = mean(pr$lrt_q < .05 & pr$HR_proxy > 1),
                     frac_q05_any = mean(pr$lrt_q < .05))
  spec$specificity_flag <- spec$frac_q05_any > cfg$max_family_pass
  write_new(r, "decision_primary.tsv"); write_new(cls, "decision_relation.tsv"); write_new(spec, "decision_family_specificity.tsv")
  print(r[, c("axis", "p_N2a", "p_N2b", "p_N3L", "p_N3R", "p_axis_holm", "nested_pass", "S4_retained", "S4_pass", "LOCO_all_same", "RE_pass", "tier")])
  print(cls[, c("axis", "HR_proxy", "lrt_p", "HR_parent_A", "HR_parent_C", "class")]); print(spec)
}

## ---------- 主流程 ----------
receipt <- function(st) {
  cat(sprintf("step=%s\ttime=%s\tscript_sha256=%s\tR=%s\tsurvival=%s\n", st, format(Sys.time(), "%FT%T%z"),
      sha256(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), getRversion(),
      packageVersion("survival")), file = out("receipt.tsv"), append = TRUE)
}
d <- load_all()
steps <- list(repro = step_repro, marginal = step_marginal, family = step_family, null = step_null,
              robust = step_robust, decide = step_decide)
todo <- if (step == "all") names(steps) else step
for (st in todo) { log_msg("== step", st); steps[[st]](d); receipt(st) }
log_msg("done")
