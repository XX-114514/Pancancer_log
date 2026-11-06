# 合成测试数据: 埋入已知答案, 用来验证 E2 脚本的判定是否正确
# 用法: Rscript make_synthetic.R <outdir>
args <- commandArgs(TRUE); o <- args[1]; dir.create(o, recursive = TRUE, showWarnings = FALSE)
set.seed(20261106)
n <- 9000; ng <- 1500; K <- 10
cancer <- sample(paste0("C", seq_len(K)), n, TRUE)
genes <- paste0("G", seq_len(ng))
fac <- matrix(rnorm(n * 8), n)
load <- matrix(rnorm(ng * 8, sd = .35), ng)
base <- runif(ng, 1, 8)
E <- t(fac %*% t(load)) + matrix(rnorm(n * ng), ng) + base
rownames(E) <- genes; colnames(E) <- paste0("S", seq_len(n))
# 指定基因
named <- c("LIG_TRUE","REC_TRUE1","REC_TRUE2","LIG_ONLY","REC_NOISE","LIG_NULL","REC_NULL",
           "PAR1","PAR2","PAR3","DESMO1","DESMO2","DESMO3","SUB_L","SUB_R","CONF_L","CONF_R","ABS_L","ABS_R","PABS1","PABS2","PABS3",paste0("PROG",1:60))
rownames(E)[seq_along(named)] <- named; genes <- rownames(E)
zc <- function(x) ave(x, cancer, FUN = function(v) (v - mean(v)) / sd(v))
L <- rnorm(n); R1 <- .5*L + rnorm(n); R2 <- .5*L + rnorm(n); RO <- rnorm(n); LO <- rnorm(n)
E["LIG_TRUE",] <- 5 + L; E["REC_TRUE1",] <- 5 + R1; E["REC_TRUE2",] <- 4 + R2
E["LIG_ONLY",] <- 5 + LO; E["REC_NOISE",] <- 5 + RO
Rs <- (zc(R1) + zc(R2)) / 2
trueAxis <- pmin(zc(L), Rs)
# 父程序: 含保护效应, 与真轴正相关(用于检验解耦)
PL <- rnorm(n)
par <- .6 * trueAxis + PL
for (g in c("PAR1","PAR2","PAR3")) E[g,] <- 5 + par + rnorm(n, sd = .4)
# 促纤维化程序与归入轴: 归入轴只是该程序的带噪版本
DS <- rnorm(n); for (g in paste0("DESMO",1:3)) E[g,] <- 5 + DS + rnorm(n, sd = .4)
E["SUB_L",] <- 5 + DS + rnorm(n, sd = .9); E["SUB_R",] <- 5 + DS + rnorm(n, sd = .9)
# 混杂轴: 两端都只是某细胞比例的读数, 风险来自该比例本身
MF <- rnorm(n); E["CONF_L",] <- 5 + MF + rnorm(n, sd = .7); E["CONF_R",] <- 5 + MF + rnorm(n, sd = .7)
# 吸收轴: 父程序没有独立效应, 它的不良关联全部来自该轴
AL <- rnorm(n); AR <- .5*AL + rnorm(n); E["ABS_L",] <- 5 + AL; E["ABS_R",] <- 5 + AR
absAxis <- pmin(zc(AL), zc(AR)); PAB <- .8*absAxis + rnorm(n)
for (g in c("PABS1","PABS2","PABS3")) E[g,] <- 5 + PAB + rnorm(n, sd = .3)
# 一批有预后的普通基因, 让匹配池有预后强度可匹配
pg <- paste0("PROG",1:60); pz <- matrix(rnorm(n*60), n); colnames(pz) <- pg
for (g in pg) E[g,] <- 5 + pz[,g]
# 风险
lp <- .35*trueAxis - .30*zc(par) + .30*zc(LO) + .30*zc(DS) + .35*MF + .30*absAxis + as.vector(pz %*% rep(c(.12,-.12), 30))
age <- rnorm(n, 60, 10); lp <- lp + .01*(age-60)
h0 <- runif(K, .02, .08)[as.integer(factor(cancer))]
t_ev <- rexp(n, h0 * exp(lp)); t_c <- runif(n, 1, 40)
os <- as.integer(t_ev <= t_c); ost <- pmin(t_ev, t_c)
t_ev2 <- rexp(n, h0 * 1.5 * exp(lp)); pfi <- as.integer(t_ev2 <= t_c); pfit <- pmin(t_ev2, t_c)
S <- colnames(E)
saveRDS(E, file.path(o, "expr.rds"))
w <- function(x, f) write.table(x, file.path(o, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(data.frame(sample = S, patient = paste0("P", seq_len(n)), cancer = cancer, age = age, OS = os, OS.time = ost, PFI = pfi, PFI.time = pfit), "clin.tsv")
sc <- function(g) rowMeans(sapply(g, function(x) zc(E[x,])))
w(data.frame(sample = S, CM_T = sc(c("PAR1","PAR2","PAR3")), CM_T_V7 = sc(c("PAR1","PAR2")), desmoplasia_myCAF = sc(paste0("DESMO",1:3)), CM_ABS = sc(c("PABS1","PABS2","PABS3"))), "parents.tsv")
w(rbind(data.frame(signature="CM_T", gene=c("PAR1","PAR2","PAR3")), data.frame(signature="desmoplasia_myCAF", gene=paste0("DESMO",1:3)), data.frame(signature="CM_ABS", gene=c("PABS1","PABS2","PABS3"))), "parent_genes.tsv")
ax <- data.frame(
  axis_id = c("AX_TRUE","AX_LIGONLY","AX_NULL","AX_CONF","AX_SUBSUMED","AX_ABSORB"),
  family = c("primary","primary","primary","primary","control","control"),
  parent = c("CM_T","CM_T","CM_T","CM_T","desmoplasia_myCAF","CM_ABS"),
  ligand_genes = c("LIG_TRUE","LIG_ONLY","LIG_NULL","CONF_L","SUB_L","ABS_L"),
  receptor_genes = c("REC_TRUE1,REC_TRUE2","REC_NOISE","REC_NULL","CONF_R","SUB_R","ABS_R"),
  direction = c("adverse","adverse","adverse","adverse","none","none"), stringsAsFactors = FALSE)
set.seed(1); pool <- setdiff(genes, named)
extra <- data.frame(axis_id = paste0("BG", 1:40), family = "priority", parent = "CM_T",
  ligand_genes = sample(pool, 40), receptor_genes = sample(pool, 40), direction = "none", stringsAsFactors = FALSE)
w(rbind(ax, extra), "axes.tsv")
lrp <- data.frame(ligand_genes = c(ax$ligand_genes, sample(pool, 30)), receptor_genes = c(ax$receptor_genes, sample(pool, 30)))
w(lrp, "lr_resource.tsv")
# 协变量: 真轴是独立生物信号, LIG_ONLY 也独立; 另设一个能解释 AX_CONF 的场景不在此测试
cv <- data.frame(sample = S, sex = sample(c("F","M"), n, TRUE), stage = sample(c("I","II","III","IV",NA), n, TRUE),
  estimate_stromal = rnorm(n), estimate_immune = rnorm(n), macrophage_fraction = MF + rnorm(n, sd = .2), caf_fraction = rnorm(n), cd8_fraction = rnorm(n), purity_absolute = ifelse(runif(n) < .3, NA, rnorm(n)))
w(cv, "covariates.tsv")
# “已发表值”: 由同样模型算出, 用于测试 repro 门
library(survival)
Rz <- rowMeans(sapply(c("REC_TRUE1","REC_TRUE2"), function(g) zc(E[g,])))
px <- pmin(zc(E["LIG_TRUE",]), Rz)
pz2 <- (age - mean(age)) / sd(age)
pp <- read.delim(file.path(o, "parents.tsv"))$CM_T
m <- coxph(Surv(ost, os) ~ pp + px + pz2 + strata(cancer))
w(data.frame(axis = "AX_TRUE", proxy = "min_AX_TRUE", endpoint = "OS", HR_proxy_in_C = exp(coef(m)[["px"]]), HR_parent_in_C = exp(coef(m)[["pp"]]), n = n), "published.tsv")
# 校准: 30 条完全随机的轴当作主轴, 考察假阳性率
set.seed(2); cal <- data.frame(axis_id = paste0("CAL", 1:30), family = "primary", parent = "CM_T",
  ligand_genes = sample(pool, 30), receptor_genes = sample(pool, 30), direction = "adverse", stringsAsFactors = FALSE)
w(cal, "axes_calibration.tsv")
cat("synthetic written to", o, "\n")
