# 合成测试配置: 参数结构与正式配置相同, 仅 B 较小以便本地快速运行
syn <- Sys.getenv("E2_SYN")
expr <- file.path(syn, "expr.rds"); clin <- file.path(syn, "clin.tsv"); parents <- file.path(syn, "parents.tsv")
parent_genes <- file.path(syn, "parent_genes.tsv"); axes <- file.path(syn, Sys.getenv("E2_AXES", "axes.tsv"))
lr_resource <- file.path(syn, "lr_resource.tsv"); covariates <- file.path(syn, "covariates.tsv")
published <- file.path(syn, "published.tsv"); out_dir <- Sys.getenv("E2_OUT")
primary_endpoint <- "OS"; endpoints <- c("OS", "PFI"); proxy_forms <- c("min", "mean")
repro_tol <- 0.01; alpha <- 0.05; equiv_lo <- 0.90; equiv_hi <- 1.10; max_family_pass <- 0.5
expr_floor <- 0; min_detect_frac <- 0.2
K_expr <- 200; K_prog <- 50; B <- as.integer(Sys.getenv("E2_B", "2000")); B_boot <- 100; seed <- 20261106; cores <- 4
min_cancer_n <- 30; min_cancer_events <- 10
stage_col <- "stage"; sex_col <- "sex"; estimate_cols <- c("estimate_stromal", "estimate_immune")
fraction_cols <- c("macrophage_fraction", "caf_fraction", "cd8_fraction"); purity_col <- "purity_absolute"
desmoplasia_col <- "desmoplasia_myCAF"; alt_parent_map <- list(CM_T = "CM_T_V7")
min_retained <- 0.5
# 测试专用: 埋入的强基因不进零分布候选池; 候选池较小, 否则会抬高 p 值下限。正式配置为 character()
null_exclude <- c("LIG_TRUE","REC_TRUE1","REC_TRUE2","LIG_ONLY","REC_NOISE","LIG_NULL","REC_NULL","PAR1","PAR2","PAR3",
  "DESMO1","DESMO2","DESMO3","SUB_L","SUB_R","CONF_L","CONF_R","ABS_L","ABS_R","PABS1","PABS2","PABS3")
