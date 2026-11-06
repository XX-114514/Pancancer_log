# E2 正式配置(第二阶段冻结)。路径用逻辑根目录, 不写真实服务器路径到公开仓库。
# 冻结后不得修改; 任何改动都要升版本并记录原因。
root <- Sys.getenv("E2_ROOT")                 # 服务器上的输入目录
expr <- file.path(root, "expr.rds")  # 基因×样本, log2(TPM+0.001), 行名为基因符号且唯一
clin <- file.path(root, "clin.tsv")           # sample, patient, cancer, age, OS, OS.time, PFI, PFI.time
parents <- file.path(root, "parents.tsv")     # sample + 各父程序打分列(直接取原输出)
parent_genes <- file.path(root, "parent_genes.tsv")  # signature, gene
axes <- file.path(root, "axes.tsv")           # axis_id, family, parent, ligand_genes, receptor_genes, direction
                                              # family: primary | control | focus | priority | background
                                              # primary: SPP1_macrophage, LGALS1_myCAF_CD8 (direction=adverse)
                                              # control: collagen_DDR1, COL1A1_only_DDR1 (direction=none)
                                              # axis_id 必须与 published 的 axis 列一致
lr_resource <- file.path(root, "lr_resource.tsv")    # ligand_genes, receptor_genes (LIANA consensus)
covariates <- file.path(root, "covariates.tsv")      # sample + 下方各协变量列
published <- file.path(root, "published.tsv")  # 公开表原样复制
out_dir <- Sys.getenv("E2_OUT")               # 必须是新目录

## 以下为冻结参数
primary_endpoint <- "OS"; endpoints <- c("OS", "PFI"); proxy_forms <- c("min", "mean")
repro_tol <- 0.01; alpha <- 0.05
min_retained <- 0.5                           # S4 至少保留 S0 的 log HR 比例
equiv_lo <- 0.90; equiv_hi <- 1.10            # 归入等效界
max_family_pass <- 0.5                        # 特异性上限
expr_floor <- 0; min_detect_frac <- 0.2       # log2(TPM+0.001)>0 约等于 TPM≥1; 至少 20% 样本检出
K_expr <- 500; K_prog <- 100; B <- 5000L; B_boot <- 2000L; seed <- 20261106L
min_cancer_n <- 50; min_cancer_events <- 20
null_exclude <- character()                   # 正式运行不额外排除

## 运行参数(不影响结果语义)
cores <- 8L

## 协变量列名已按实际表填写
stage_col <- "stage"; sex_col <- "sex"
estimate_cols <- c("estimate_stromal", "estimate_immune")
fraction_cols <- c("leukocyte_fraction", "macrophage_fraction", "caf_fraction", "cd8_fraction")
purity_col <- "purity_absolute"
desmoplasia_col <- "desmoplasia_myCAF"
alt_parent_map <- NULL  # V7 bulk 打分尚未构建: S7 记不可估计, 不替代
