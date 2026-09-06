# CM 模块 → 具体 ligand–receptor 通讯 → 空间支持 → 受体端功能 → 临床关联：完整证据链审计

生成日期：2026-09-06。工作目录 `communication_evidence_chain_20260906/`（全部为新输出，**未覆盖**任何 Phase B 结果）。
计算：本地 V100 机；复用已完成的逐样本 LIANA、Visium 空间置换、TCGA toil、ICI bulk 队列，未重跑数百万细胞。

---

## 0. 重要：CM 编号对齐（用户 brief 用的是旧编号）

用户描述里的 CM03 / CM06 / CM09 是 **v5 之前**的 NMF 编号。当前冻结契约
（`final_annotation_v6_identity_v5_malignancy_v3`，20260825 版 CoVarNet K=9）编号不同。
按 `frozen/reports/V5_CM_IMMUNE_COMMUNICATION_COMPARISON_CN.md` 的旧→新对齐表（样本活动相关 65% + top-10 谱系组成余弦 35%）：

| 用户说的 | = 当前 CM | 对齐强度 | 模块本质 |
|---|---|---|---|
| **CM03**（myCAF→malignant，COL1A1–DDR1）| **CM07** | strong (0.86) | 促纤维化 / myCAF–恶性上皮结构性通讯 |
| **CM06**（myCAF→T / myeloid）| **CM08** | strong (0.84) | 炎症单核/巨噬 + myCAF→CD8 抑制界面 |
| **CM09**（FCN1 monocyte→SPP1 macrophage→cytotoxic T/NK）| **CM01** | strong (0.76) | 血管周 SPP1-TAM / 髓系–内皮–淋巴界面 |

本报告一律用**当前编号**。注意：**当前 CM03** 是"严格恶性分数"冷肿瘤模块（strict malignant ρ=0.70，
immune-hot proxy ρ=−0.45），top-10 里**没有 myCAF**，且"没有通过研究级过滤的非背景跨细胞关系"——
用户记忆中的"CM03 myCAF→malignant"确实存在，但它在当前契约里叫 **CM07**。

---

## 1. 审计：现有 consensus 通讯是怎么从单样本 LIANA 汇总出来的

冻结 LIANA run：`Project_v3/runs/.../20260825_v5_final_v1/`。

### 1.1 逐样本 LIANA（脚本 07）
- 每个**原始样本独立**跑 `liana.mt.rank_aggregate`（RRA 聚合，resource=`consensus`）。
- 表达输入：counts → 样本内 `normalize_total(1e4)` → `log1p`。
- 分组：`final_annotation_with_malignancy_v6`（细分细胞状态），每组 **≥20 细胞**才纳入。
- `expr_prop=0.10`（配体在发送组 ≥10% 细胞表达、受体在接收组 ≥10% 细胞表达才计），`seed=1337`。
- **within-CM 过滤**：只保留 source 和 target **同时**属于某个 CM 的 CoVarNet top-10 成员的互作，
  打上 `shared_cm` 标签。
- 结果：1322 样本，1301 成功，21 因 <2 组跳过。聚合表 `all_within_cm_interactions.tsv.gz` 共 **1121 万行**
  （每样本 × 每有序 (source,target) × 每 LR × 每方向）。

### 1.2 consensus 汇总（脚本 08 `build_communication`）
1. **事件过滤**：逐样本逐 LR 逐方向，保留 `magnitude_rank ≤ 0.05 且 specificity_rank ≤ 0.05` 的行
   → 38.2 万 "rank05 事件"（`communication_rank05_events.tsv.gz`）。
2. 按 `shared_cm` 展开成单 `cm`。
3. 按 **(cm, source, target, ligand_complex, receptor_complex)** 分组，统计：
   - `n_samples` = 去重样本数；`n_gse` = 去重 GSE 数；`n_cancer_types` = 去重癌种数
   - `mean_magnitude_rank`、`mean_specificity_rank`、`mean_lr_means`、`mean_lrscore`（对**通过过滤的样本**取均值）
   - `replicate_rank_score = log1p(n_samples) × −log10( sqrt(mean_mag_rank × mean_spec_rank) )`
   → **29 605 条 consensus**（`communication_consensus_global.tsv.gz`）。
4. `communication_research_priorities.tsv`（139 条）再加：**可检验分母** `n_eligible_samples`
   （source 和 target 两组都 ≥20 细胞的样本数）、`support_fraction = n_samples / n_eligible_samples`
   及其 Wilson 95% CI、`mechanism_class`、`research_priority_score`；并剔除自通讯、复核端点、
   B2M–KLRD/KLRC 高频背景轴，要求 ≥10 可检验样本、≥3 GSE、≥3 癌种。

**口径限制（原作者已声明，本报告沿用）**：consensus 的 n_samples/n_gse/n_cancer 是**复现广度**证据，
**不是**因果通讯；比较单位必须是样本，不能用 pooled-cell P 值替代生物学重复。

---

## 2. 重点 CM 的 communication-level 复分析（sample-level，discovery 复用）

方法：从 `all_within_cm_interactions.tsv.gz` 流式抽出每条目标 axis 的**全部逐样本行**（不只 rank≤0.05 的），
建立 `focus_axes_sample_level.tsv.gz`；分母 = source 与 target 两个细胞状态都在该样本出现过的样本集。
三档阈值：strict (rank≤0.01)、standard (≤0.05)、loose (≤0.10)。

### 2.1 结果总表（standard 阈值）

| axis | CM | n_samp | n_GSE | n_cancer | support | 最大单 GSE 占比 | strict/standard 稳健比 | spec_rank 中位 |
|---|---|---|---|---|---|---|---|---|
| **COL1A1→DDR1  myCAF→恶性上皮** | CM07 | 115 | 28 | 22 | 0.57 | 0.12 | **0.57** | 6e-6 |
| COL1A1→DDR2  myCAF autocrine | CM07 | 238 | 34 | 28 | 0.62 | 0.10 | 0.48 | 4e-7 |
| COL1A1→CD44  myCAF→单核 | CM07 | 145 | 34 | 27 | 0.92 | 0.09 | **0.91** | 3e-4 |
| **LGALS1→PTPRC  myCAF→CD8T** | CM08 | 291 | 41 | 33 | 0.92 | 0.10 | **0.42** | 1e-2 |
| MYL9→CD69  myCAF→CD8T | CM08 | 307 | 39 | 31 | 0.97 | 0.10 | 0.79 | 3e-4 |
| COL1A1→CD44  myCAF→CD8T | CM08 | 265 | 39 | 31 | 0.84 | 0.10 | 0.65 | 4e-3 |
| S100A9→CD68  单核→巨噬 | CM08 | 317 | 41 | 34 | 0.99 | 0.11 | **0.92** | 1e-4 |
| S100A8→CD69  单核→CD8T | CM08 | 324 | 41 | 34 | 0.95 | 0.13 | 0.87 | 6e-4 |
| **S100A9→CD68  FCN1单核→SPP1巨噬** | CM01 | 66 | 23 | 22 | 1.00 | 0.15 | **0.95** | 5e-5 |
| C1QA→CD93  SPP1巨噬→内皮 | CM01 | 121 | 30 | 25 | 1.00 | 0.15 | **0.99** | 8e-7 |
| APP→CD74  内皮→SPP1巨噬 | CM01 | 121 | 30 | 25 | 1.00 | 0.15 | 0.79 | 5e-3 |

（完整逐样本表、per-GSE / per-cancer 拆分见 `tables/focus_axes_*.tsv*`）

### 2.2 关键判断

- **COL1A1–DDR1（用户主假设）**：确实**广泛复现**（115 样本 / 28 GSE / 22 癌种，乳腺/肺/胃/CRC/卵巢/
  食管/膀胱/HNSCC/间皮/BCC/panNET 都有），**不是**单个超大队列驱动（最大 GSE 仅占 12%）；
  specificity 极强（中位 6e-6，一旦检出就高度细胞类型特异）。**但**：(a) magnitude_rank 中位仅 ~0.01
  （胶原表达高、肿瘤 DDR1 表达偏低，处于 magnitude 边界）；(b) **strict 阈值下只有 57% 的样本保留**——
  是所有 CM07 axis 里稳健性最差的一条。相比之下 COL1A1→CD44（myCAF→单核）稳健比 0.91、support 0.92。
- **LIANA consensus resource 里只有 COL1A1–DDR1，没有 COL1A2–DDR1 / COL3A1–DDR1**（collagen-DDR 扫描
  `collagen_ddr_sweep_summary.tsv` 里 DDR1 只与 COL1A1 配对）。所以"collagen–DDR1"在单细胞层面
  实际上=**COL1A1–DDR1**。
- **FCN1→SPP1（CM01）**：66/66 可检验样本全部支持，23 GSE，稳健比 0.95——极干净，但可检验分母小
  （FCN1 单核 + SPP1 巨噬同时 ≥20 细胞的样本不多）。下游 SPP1巨噬→内皮 C1QA–CD93 121/121、稳健比 0.99。
- **myCAF→CD8 抑制轴（CM08）**：LGALS1–PTPRC / MYL9–CD69 / COL1A1–CD44 都高 support，但
  **LGALS1–PTPRC 在 strict 阈值下只剩 42%**（PTPRC 是泛免疫标记，特异性弱）；MYL9–CD69 更稳（0.79）。
- **纯组成驱动的**：S100A8/A9–CD68（单核→巨噬，support 0.99）——本质就是"髓系细胞多"。

---

## 3. 空间支持

### 3.1 Visium 跨切片置换（复用，5 数据集 61 切片：ccRCC×2、乳腺、HCC、HNSCC）

`spatial_expansion_20260903/runs/20260904_all_visium_cm_region_communication_v1/`。
置换 z>0 = 配体-发送态活动与受体-接收态活动在空间邻域**共富集**。

| axis | 5 数据集里 z>0 的数 | 中位 z | 判断 |
|---|---|---|---|
| VCAN→ITGB1  单核↔myCAF (CM07) | 5/5 | 4.4–10.5 | **强、泛癌一致** |
| S100A9→CD68  单核→巨噬 (CM08) | 5/5 | 1.7–4.2 | **一致** |
| C1QA→CD93  SPP1巨噬→内皮 (CM01) | 4/5 | 2.6 | 基本一致 |
| S100A8→CD69  单核→CD8T (CM08) | 4/5 | 0.5 | 弱 |
| COL1A1→CD44  myCAF→单核 (CM07) | 2/5（HCC + 1 ccRCC 阳，乳腺/HNSCC/另一 ccRCC 阴）| −0.5 | **混杂** |
| LGALS1→PTPRC  myCAF→CD8T (CM08) | 2/5（HCC + 乳腺阳，2 ccRCC + HNSCC 阴）| −1.4 | **混杂——与"排斥"相符（抑制 = 不共定位）** |
| APP→CD74  内皮→SPP1巨噬 (CM01) | 2/5 | −5.1 | 背景/方向不稳 |

**注意：Visium run 没有测 DDR1 轴**——COL1A1–DDR1 的空间证据只有下面的 CosMx。

### 3.2 ccRCC CosMx（Manley lab，新分析，3 组织 199,112 细胞，978 基因，FOV 坐标）

`scripts/02b_spatial_cosmx_col_ddr1.py`；k=12 邻居；999 次标签置换。
胶原源 = Fibroblast+Myofibroblast（COL1A1 CPM 298 / 62）；DDR1 受体端 = Tumor 区室的上皮性细胞；
CD8 = CD8.T.cell + NKT。

| 检验 | 结果（RCC3 / RCC4 / RCC5） | 解读 |
|---|---|---|
| **T1 直接毗邻**：DDR1⁺肿瘤细胞 vs DDR1⁻，胶原高 CAF 邻居数 | obs_diff **−0.013 / −0.21 / −0.013**，perm p 0.004 / 0.001 / 0.042 | **DDR1⁺肿瘤细胞邻居里胶原高 CAF 更少**——ccRCC 里**不是**近距 juxtacrine（见下方 caveat）|
| **T3 CD8 排斥**：肿瘤邻域 CD8 数 ~ 局部胶原暴露 + 自身 DDR1 | 局部胶原 β −0.06/−0.13/−0.16（p 1e-12～1e-90）；DDR1 β −0.07/−0.07/−0.10（p 1e-15～1e-30）；交互项不显著 | **胶原富集邻域 与 DDR1 高肿瘤细胞 都独立地 CD8 更少**（加性，非乘性）|
| **T3b CD8 定位**：CD8 到"interface(胶原×DDR1)高"肿瘤 vs "低"肿瘤的距离 | 中位 320 vs 76 μm（RCC3）、556 vs 100 μm（RCC5），3/3 p<1e-30 | **CD8 被空间排斥在胶原×DDR1 高的肿瘤区域之外** |
| **T4 受体端后果**：肿瘤 EMT/MMP 程序 ~ 局部胶原 + 自身 DDR1 | 胶原 β +0.04（p 1e-35～1e-143）；DDR1 β +0.02（p<1e-8），3/3 | **邻近胶原 / 自身 DDR1 高的肿瘤细胞更偏间叶/MMP 表型**——与 DDR1 已知生物学一致 |

**ccRCC caveat（关键）**：ccRCC 恶性细胞保留近端小管身份，正常近端小管本身 DDR1⁺（Proximal.tubule CPM 20，
Epithelial.progenitor 21）——这正是记忆里记录的"CM07 在 ccRCC 有上皮污染"。所以 **T1 的阴性不能推广**：
ccRCC 是检验肿瘤 DDR1 轴最差的组织。可迁移的结论是 **T3/T3b/T4**——
"胶原 + DDR1 高的肿瘤区域 = CD8 排斥 + 间叶化"，方向在 3/3 组织一致且极显著。

---

## 4. Patient-level bulk proxy 模型（Model A/B/C）

**严格 caveat（贯穿全节）**：bulk RNA **无法**恢复 source→target。配体+受体同高只是 **proxy**，
**不是 LIANA 通讯验证**。这里只回答一个问题：**LR proxy 是否在 parent CM/signature 之外提供额外信息**。

proxy 定义：队列内（TCGA 为癌种内）z-score 后，`L_score`=配体基因均值 z、`R_score`=受体基因均值 z、
`min`=两者取小（"通讯需两端都高"）、`mean`=两者均值、`prod`=乘积。
- Model A: `outcome ~ parent_signature (+age, TCGA 按癌种分层)`
- Model B: `outcome ~ LR_proxy`
- Model C: `outcome ~ parent_signature + LR_proxy` → 报告 proxy 偏系数 + LRT(C vs A)
- Model D: `~ parent + L_score + R_score` → 信号是配体还是受体端

### 4.1 TCGA 泛癌 OS（分层 Cox，9034 例 / 2664 事件；`bulk_proxy_TCGA_models.tsv`）

| axis | parent | proxy(min) 在 C 里 HR (p) | parent 在 C 里 HR (p) | LRT C vs A | 结论 |
|---|---|---|---|---|---|
| **SPP1 macrophage** (CM01) | CM01 sig | **1.23 (5e-21)** | 0.90 (6e-7, 转为保护) | **1.8e-21** | proxy 远强于 parent，独立不良预后 |
| **LGALS1 myCAF→CD8** (CM08) | CM08 sig | **1.10 (8e-4)**；mean 1.15 (2e-5) | 0.99 (0.64，parent 归零) | **1.7e-5** | galectin-1 proxy 吸收了 CM08 的预后信号 |
| collagen–DDR1 (min) | desmoplasia_myCAF | 0.95 (0.028，条件后**翻为保护**) | 1.22 (2e-17) | 0.030 | 与 parent 共线，**符号翻转 = 非独立** |
| COL1A1–DDR1 (min) | desmoplasia_myCAF | 0.96 (0.062，ns) | 1.21 (9e-17) | 0.064 | **不在 desmoplasia 之外提供稳健信息** |

- Model D：SPP1 轴 p_L(SPP1)=9e-9、p_R(CD44/ITGAV/ITGB1)=1e-10（两端都贡献）；
  LGALS1 轴 p_L(LGALS1)=8e-12 >> p_R(PTPRC/CD69)=2e-5（**配体端 galectin-1 主导**）。

### 4.2 TCGA 受体端后果：`CYT ~ L × R`（每癌种 OLS，`bulk_proxy_TCGA_receiver_consequence.tsv`）

- **DDR1（受体）本身是稳健的泛癌"免疫冷"标记**：beta_R 在多数癌种显著为负——
  LUAD −0.22 (p 5e-8)、LUSC −0.30 (p 1e-14)、STAD −0.22 (p 2e-7)、SARC −0.27、ESCA −0.29、PAAD −0.34；
  IMvigor210 里也是 −0.35 (p 1e-13)。胶原（配体）beta_L 多数为正（胶原多 = 基质多 = 免疫细胞也多）。
- **胶原 × DDR1 的协同免疫排斥交互项（beta_LxR < 0，BH q<0.1）**只在一个子集显著：
  **THCA、LUAD、ESCA、HNSC、THYM、GBM、SKCM(COL1A1-only)、PAAD(名义)**——
  即**强促纤维化癌 + 部分上皮癌**，**不含 ccRCC/KIRC**（DDR1 在肾被正常上皮身份混杂），
  也不含 PRAD/OV（那里交互为正）。

### 4.3 ICI 队列（response = logistic + LOO-AUC + 置换；OS/PFS = Cox）

| axis / proxy | Riaz2017 黑素瘤 (n=71) | Rose2021 尿路 (n=88 OS) | IMvigor210 尿路 (n=347) |
|---|---|---|---|
| **SPP1 macrophage** min | 应答 OR 0.52 p0.04（NR↑），AUC_loo 0.71 (perm p 0.001)，LRT vs CM01 p0.03 | OS 无 | **OS HR 1.21 p0.006 beyond CM01，LRT p0.006**；应答 NR↑ 趋势 p0.07 |
| **LGALS1 myCAF→CD8** mean | 应答方向不稳 | **OS HR 1.98 p0.003 beyond CM08，LRT p0.002** | OS HR 1.28 p0.07（同向趋势）；应答 NR↑ 趋势 p0.06 |
| **collagen–DDR1 / COL1A1–DDR1** min/mean | 应答 mean OR 0.49 p0.044（NR↑），LRT p0.036，AUC_loo 0.63 (perm p 0.039) | OS 无 (HR 1.17 p0.28) | **应答无 (OR 1.24 p0.21)**；OS 无 |

- Hugo2016 (n=26) 的 collagen–DDR1"显著"（OR 6.7–271）是 **n=26 过拟合**，不可信，忽略。
- **collagen–DDR1 在最大、最规范的 IMvigor210 里对应答和 OS 都无效应**——和 desmoplasia_myCAF 一样。
  Riaz 黑素瘤里那点 NR 方向的弱信号（LRT p~0.04）**没有在 IMvigor210 复现**。
- ICI 受体端后果：IMvigor210 里 DDR1 beta_R = −0.35 (p 1e-13，免疫冷)，但**胶原×DDR1 交互为正**
  （+0.18，p 1e-5）——即尿路上皮里胶原高反而**削弱** DDR1 的免疫冷效应，**没有协同排斥**。
  与 4.2 的 LUAD/ESCA/THCA 相反 → **协同排斥是癌种特异的（肺/食管/甲状腺/胰腺），不在尿路上皮/黑素瘤**。

---

## 5. 证据矩阵（`tables/EVIDENCE_MATRIX.tsv`）

| axis | CM | 单样本复现 | 跨 GSE/癌种 | 空间支持 | 受体端后果 | ICI 应答 | OS 关联 | 独立于 parent CM？ | 文献 | 干预价值 |
|---|---|---|---|---|---|---|---|---|---|---|
| **FCN1单核→SPP1巨噬→内皮/淋巴** | CM01 | 66/66 + 121/121，稳健 0.95–0.99 | 23–30 GSE / 22–25 癌种 | Visium C1QA–CD93 4/5 阳 | （未单独测）| Riaz AUC 0.71；IMvigor NR 趋势 | **TCGA OS HR 1.23 (5e-21)；ICI OS HR 1.21 (0.006)** | **是**（加进模型后 CM01 残差转保护）| SPP1-TAM 泛癌不良预后，Cheng 2021/Bill 2023 | **高**：SPP1/CD44、TREM2、CSF1R |
| **myCAF galectin-1 (LGALS1) → CD8T 排斥** | CM08 | 291 样本 support 0.92，但 strict 稳健仅 0.42 | 41 GSE / 33 癌种 | Visium myCAF–CD8 共定位**阴/混杂**（与"排斥"相符）| （bulk：配体端主导）| Rose/IMvigor NR 趋势 p0.06 | **TCGA OS HR 1.10–1.15 (p<1e-4，令 CM08 归零)；Rose2021 ICI-OS HR 1.98 (0.003)** | **是**（LRT p<1e-4；配体端 p_L 8e-12）| CAF 来源 galectin-1 促 T 凋亡/排斥，Nambiar 2019 | **中-高**：galectin-1 抑制剂（临床期）|
| **myCAF collagen → 肿瘤 DDR1 免疫排斥轴** | CM07 | 115 样本 / 28 GSE / 22 癌种，但 strict 稳健仅 0.57（CM07 里最差）；LIANA 只有 COL1A1–DDR1 | 22 癌种 | **ccRCC CosMx**：非 juxtacrine（ccRCC caveat），但胶原+DDR1 高区 CD8 排斥 3/3 (p<1e-30)、间叶化 3/3 (p<1e-8)；Visium 未测 DDR1 | **DDR1 = 泛癌免疫冷标记**（TCGA + IMvigor p<1e-13）；胶原×DDR1 协同排斥仅 THCA/LUAD/ESCA/HNSC/THYM/PAAD | **无**（IMvigor210 应答/OS 均 null；Riaz 弱 NR 信号未复现）| **bulk 不独立于 desmoplasia_myCAF**（共线，符号翻转）| collagen-DDR1-SRC 致 CD8 排斥，Sun 2021 Nature (PDAC/NSCLC) | **中**：DDR1 激酶抑制剂 / 抗 DDR1 抗体——但仅限肺/食管/胰腺 |
| S100A8/A9–CD68 单核→巨噬 | CM08 | 317 样本 support 0.99 | 41/34 | Visium 5/5 阳 | — | — | — | **否**（纯髓系丰度）| — | 低（组成信号）|
| VCAN–ITGB1 单核↔myCAF | CM07 | — | — | **Visium 5/5 最强** | — | — | — | 未测 | 基质-髓系重塑 | 待评估 |

---

## 6. 三个问题的回答

### Q1. 哪些 CM 的临床效应能被一条/少数几条具体 communication axis 解释？

- **CM01**：能。CM01 作为 signature 本身预后效应很弱（TCGA OS HR 0.97，ns），但它的临床相关性
  几乎完全由 **FCN1单核→SPP1巨噬 轴**承载——SPP1 proxy 的 OS HR 1.23（p 5e-21）远超 signature，
  且加进模型后 CM01 残差**翻为保护**。可以把 CM01 重新表述为"SPP1-TAM 血管生成 niche"。
- **CM08**：部分能。CM08 的泛癌不良预后（OS HR 1.06，p 0.002）在加入 **myCAF galectin-1 (LGALS1) 轴** 后
  **完全归零**（HR 0.99, p 0.64），而 LGALS1 proxy 保持 HR 1.10–1.15（p<1e-4）——CM08 的预后信号
  可被 galectin-1 这一条轴吸收。ICI 场景 CM08 的信号（Rose2021 保护）也在加 LGALS1 后被拆开。
- **CM07**：**不能**用 COL1A1–DDR1 单独解释。CM07/desmoplasia 的强泛癌不良预后是**整个策展胶原程序**
  （desmoplasia_myCAF 21 基因）承载的；COL1A1–DDR1 proxy 与它共线、不独立。DDR1 单独作为免疫冷标记
  在肺/食管/胰腺有额外价值，但那是癌种特异的子结论。

### Q2. 哪些 LR 只是 cell-type abundance/composition 驱动，哪些控制 parent CM 后仍有独立信息？

- **纯组成/丰度驱动**（控制 parent 后无独立信息）：
  - S100A8/A9–CD68、S100A8–CD69（单核→巨噬/CD8）——就是"炎症髓系多"
  - MYL9–CD69、COL1A1–CD44（myCAF→单核/CD8）——就是"myCAF 多"
  - **collagen–DDR1 的 joint proxy**（bulk 里）——与 desmoplasia_myCAF 共线，条件后符号翻转
- **控制 parent CM 后仍有独立信息**：
  - **SPP1 macrophage 轴**：TCGA OS LRT p 5e-21、ICI-OS LRT p 0.006，两端（SPP1 配体 + CD44/整合素受体）都贡献
  - **myCAF galectin-1 (LGALS1) 轴**：TCGA OS LRT p 1.7e-5（parent 归零）、Rose2021 ICI-OS LRT p 0.002；
    信号在**配体端**（p_L 8e-12 >> p_R 2e-5）
  - **肿瘤 DDR1（受体端本身）**：作为免疫冷标记，在胶原之外独立（TCGA + IMvigor p<1e-13）；
    但胶原×DDR1 的**协同**只在 THCA/LUAD/ESCA/HNSC/THYM/PAAD 独立显著

### Q3. 只能选 1–3 条机制进主线 + 后续实验，选哪些？

1. **FCN1 monocyte → SPP1⁺ macrophage → 内皮/淋巴抑制（CM01）** — 首选。
   理由：单细胞复现最干净（66/66、121/121，跨 23–30 GSE）；是唯一在 **TCGA 泛癌 OS（HR 1.23, p 5e-21）
   和 ICI-OS（IMvigor210 HR 1.21, p 0.006）both 独立于 parent CM** 的轴；文献极强（SPP1-TAM 是公认的
   泛癌不良预后 / ICI 抵抗 TAM 状态）；干预靶点明确（SPP1–CD44、TREM2、CSF1R）。
   实验验证：SPP1 或 CD44 阻断 + 空间验证 FCN1→SPP1 的 niche 转变。

2. **myCAF galectin-1 (LGALS1) → CD8 T 排斥（CM08）** — 次选。
   理由：能**吸收 CM08 的全部预后信号**（加入后 parent 归零）；ICI-OS 在 2 个尿路上皮队列同向
   （Rose2021 HR 1.98 p0.003，IMvigor210 HR 1.28 p0.07）；空间上 myCAF–CD8 **不共定位**（与"排斥"机制自洽）；
   galectin-1 抑制剂已在临床。信号定位在**配体端**，机制清晰。
   实验验证：galectin-1 中和 / LGALS1 敲低的 CAF-CD8 共培养 + 类器官。

3. **myCAF collagen → 肿瘤 DDR1 免疫排斥轴（CM07）** — 有条件选，定位为**肺/食管/胰腺特异的
   desmoplasia 精细化**，不作为泛癌 ICI biomarker。
   理由：这正回答用户的核心疑问——"整体 fibrosis（desmoplasia_myCAF）在 ICI 队列几乎无预测力"是对的，
   而 **DDR1 这个更具体的受体端节点确实是稳健的泛癌免疫冷标记**（TCGA + IMvigor p<1e-13），
   胶原×DDR1 的协同排斥在 THCA/LUAD/ESCA/HNSC/THYM/PAAD 独立于总纤维化。CosMx 显示胶原+DDR1 高区
   CD8 空间排斥 + 肿瘤间叶化（3/3 组织 p 极小）。**局限**：bulk 里 joint proxy 不独立于 desmoplasia；
   ICI 预测性未在 IMvigor210 复现；空间验证只有 ccRCC（最差组织）。
   实验验证：**在肺癌 / PDAC 空间数据（非 ccRCC）**重做 CosMx 的 T1/T3/T4；DDR1 激酶抑制剂 + 抗 PD-1
   在纤维化肺癌模型。

---

## 7. 局限性

1. **Discovery 与 validation 未完全分离**：单样本 LIANA（discovery）用的是同一批泛癌单细胞；
   Visium / CosMx / TCGA / ICI 是独立数据，算 validation。但 CM 成员定义来自同一 CoVarNet，
   axis 选择基于 discovery 结果——存在选择效应，报告的 p 值不做多假设全局校正的"发现级"解释。
2. **CosMx 空间验证只有 ccRCC**，而 ccRCC 是 DDR1 轴最差的测试组织（正常近端小管 DDR1⁺ 混杂）。
   T1（非 juxtacrine）不能推广；可迁移的是 T3/T3b/T4（CD8 排斥 + 间叶化）。
   带 ICI 标签的 ccRCC CosMx（ZENODO DOI 10.5281/zenodo.16833780）仍 404。
3. **bulk proxy ≠ LIANA 验证**：全节反复强调。bulk 里 L 和 R 同高不能说明 source→target。
4. **LIANA consensus resource 限制**：DDR1 只与 COL1A1 配对，无法评估 COL1A2/COL3A1→DDR1。
5. **ICI 队列小**：Riaz n=71、Rose n=88、Hugo n=26（Hugo 结果已判为过拟合不可信）。
   只有 IMvigor210 (n=347) 够稳，而它对 collagen–DDR1 和多数 signature 都是 null。
6. **未做 NicheNet / LIANA+ 的严格 downstream target 推断**：受体端后果用的是 EMT/MMP 程序关联
   （CosMx）和 CYT/CD8 关联（bulk），不是配体特异的 regulon 活性。这是"功能一致性"证据，
   不是"通讯已发生"的证明。
7. **CM07 frozen projection 已知有 ccRCC 上皮污染**（记忆 [[phase-b-tcga-ici-progress]]），
   本报告统一用策展 `desmoplasia_myCAF` 作为纤维化背景，未把 CM07 frozen 投影当 myCAF readout。

---

## 8. 产物索引

| 文件 | 内容 |
|---|---|
| `tables/focus_axes_sample_level.tsv.gz` | 12 条目标 axis 的**逐样本** LIANA 行（含三档阈值 pass 标记）|
| `tables/focus_axes_replication_summary.tsv` | 每 axis 的复现 / 稳健性汇总（第 2.1 表来源）|
| `tables/focus_axes_by_gse.tsv` / `_by_cancer.tsv` | per-GSE / per-cancer 拆分 |
| `tables/collagen_ddr_sweep_summary.tsv` | 所有 collagen×DDR 配对扫描（证明只有 COL1A1–DDR1）|
| `tables/cosmx_cells_meta.tsv.gz` / `_counts_subset.tsv.gz` | ccRCC CosMx 细胞级坐标 + 46 基因 |
| `tables/cosmx_col_ddr1_T1..T4*.tsv` | CosMx 空间邻近 / CD8 排斥 / 受体端 EMT 检验 |
| `tables/bulk_proxy_TCGA_models.tsv` | TCGA Model A/B/C/D（OS + PFI）|
| `tables/bulk_proxy_TCGA_receiver_consequence.tsv` | 每癌种 `CYT/CD8 ~ L×R` |
| `tables/bulk_proxy_ICI_response_models.tsv` / `_survival_models.tsv` / `_receiver_consequence.tsv` | ICI 队列 |
| `tables/EVIDENCE_MATRIX.tsv` | 第 5 节证据矩阵（机器可读）|
| `reports/01_audit.json` / `02b_cosmx_spatial.json` | 方法元数据 |
| `scripts/01..04*` | 全部分析脚本 |
