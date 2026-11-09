# TME 注释 v2：谱系两轮分析、组织池重分配、低深度标签转移与多种子共识（2026-11-07）

本文件记录 [Lodi 式注释流程 v2](tme-annotation-v2-lodi-workflow-20261105.md) 在谱系阶段（B2–B6）实际采用的策略、参数、代码、结果和遇到的问题。v2 流程文件仍是设计依据；凡是这里的做法与它不同，以本文件为准，并在 §6 列出。运行细节见 [运行记录](../runs/20261107_annotation_v2_lineage_r2_r3_v1.md)。

标记约定（AGENTS §3）：**事实**＝有文件或日志可查；**解释**＝基于事实的推断；**假设**＝未验证；**建议**＝尚未执行。时间为 UTC，括号内为墨尔本时间（UTC+11）。

当前状态（2026-11-07 11:30 UTC）：
- T_NK、B/浆细胞、基质的第二轮命名和第二轮 QC 已审阅。
- 组织池已重分配；Myeloid 和 Endothelial 已合并重跑（第三轮），新簇尚未审阅。
- T_NK、B/浆细胞、基质的低深度标签转移已完成。
- 多种子共识第一步正在计算。
- **所有亚型标签都还不是最终版本**：§2.8 发现的聚类不稳定问题要靠多种子共识解决。

## 1. 总流程和细胞流向

```text
TME v1 run v2（3,347,659 个细胞）
 ├─ A0 统一双细胞旁证分数（Scrublet，只作簇级证据）
 ├─ A2/A3 各癌种大类（27 个癌种、642 个簇，51 处人工改写）
 ├─ B2 第一轮（QC 轮）：5 个谱系池
 │    → 移除双细胞样 102,770、低质量 8,833、错误谱系 121；搁置 20,893
 │    → 进入第二轮 2,676,921（其中低深度标记 553,641）
 ├─ B5 第二轮：5 个谱系池 + 4 个组织特异池（在 A4 中识别）
 ├─ 第二轮审阅：命名、第二轮 QC、组织池重分配（lineage_r2_reviewed.tsv，183 个簇）
 │    → 再移除 14,595，再搁置 10,072
 ├─ 第三轮：Myeloid（656,808）和 Endothelial（137,014）合并重跑
 ├─ 低深度标签转移（用户决定）
 └─ 多种子共识（用户决定；进行中）
```

| 阶段 | 输入细胞 | 去向 | 证据 |
|---|---|---|---|
| 第一轮 | 2,809,538（5 个池） | 见上 | [r1_action_summary.tsv](../inventories/annotation_v2/r1_action_summary.tsv) |
| 第二轮审阅 | 2,676,921 + 组织池 128,531 | 按簇判定 | [r2_action_summary.tsv](../inventories/annotation_v2/r2_action_summary.tsv) |
| 第三轮之后 | T_NK 1,409,873；Myeloid 656,808；B_Plasma 293,436；Stromal 236,895；Endothelial 137,014；Microglia 36,467；Alveolar_macrophage 8,048；Kupffer 857；Sinusoidal_EC 387 | — | `cell_labels_r3.tsv.gz`（§7） |

## 2. 各步骤：方法、参数和结果

所有聚类步骤共用 [`scripts/annotation_v2/major_by_cancer_v1.py`](../scripts/annotation_v2/major_by_cancer_v1.py)。设置环境变量 `LINEAGE_POOL` 和 `CELL_LABELS` 后，它在合并后的谱系池上运行。

### 2.1 公共的嵌入和聚类（事实，来自代码和日志）

| 步骤 | 实现 | 第一轮 | 第二、三轮 |
|---|---|---|---|
| 读取 | 按数据集 shard 读连续行块；基因为 run v2 的 HVG 公共集、注册表基因和状态基因（12,782 个） | 同左 | 同左 |
| HVG | 原始计数上做 seurat_v3，按来源分批后合并排名；排除 MT-、RPS/RPL、MALAT1、血红蛋白、IG V、TCR V 基因。loess 奇异时把小来源合并（`MIN_SRC` 200 起，每次翻倍）后重试 | 2,000 | 3,000（`N_HVG=3000`） |
| 状态回归 | 用 scanpy `score_genes` 打分后做 OLS 残差。S/G2M 用 Seurat cc.genes 2019；IFN 用 MSigDB BROWNE；stress 用 V7 `Stress_response`（Lodi 的 stress 基因表不公开，以此替代）；hypoxia 用 Lodi 15 基因。不回归样本和癌种 | 按谱系，见下 | 同左 |
| 缩放和 PCA | z 分数截断在 ±10；PCA 50 维，arpack | 同左 | 同左 |
| PC 数 | 方差比曲线的 Kneedle 拐点，截断在 `[PC_MIN, 40]` | `PC_MIN=15`，所有池都落在 15 | `PC_MIN=30`，所有池都落在 30 |
| 批次校正 | harmonypy 2.0.0，按 `source` 校正，θ=2，λ=1，`max_iter_harmony=20`，其余用默认值（epsilon_cluster 0.001、epsilon_harmony 0.01，**会提前停止**，见 §5 问题 5） | 同左 | 同左 |
| 近邻和聚类 | 在 X_emb 上做 15-NN；Leiden（igraph，`n_iterations=2`）；UMAP min_dist 0.3 | 分辨率 0.3/0.6/1/1.5 | 分辨率 0.3/0.5/0.7/1/1.5/2/3 |
| 主分辨率 | **固定为 1**（`MAIN = 1.0`），只在主分辨率上算 DE 和证据表 | 同左 | 同左 |
| DE | Seurat 式 Wilcoxon，每簇最多抽 2,000 个细胞，Bonferroni 校正 | 同左 | 同左 |
| 证据 | 只在测到该基因的细胞中计算检出率（`expression.tsv`）；V7 证据引擎在合并面板 v2 上逐簇给出所有候选（不选赢家）；附 QC、双细胞 z 分数、恶性判定、来源构成 | 同左 | 同左 |
| 随机种子 | — | 2026110603 | 第二轮 2026110701；第三轮 2026110702 |

各谱系池的状态回归变量（事实，作业日志第一行）：

| 池 | 回归变量 |
|---|---|
| T_NK、Stromal、Endothelial、Sinusoidal_EC | S、G2M、IFN、stress |
| Myeloid、Kupffer、Microglia、Alveolar_macrophage | S、G2M、IFN、hypoxia |
| B_Plasma | S、G2M、IFN、stress、hypoxia |

### 2.2 第一轮（QC 轮）

逐簇判定见 [lineage_r1_reviewed.tsv](../inventories/annotation_v2/lineage_r1_reviewed.tsv)，应用脚本是 [apply_lineage_r1_review.py](../scripts/annotation_v2/apply_lineage_r1_review.py)（commit 84c00103）。判据如下：
- **低质量：** 相对同池其他簇判断，不用固定阈值。
- **双细胞样：** 同时要有其他谱系 marker 和双细胞 z 分数偏高；只满足一项的记为 `set_aside_review`。
- **低深度：** MT 和污染都正常、只是没有 marker 的簇，保留并标记 `keep_flag_low_depth`。

### 2.3 第二轮（干净对象重建）

参数见 §2.1。除 5 个谱系池外，A4 识别出的组织特异群体（Kupffer、Microglia、Alveolar_macrophage、Sinusoidal_EC）作为独立池首次聚类。各池簇数（分辨率 1）：T_NK 17、B_Plasma 19、Myeloid 24、Stromal 29、Endothelial 16、Alveolar_macrophage 18、Microglia 24、Kupffer 18、Sinusoidal_EC 16。

### 2.4 第二轮审阅（命名、第二轮 QC、组织池重分配）

逐簇的判定和理由都在 [lineage_r2_reviewed.tsv](../inventories/annotation_v2/lineage_r2_reviewed.tsv)，由 [apply_lineage_r2_review.py](../scripts/annotation_v2/apply_lineage_r2_review.py) 写到细胞上，输出 `cell_labels_r3.tsv.gz`。每个细胞保留 `origin_pool`、`r2_cluster`、`r2_label`、`r2_action`，所以每一步都可以撤回。

判定规则（审阅前确定，对所有池相同）：

| 动作 | 规则 |
|---|---|
| `keep` | 有该谱系明确的身份 marker |
| `keep_state` | 主要由状态程序驱动（增殖、IFN、热休克、即刻早期、金属硫蛋白、缺氧）；状态单独记录，不作为亚型 |
| `keep_flag_low_depth` | 低深度、没有 marker，MT 和污染正常；进入标签转移（§2.6） |
| `keep_flag_source_restricted` | 来自 1–2 个来源，但不属于单样本 |
| `keep_flag_sample_specific` | **单个样本占比 ≥0.70**；计入谱系总数，但不作为泛癌亚型 |
| `remove_doublet_like` | 有其他谱系 marker，且双细胞 z ≥1.5（第一轮接受的最低值约为 1.7） |
| `remove_low_quality` | 相对同池 MT 偏高或只有核内长转录本，又没有身份程序 |
| `set_aside_review` | 只有其他谱系基因、没有双细胞证据；或怀疑是残留肿瘤（如骨肉瘤的成骨样基质） |
| `keep_tissue` | 组织特异池中**有该组织特异 marker** 的簇（见下表） |
| `move_to_<pool>` | 身份属于另一个池；目标池合并重跑，或者记入目标谱系总数但不定亚型 |

组织身份的判据（事实，计数层检出率）。这些池在 A2/A3 时是按癌种的整个簇命名的，结果命名偏宽（**解释**）：

| 池 | 特异 marker | 保留 | 占池比例 | 移回主池 |
|---|---|---|---|---|
| Kupffer | TIMD4 0.48、CETP 0.76、MARCO 0.69、VCAM1 0.93（只有一个簇 TIMD4 >0.1） | 857 | 5.6% | Myeloid 13,593 |
| Sinusoidal_EC | CLEC4G 0.85、CLEC4M 0.67、CLEC1B 0.80、FCN2 0.87（其他簇都 <0.3） | 387 | 7.3% | Endothelial 4,695 |
| Alveolar_macrophage | MARCO 0.97、FABP4 0.74、PPARG 0.83、INHBA 0.70（只有一个簇 FABP4 >0.3） | 8,048 | 13.2% | Myeloid 51,381 |
| Microglia | P2RY12 ≥0.35 且 TMEM119 ≥0.35，同时 TGFBI、MRC1 不高；低检出簇留在池内转移 | 24,106 + 低深度 12,361 | 79.6% | Myeloid 8,798（TGFBI 0.72–0.91、F13A1、MRC1 高） |

AM 的部分符合簇（14、8、9）和 Kupffer 样 VCAM1⁺ 簇（5、3、1）没有达到判据，也回了 Myeloid。它们的 `origin_pool` 保留，可以撤回。

第二轮的主要命名（细胞数为该簇在第二轮对象中的细胞数）：

| 池 | 亚型（第二轮，分辨率 1） |
|---|---|
| T_NK | CD8 效应/耗竭样 201,638；CD4 记忆 IL7R 样 197,669；T naive/CM 样 177,706；Treg 样 152,968；CD8 记忆 ANXA1/GZMK 样 153,063；NK 细胞毒样 105,967；NK XCL1/KLRC1 样 73,106；CD8 GZMK 样 61,390；CD4 CXCL13/Tfh 样 60,285；增殖 T 30,853；低深度 177,709 |
| B_Plasma | B 记忆活化样 53,890；B naive 样 40,902；浆细胞 39,074 + 34,841；B 非典型记忆样 21,118；生发中心 B 6,726；生发中心增殖 6,745；低深度 B 48,136、浆细胞低深度 19,765、B 带环境 RNA 6,397 |
| Stromal | CAF 基质样 33,012；周细胞 RGS5/STEAP4 29,641；血管平滑肌 MYH11 26,840；周细胞 COL4A1/NOTCH3 20,664；炎症缺氧 CAF 14,107；外膜 APOD 样 13,978；PI16 通用成纤维 13,852；CCL19/C3 成纤维 12,463；上皮下成纤维 8,261；CAF 基质 COL12A1 5,662；内脏平滑肌 ACTG2 1,815；低深度 29,321 + 9,861 + 4,364 |

第二轮新增的移除和搁置（事实）：
- **移除：** 双细胞样 12,282；低质量 2,313。
- **搁置：** 10,072。
- **移到其他池：** Myeloid 的 B 细胞簇 2,384 个记入 B_Plasma 总数（B 池不重跑，所以没有亚型）；B 池里的 mregDC 1,014 个移到 Myeloid。

### 2.5 第三轮：Myeloid 和 Endothelial 合并重跑

参数与第二轮相同（§2.1），种子 2026110702，代码 commit 0764824a（worktree），作业 30677658（Myeloid）、30677659（Endothelial）。
- **Myeloid：** 656,808 个细胞，30 个 PC，22 簇，耗时 3 小时 04 分。
- **Endothelial：** 137,014 个细胞，30 个 PC，19 簇，耗时 2 小时 09 分。

新簇**尚未审阅**。组织池移来的细胞按 `origin_pool` 可查。

### 2.6 低深度标签转移（用户决定，2026-11-07）

用户的决定：低深度细胞保留在谱系总数里，通过标签转移得到亚型，并标记为推断；这些亚型不计入主分析的亚型比例，只用于敏感性分析。

代码：[transfer_low_depth_v1.py](../scripts/annotation_v2/transfer_low_depth_v1.py)。

| 项 | 取值 |
|---|---|
| 查询细胞 | 第一轮标记低深度的细胞，加上第二轮 `keep_flag_low_depth` 簇中的细胞；已移出本池的不算 |
| 参考细胞 | 同池中 `keep` 或 `keep_tissue` 簇里、没有低深度标记的细胞；不用状态簇、单来源簇、单样本簇 |
| 空间 | 第二轮对象的 X_emb（Harmony，30 维） |
| 近邻 | pynndescent，只在参考细胞上建索引，k=30；查询细胞之间不互为近邻 |
| 投票 | 距离倒数加权；置信度 = 得票最多的标签所占权重 |
| 判定 | 置信度 ≥0.6 记为 `transferred`，否则记为 `<pool>_unresolved_low_depth` |
| 验证 | 每个标签中基因数最低的 10% 参考细胞留出，用其余参考细胞预测。留出细胞仍比真实查询细胞深，所以**准确率是上限** |
| 种子 | 2026110704 |

结果（事实，见 `inventories/annotation_v2/r3_checks/transfer_*`）：

| 池 | 参考 | 查询 | 未定 | 留出验证宏平均准确率 | 置信细胞上的准确率 | 查询细胞与自身身份簇一致 |
|---|---|---|---|---|---|---|
| T_NK | 1,007,305 | 364,196 | 19.3% | 0.919（留出细胞中位基因数 629） | 0.966 | 0.860 |
| B_Plasma | 162,338 | 111,088 | 9.5% | 0.907（647） | 0.937 | 0.964 |
| Stromal | 170,350 | 53,971 | 21.6% | 0.888（903） | 0.939 | 0.812 |

最弱的标签：B 非典型记忆样 0.703，Stromal 内脏平滑肌 0.725、上皮下成纤维 0.741。T_NK 的查询细胞中最多的去向是 T naive/CM 样，有 155,680 个。

Myeloid、Endothelial 等第三轮审阅完成后再做；Microglia 在自己的对象内做，参考集是小胶质和胶质瘤巨噬两类。

### 2.7 批次检查：LISI 和 Harmony 提前停止

代码：[lineage_batch_and_rare_checks_v1.py](../scripts/annotation_v2/lineage_batch_and_rare_checks_v1.py)，作业 30677660。
- **LISI：** 按 Korsunsky 2019 的定义，困惑度 30，在 200,000 个随机细胞上取 90-NN。
- **Harmony 对照：** 设 epsilon 为 −inf，跑满 20 次迭代，不提前停止。

| 池 | 来源 LISI 中位数：PCA → 当前 Harmony（→ 不提前停止） | 样本 LISI 中位数：PCA → 当前 Harmony（→ 不提前停止） | 来源数 |
|---|---|---|---|
| T_NK | 2.64 → 7.73 | 6.07 → 10.38 | 70 |
| Myeloid（第二轮） | 2.08 → 6.12（→ 6.19） | 4.48 → 8.52（→ 8.59） | 70 |
| B_Plasma | 1.92 → 6.57 | 4.12 → 9.00 | 68 |
| Endothelial（第二轮） | 1.46 → 6.26（→ 6.35） | 3.37 → 8.28（→ 8.35） | 64 |
| Stromal | 1.17 → 4.86 | 2.39 → 6.62 | 69 |

**解释：**
- 不提前停止对批次混合的改善在 0.1 以内。
- B_Plasma 的 LISI 不比其他池低。所以 Lodi 为 B 细胞改用 CCA 的理由，在本数据上暂时不成立（只按 LISI 判断）。
- Stromal 的 LISI 最低，和它有很多单样本簇是一致的（§5 问题 10）。

### 2.8 聚类稳定性（关键发现）

代码：[leiden_stability_check_v1.py](../scripts/annotation_v2/leiden_stability_check_v1.py)，作业 30678135。做法是在第二轮对象上重算近邻和分辨率 1 的 Leiden，再与运行中的簇比较。比较分两层：簇层面用 ARI；标签层面先把每个新簇映射为多数审阅标签，再看细胞标签一致的比例。

| 池 | 变体 | 簇数 | 簇 ARI | 标签 ARI | 标签一致的细胞 |
|---|---|---|---|---|---|
| Endothelial | 同一嵌入、同一种子重算 | 17 | 0.854 | 0.898 | 0.951 |
| Endothelial | 同一嵌入，换种子（2 次） | 16、18 | 0.507、0.442 | 0.501、0.513 | 0.682、0.688 |
| Endothelial | Harmony 不提前停止（2 个种子） | 18、16 | 0.401、0.471 | 0.478、0.473 | 0.659、0.641 |
| Stromal | 同一嵌入、同一种子重算 | 29 | 0.889 | 0.889 | 0.929 |
| Stromal | 同一嵌入，换种子（2 次） | 30、28 | 0.520、0.498 | 0.549、0.549 | 0.696、0.699 |
| Stromal | Harmony 不提前停止（2 个种子） | 32、30 | 0.469、0.491 | 0.504、0.501 | 0.649、0.653 |

**解释：**
- 只换 Leiden 种子，就有约 30% 的细胞亚型标签会变。原因是分辨率 1 的分区本身不稳定，不在 Harmony。
- 同一种子重算也只有 0.93–0.95，说明近邻图本身也有随机性（pynndescent）。
- Harmony 不提前停止的标签一致性比只换种子低约 4–5 个百分点，可能有一点额外影响，但落在种子波动的范围内。

**结论：** Harmony 保留当前设置（建议，见 §6）。单次聚类得到的亚型不能直接定稿，改用多种子共识（§2.10）。

### 2.9 稀有亚型（T_NK，分辨率 3）

数据：`r3_checks/rare_gene_detection.tsv`、`r3_checks/resolution_crosswalk.tsv`（计数层检出率）。

| 分辨率 3 簇 | 主簇 | 证据 | 拟用名称 |
|---|---|---|---|
| 16（56,658） | 11 NK 细胞毒 | CD3D 0.75、CD8A 0.45、FGFBP2 0.66 | CD8 TEMRA 样，从 NK 中拆出 |
| 36（51,185） | 11 | FGFBP2 0.76、KLRF1 0.67、CD3D 0.12 | NK 细胞毒 |
| 4（53,355） | 7 | KLRB1 0.81、SLC4A10 0.14、ZBTB16 0.15、RORC 0.12、IL23R 0.11 | MAIT/Th17 样；SLC4A10 检出低，不单称 MAIT |
| 38（4,167） | 10 | KIT 0.47、IL23R 0.23、RORC 0.19、KLRB1 0.85、CD3D 0.08 | ILC3 样 |
| 37（815） | 0 | CCR7 0.85、SELL 0.73、CD8A 0.79 | CD8 naive 样 |
| 42（2,383） | 10 | TRDC 0.83、MKI67 0.82 | 增殖 NK（状态） |

做不到的细分（事实）：
- **γδ T：** TRDC 在 NK 簇里也检出（0.62），而能区分它们的 TRDV/TRGV 基因在选基因时被排除了。
- **Tfh：** CXCR5 在所有簇中检出都 ≤0.14，只能并称 CXCL13⁺ Tfh/Tph 样。

### 2.10 多种子共识（进行中）

第一步：[consensus_seed_clusters_v1.py](../scripts/annotation_v2/consensus_seed_clusters_v1.py)。
- **做法：** 每个池用 10 个种子（2026110800–2026110809），每个种子都重算 15-NN 和分辨率 1 的 Leiden，只保存簇编号。
- **作业：** 30678224（T_NK）、30678225（B_Plasma）、30678226（Stromal）、30678227（Endothelial，第三轮对象）、30678229（Myeloid，第三轮对象）。

第二步的规则在看到结果之前写定（**建议，尚未执行**）：
1. 每个种子的簇映射为多数审阅标签。
2. 两个标签之间互换的细胞比例取双向平均，≥0.25 的就合并成一个更粗的亚型（用并查集）。
3. 合并后，每个细胞的稳定性等于与众数标签一致的种子比例；≥0.7 为稳定。不稳定的细胞和低深度细胞同样处理：计入谱系总数，不计入主分析的亚型比例。
4. 稀有亚型按分辨率 3 的簇覆盖，只用于观测细胞（不用于转移细胞）。

## 3. 参数总表

| 参数 | 值 | 位置 |
|---|---|---|
| HVG 数 | 第一轮 2,000；第二、三轮 3,000 | `N_HVG` |
| PC 下限和上限 | 第一轮 15–40；第二、三轮 30–40 | `PC_MIN`，代码中上限 40 |
| Harmony | 按来源，θ=2，λ=1，最多 20 次，默认 epsilon | 代码 |
| 近邻 | 15 | 代码 |
| 分辨率 | 第一轮 0.3、0.6、1、1.5；第二、三轮 0.3、0.5、0.7、1、1.5、2、3；主分辨率 1 | `RESOLUTIONS`、`MAIN` |
| 双细胞阈值 | z ≥1.5，且有其他谱系 marker | 审阅规则 |
| 单样本阈值 | 最大样本占比 ≥0.70 | 审阅规则 |
| 标签转移 | k=30；置信度 ≥0.6；留出最低 10% 基因数 | `transfer_low_depth_v1.py` |
| 共识 | 10 个种子；合并阈值 0.25；稳定阈值 0.7 | §2.10 |
| 资源 | 实测内存峰值：T_NK 111 G、Myeloid 56 G、B 28 G、Stromal 27 G、Endothelial 15 G | sacct |

## 4. 代码清单

| 脚本 | 作用 | 用法 |
|---|---|---|
| [doublet_scores_v1.py](../scripts/annotation_v2/doublet_scores_v1.py) | A0：每个单位统一的 Scrublet 分数 | 见脚本头 |
| [major_by_cancer_v1.py](../scripts/annotation_v2/major_by_cancer_v1.py) | A2 和所有谱系轮次的嵌入、聚类、证据 | `RUN_DIR CANCER OUT_ROOT REGISTRY_DIR DOUBLET_DIR REPO SEED`，环境变量 `LINEAGE_POOL`、`CELL_LABELS`、`STATE_REGRESS`、`N_HVG`、`PC_MIN`、`RESOLUTIONS` |
| [propose_major_labels_v1.py](../scripts/annotation_v2/propose_major_labels_v1.py)、[review_major_labels_v1.py](../scripts/annotation_v2/review_major_labels_v1.py) | A3 大类提议和审阅 | 见脚本头 |
| [apply_lineage_r1_review.py](../scripts/annotation_v2/apply_lineage_r1_review.py) | 第一轮判定写到细胞上 | `R1_DIR REVIEW_TSV CELL_MAJOR_LABELS OUT_CELLS OUT_SUMMARY` |
| [apply_lineage_r2_review.py](../scripts/annotation_v2/apply_lineage_r2_review.py) | 第二轮判定和组织池重分配写到细胞上 | `R2_DIR REVIEW_TSV CELL_LABELS_R2 OUT_CELLS OUT_SUMMARY` |
| [transfer_low_depth_v1.py](../scripts/annotation_v2/transfer_low_depth_v1.py) | 低深度标签转移和留出验证 | `OBJ_H5AD REVIEW_TSV POOL CELL_LABELS OUT_PREFIX SEED`，环境变量 `COLLAPSE`、`EXTRA_REF_ACTIONS` |
| [lineage_batch_and_rare_checks_v1.py](../scripts/annotation_v2/lineage_batch_and_rare_checks_v1.py) | LISI、Harmony 不提前停止对照、高分辨率基因检出 | `OBJ_ROOT OUT_DIR SEED POOL...`，环境变量 `HARMONY_POOLS` |
| [leiden_stability_check_v1.py](../scripts/annotation_v2/leiden_stability_check_v1.py) | Leiden 种子和 Harmony 稳定性对照 | `OBJ_ROOT REVIEW_TSV OUT_TSV POOL...` |
| [consensus_seed_clusters_v1.py](../scripts/annotation_v2/consensus_seed_clusters_v1.py) | 多种子共识第一步 | `OBJ_H5AD OUT_TSV_GZ N_SEEDS BASE_SEED` |

环境：`${PC_ROOT}/envs/conda/scanpy`（Python，scanpy、harmonypy 2.0.0、pynndescent 0.5.13）。作业在 WEHI Milton 的 bigmem 分区运行，sbatch 由 `export` 后的环境变量驱动，因为 `sbatch --export` 会按逗号拆值。

## 5. 遇到的问题与处理

| # | 问题 | 证据 | 处理 | 状态 |
|---|---|---|---|---|
| 1 | R 环境没有 DoubletFinder | A0 记录 | 改用 Scrublet，按 10x 预期率；只作簇级旁证 | 已记录的偏离 |
| 2 | A2/A3 的组织特异群体命名偏宽：特异 marker 只支持 5.6–28% 的细胞（第二轮审阅前的初始比例；Microglia 规则写定后保留更多，见 §2.4） | §2.4 | 按组织特异 marker 保留，其余回主池并合并重跑 | 已处理；A4 的判据应改为按簇的特异 marker（建议） |
| 3 | 第一轮 Kneedle 拐点在所有池都落在下限 15，亚型分辨率不够 | 作业日志 `n_pc` | 第二轮下限提到 30；第二、三轮又全部落在 30 | 部分解决：拐点法在这里没有区分度（解释） |
| 4 | seurat_v3 loess 在小来源上奇异 | 第一轮 B 失败（30673484）；第二轮 T_NK 和 B 用到 `min_src=400` | 合并小来源后重试，阈值写入日志（0701cbca） | 已解决 |
| 5 | v2 方法文件 B5 写“Harmony 不提前停止”，代码用了默认 epsilon，实际 2–4 次迭代就停了 | Harmony 日志 | 做了对照（§2.7、§2.8）：批次混合几乎相同，稳定性不更好 | 已检查；建议保留当前设置并修正方法文件 |
| 6 | 原流程没有计算 LISI | 代码 | 补算（§2.7） | 已解决 |
| 7 | 主分辨率写死为 1，没有按相邻 ARI 选 | 代码 `MAIN = 1.0` | 由多种子共识和稀有亚型高分辨率补标代替 | 进行中 |
| 8 | **分辨率 1 的 Leiden 不稳定：换种子约 30% 细胞标签会变** | §2.8 | 多种子共识（用户决定） | 进行中 |
| 9 | 低深度细胞多：第一轮标记 553,641 个（20.7%），第二轮又形成无 marker 簇（如 Myeloid 簇 12，占 9.5%） | r1/r2 汇总 | 标签转移，不计入主分析亚型比例（用户决定） | T/B/Stromal 已完成 |
| 10 | 单样本簇：Harmony 只按来源校正，管不到同一来源内的样本差异（如 Kupffer 池 4 个簇 70–89% 来自 GSE299340 的单个样本；Stromal 有 9 个簇单样本占比 ≥0.74） | `clusters.tsv`、obs | 单样本规则（≥0.70）：计入谱系总数，不作为亚型 | 规则已用；是否改用样本层面的批次变量未决 |
| 11 | 环境 RNA 抬高检出率（Kupffer 池 ALB 在多数簇 >0.6，AM 池 SFTPC，Microglia 池 SPP1 近 1） | §2.4 的计数表 | 检出率只作簇间相对比较；命名要求多个特异 marker 同时出现 | 已在判据中考虑 |
| 12 | γδ T 不可分：TR V 基因被排除；MAIT 的 SLC4A10 检出低；Tfh 的 CXCR5 检出低 | §2.9 | 只给保守名称 | 未解决；可选做法是对 γδ 单独读回 TRDV/TRGV（建议） |
| 13 | 只有核内长转录本的簇处理不一致：第一轮 T_NK 簇 1 作为低质量删除，第二轮 Kupffer 簇 16、1、11 和 Sinusoidal 簇 9 跟随所在池 | 审阅表 | 第二轮中这些簇多为单样本，随单样本规则处理 | 未统一；需要按来源的 assay 类型（单核 vs 单细胞）复核（建议） |
| 14 | 标签转移脚本第一次运行时读错标签列（`iloc[:, 0]` 取到了 `pool`），留出验证得到全 1.0 的假准确率 | 作业 30677663 | 发现后取消，改为按列名 `r2_label` 读取，重跑（30677673） | 已解决；§2.6 的结果来自修复后的版本 |
| 15 | 第二轮作业申请 480 G 内存，同一用户的配额（QOSMaxMemoryPerUser）让 3 个作业排队；实测峰值只有 15–111 G | sacct | 第三轮改为 160 G、检查作业 100–240 G | 已解决 |
| 16 | Myeloid 中错放的 B 细胞（2,384）计入 B 总数，但 B 池不重跑，所以没有亚型 | 审阅表 | 标记 `B_misplaced` | 已接受的限制 |
| 17 | `keep_flag_low_depth` 簇中的细胞，以及身份簇中带第一轮低深度标记的细胞，查询时没有区分：后者也会被转移 | 代码 | 转移结果与自身身份簇的一致率是 0.81–0.96，作为旁证 | 已记录 |

## 6. 偏离 Lodi 设计和 v2 方法文件的地方

| 项 | v2 方法文件 | 实际 | 理由 |
|---|---|---|---|
| Harmony 停止条件 | 不提前停止 | 默认 epsilon，会提前停止 | §2.7 和 §2.8：批次混合和稳定性都没有因此变差。**建议**修正 B5 的写法 |
| B 细胞 CCA | LISI 不足时改用 CCA | 没有改 | B 的来源 LISI 中位数 6.57，在各池中属于高的 |
| 主分辨率 | 按相邻 ARI 和可解释性选 | 固定为 1，加多种子共识和高分辨率补标 | §2.8 |
| 组织特异群体 | A4 在各癌种内识别后移出 | 加了第二轮审阅：只按特异 marker 保留 | §5 问题 2 |
| 低深度细胞 | 没有规定 | 标签转移，计入谱系总数，不计入主分析亚型比例 | 用户决定（2026-11-07） |
| 单样本簇 | 没有规定 | 最大样本占比 ≥0.70 的簇不作为亚型 | §5 问题 10 |

## 7. 产物位置（不入库；以 `${PC_ROOT}` 为根）

`${PC_ROOT}/tme_v1` 目前是指向 scratch 的符号链接，scratch 上的数据 30 天不访问会被清理（见 `TEMP_LINKS_TO_SCRATCH.txt`）。

| 产物 | 路径 | SHA256 前 16 位 |
|---|---|---|
| 第二轮细胞标签 | `tme_v1/annotation_v2/lineage_r2/cell_labels_r2.tsv.gz` | 5f876fed3c1e8426 |
| 第三轮细胞标签 | `tme_v1/annotation_v2/lineage_r3/cell_labels_r3.tsv.gz` | 74b0973cf492e2fc |
| 第二轮对象（9 个池） | `tme_v1/annotation_v2/lineage_r2/<pool>/object.h5ad` | — |
| 第三轮对象 | `tme_v1/annotation_v2/lineage_r3/{Myeloid,Endothelial}/object.h5ad` | — |
| 检查和转移结果 | `tme_v1/annotation_v2/lineage_r3/{checks,transfer,consensus}/` | 汇总表已复制到 [inventories/annotation_v2/r3_checks/](../inventories/annotation_v2/r3_checks/) |

第二轮的结果表和细胞标签已经同步到旧服务器（`object.h5ad` 按用户决定不传），见运行记录。

## 7A. 四个组织池的十种子共识（2026-11-09，用户决定）

用 2026-11-08 接续方法的同一 `consensus()`、同一 10 个种子（2026110800–809）、合并阈值 0.25、稳定阈值 0.7，作用于第二轮对象和第二轮审阅标签。作业 30716245–248，全部 COMPLETED。脚本：[consensus_tissue_pools_v1.py](../scripts/annotation_v2/consensus_tissue_pools_v1.py)；汇总：[tissue_consensus_v1](../inventories/annotation_v2/tissue_consensus_v1/)。

| 池 | 对象细胞 | 保留 | 保留且稳定 | 符合主分析条件 |
|---|---:|---:|---:|---:|
| Kupffer | 15,305 | 857 | 848 | 846 |
| Microglia | 45,837 | 36,467 | 34,196 | 21,982 |
| Alveolar_macrophage | 61,053 | 8,048 | 7,809 | 7,788 |
| Sinusoidal_EC | 5,336 | 387 | 386 | 386 |

- 工作版里 `observed_tissue_no_seed_consensus` 的33,398个细胞中，31,002个符合条件；工作版文件尚未更新。
- 已观察：Microglia_proliferating（355个）只有14个保持原标签且都不稳定；Microglia_S100A8 有728个（约26%）众数归入 Microglia。
- 边界：Kupffer、AM、Sinusoidal_EC 各只有一个保留标签，这里只检验该簇边界，不验证其内部亚型。这是工作结果，不是最终注释。

## 8. 未完成和下一步（建议）

1. 多种子共识第二步：按 §2.10 写定的规则合并亚型、计算每个细胞的稳定性。
2. 审阅第三轮 Myeloid 和 Endothelial 的新簇，写 `lineage_r3_reviewed.tsv`；然后做这两个池和 Microglia 的标签转移。
3. 稀有亚型覆盖（§2.9），生成最终细胞标签，并写明每个标签的状态：观测、转移、未定、状态、单样本。
4. 修正 v2 方法文件 B5 的 Harmony 说明。
5. E1：在有十月标签的细胞上，与十月注释做列联表和 ARI/NMI 对照。
6. 把第三轮结果和最终标签同步到旧服务器，并把 `tme_v1` 从 scratch 移到 labshare（复制前先查 inode 配额）。
