# 运行记录：20260731_copykat_infercnv_intersection_v1

## 运行身份

| 字段 | 值 |
| --- | --- |
| Scope | `Project_v3` 恶性细胞识别、CoVarNet 与细胞通讯扩展 |
| Run type | production / retry / append-only attempts |
| Initial start | 2026-07-31 18:06 +08:00 |
| Active v2 restart | 2026-08-05 05:11 +08:00 |
| Snapshot | 2026-08-12 10:06 +08:00 |
| Status | `interrupted_stalled_nonterminal_stage01` |
| Source commit | `not_recorded`；Project_v3 路径未检测到 Git 元数据 |
| Environments | `scanpy`、`scrna_r` |
| Active manifest | `samples_analysis_v2.tsv` |

逻辑根 `${CM_RUN_ROOT}` 指向：
`${RUN_ROOT}/09_malignancy_communication/20260731_copykat_infercnv_intersection_v1`。

## 目标

在最终证据注释对象和原始计数 Zarr 上按生物学样本运行 inferCNV、CopyKAT 与 scMalignantFinder，以 CopyKAT 和 inferCNV 的严格交集定义恶性细胞；随后重建恶性状态输入的 CoVarNet，按样本和 CM 运行 LIANA，并用本地 NicheNet v2 ligand-target 模型验证优先 ligand-receptor 对。

## 输入与分析分母

源对象：31 个 GSE、1,963,745 个细胞、773 个原始 sample units、11,063 个基因。

当前分析分母经过两项有审计的范围处理：

1. 根据 2026-08-03 的既有用户授权，下游排除 GSE155698 和 GSE116256：共 75 个原 sample units、114,332 个细胞；原始表达对象与全量注释对象未删除。
2. GSE166555 从 2 个技术性 `orig.ident` 分组重组为 25 个生物学样本。该操作保持 48,819 个细胞和 global positions 完整不变，使 active samples 从 698 变为 721。

因此当前 v2 清单为 721 个生物学样本、1,849,413 个分析细胞。GSE166555 重组清单 SHA256 已记录在源运行审计中。

## 候选与参考细胞规则

- 默认实体瘤：Epithelial 为肿瘤候选。
- Osteosarcoma：Epithelial 与 Fibroblast 为候选。
- Melanoma：Epithelial、Neural 与 Fibroblast 为候选。
- Acute myeloid leukemia：HSPC、Monocyte、Dendritic、Granulocyte、Erythroid 为候选。
- 参考大类为 Bcell、Dendritic、Endothelial、Macrophage、Mast、Monocyte、NK、Plasma、Tcell。
- 分析必须按样本使用相同原始 global positions 对齐 counts、注释和逐样本证据。

## 分析流程、方法与参数

### Stage 00：清单与输入审计

- 从最终注释 H5AD 构建逐样本 manifest 和每样本 NPZ 位置清单。
- 记录 source/analysis/excluded 分母、GSE 排除、候选/参考计数和 hash。
- GSE166555 v2 重组要求 48,819/48,819 个 cell ID 可解析、global positions 唯一、细胞覆盖不变，并与原 obs 生物学 sample ID 计数一致。

### Stage 01：逐样本恶性判定

每个样本按以下顺序处理：

```text
从 counts Zarr 提取 raw counts
→ infercnvpy CNV score
→ CopyKAT aneuploid/diploid
→ scMalignantFinder 辅助概率
→ 严格 inferCNV-positive ∩ CopyKAT-aneuploid 交集
→ append-only attempt 记录与完整性检查
```

inferCNV：

- engine：infercnvpy。
- gene annotation：hg38 GENCODE v27 位置表。
- window size：250。
- reference quantile：0.99。
- CNV-positive threshold：每样本参考细胞 CNV score 的 99% quantile，因此阈值随样本变化。
- 候选或参考细胞不足时记录 `skipped_insufficient_cells`，不伪造判定。

CopyKAT：

- genome：`hg20`，按当前实际命令原样记录；正式解释前应核对该值是否为项目有意约定。
- ID type：symbol (`S`)。
- `ngene.chr=5`。
- `win.size=25`。
- `KS.cut=0.1`。
- distance：Euclidean。
- cores：2；OMP/OpenBLAS/MKL/NumExpr 默认限制为 2 threads。
- 每 60 秒记录进程状态、CPU、RSS、swap、线程和最后一条已刷新输出。

scMalignantFinder：

- 只评估候选细胞。
- library-size normalization 10,000，不做 log1p。
- 模型特征 2,707；缺失特征补零，并严格按 `ordered_feature.tsv` 排序。
- 作为独立参考，不改变 CopyKAT 与 inferCNV 的严格交集。

恢复策略：

- 旧 attempt 的原生 CopyKAT prediction 只有在 cell ID 行数相等、零缺失、零额外、零重复并通过 hash/覆盖审计后才复用。
- 失败与中断 attempt 原样保留；新 attempt 不覆盖旧证据。

### Stage 02：合并、映射与 CoVarNet 导出

- 对 1,849,413 个 active cells 合并逐样本判定，并保持 source/global position 对齐。
- 生成严格恶性状态、scMalignantFinder 辅助状态、UMAP、逐细胞/逐样本汇总和精简审阅对象。
- 硬门槛包括零未映射、零冲突、排除 GSE 不得泄漏。
- 导出 CoVarNet `meta.tsv` 与 `sample_metadata.tsv`。

### Stage 03–05：CoVarNet 与模块网络

| Parameter | Value |
| --- | --- |
| min cells/sample | 100 |
| abundance normalization | min-max |
| correlation | Pearson |
| target K | 9 |
| rank survey | 2–20；`vp` |
| final NMF | nsNMF；`nrun=30`；options `v` |
| top states | 10；模块解释 top 20 |
| network correlation threshold | 0.2 |
| FDR threshold | 0.05 |

Stage 04 输出全局和逐样本 Top CM 与 top-20 细胞状态解释；Stage 05 基于已完成网络表渲染全局和模块图，不重新拟合 NMF。

### Stage 06：LIANA sample × CM 通讯

- 使用 LIANA `rank_aggregate` 共识结果。
- 每个 CM 使用 top 10 members。
- sender/receiver group 最少 20 cells。
- expression proportion threshold 0.10。
- 逐样本、逐 CM 运行，输出完整 coverage table、Top LR、receiver marker/background、聚合 source-target heatmap 和逐样本 Top-LR heatmap。
- 失败 attempt 隔离保存，只有全部计算和合同检查通过后才发布 canonical result。

### Stage 07：NicheNet 验证

- 只使用本地 NicheNet v2 `ligand_target_matrix_nsga2r_final.rds` 和本地 nichenetr，不依赖运行中联网安装。
- 输入为 LIANA Top LR 与 receiver expression/signature background。
- 输出 ligand activity、candidate ligand、target links、跳过原因和热图。
- 零 ligand activity 或零 ligand-target links 视为失败，不得发布空结果冒充成功。

### Stage 08：最终审计

最终硬门槛覆盖：

- 注释与 active/excluded 分母；
- GSE166555 生物学样本重组；
- CNV 逐样本 terminal status、零映射冲突；
- CopyKAT/inferCNV 严格交集和 scMalignantFinder 辅助汇总；
- CoVarNet W/H、frequency、network、CM abundance、RDS、Top CM；
- 完整 sample × CM LIANA 组合与两类热图；
- NicheNet key/order 一致性和非空证据。

只有 `${CM_RUN_ROOT}/FINAL_AUDIT.json` 通过后，本轮 v2 才能标记完成。

## 2026-08-06 05:29 历史进度快照

### 已观察事实

- active v2 逐样本队列目标为 721 个生物学样本。
- 当前顺序已完成前 316/721 个样本，即 43.8%；第 317 个样本正在运行 CopyKAT。
- 当前样本为 GSE178341 的一个 colorectal cancer 生物学样本，2,925 cells、46 candidate、2,872 reference；inferCNV 已成功，阈值约 0.0782，CopyKAT R 进程持续计算。
- 主 runner 与 Python worker 已连续运行约 24 小时，日志仍按 60 秒更新，未观察到僵尸或重复锁持有者。
- 结果树可见 312 个规范 CopyKAT prediction 文件；少于顺序完成数的差异来自允许的跳过/历史 attempt 结构，不能用单一文件计数代替最终 terminal audit。
- Stage 02–08 的当前 canonical expected artifacts 均不存在，因此 active runner 在 Stage 01 完成后会重算而不是复用旧分母产物。
- 全链路 `FINAL_AUDIT.json` 尚未生成。

### 历史记录解释

`status/stage_history.tsv` 保留了旧清单上的早期 Stage 01–06 成功和一次 NicheNet 失败。由于当前 active manifest 已变为 v2 的 721 个生物学样本，旧阶段只能作为历史尝试，不能作为本轮完成证据。当前 canonical expected artifacts 已不在结果树中，后台入口会从 active Stage 01 继续顺序重建。

## 已知失败、修复和风险

- 早期 NicheNet attempt 因依赖/非空合同失败，后续改为本地模型并增加非空发布门槛。
- 长时间无主日志更新的根因曾是同步等待 CopyKAT 与 R 输出缓冲；当前已加入 60 秒 heartbeat，并能区分持续计算与挂起。
- GSE166555 原先按 2 个技术来源分组会破坏生物学样本单位；v2 已重组为 25 个样本，细胞覆盖保持不变。
- 当前运行速度受单样本细胞数及 CopyKAT baseline.GMM 路径影响，样本耗时差异很大，不给出固定 ETA。
- `genome=hg20` 是实际执行参数，但命名异常；应在最终生物学解释和发布前独立确认。

## 输出与证据

| Artifact | Logical path | Snapshot status |
| --- | --- | --- |
| active manifest | `${CM_RUN_ROOT}/manifests/samples_analysis_v2.tsv` | 721 samples / 1,849,413 cells |
| regroup audit | `${CM_RUN_ROOT}/manifests/GSE166555_SAMPLE_REGROUP_AUDIT.json` | pass |
| exclusion audit | `${CM_RUN_ROOT}/manifests/GSE_EXCLUSION_AUDIT.json` | pass |
| Stage 01 log | `${CM_RUN_ROOT}/logs/stages/01_sample_cnv.log` | active |
| stage history | `${CM_RUN_ROOT}/status/stage_history.tsv` | append-only |
| merged malignancy audit | `${CM_RUN_ROOT}/results/MALIGNANCY_MERGE_AUDIT.json` | pending current v2 |
| CoVarNet record | `${CM_RUN_ROOT}/results/covarnet_r/Covarnet_R_malignancy_intersection_v1/RUN_RECORD.md` | pending current v2 |
| LIANA record | `${CM_RUN_ROOT}/results/liana_cm/LIANA_RUN_RECORD.json` | pending current v2 |
| NicheNet record | `${CM_RUN_ROOT}/results/nichenet/NICHENET_RUN_RECORD.md` | pending current v2 |
| final audit | `${CM_RUN_ROOT}/FINAL_AUDIT.json` | pending |

## 2026-08-06 后续计划（历史）

1. 继续完成剩余 405 个 v2 样本的 Stage 01 terminal results。
2. 自动执行 Stage 02 合并并核验零未映射、零冲突和排除/重组合同。
3. 重算 CoVarNet K9、LIANA sample × CM 和本地 NicheNet 验证。
4. 只有最终审计通过后更新为 completed；若 `hg20` 不是预期值，应先形成参数纠正决定并使用新 attempt 重跑受影响阶段。

## 2026-08-08 06:16 审计快照

- active v2 分母保持 721 个生物学样本、1,849,413 个细胞；没有用旧 attempt 改写当前分母。
- Stage 01 已顺序完成 367/721（50.9%），尚余 354 个未完成。第 368 个为 GSE183904 的 GSM5573499/sample34：8,782 cells、669 candidates、7,562 references。
- 该样本 inferCNV 已成功，实际阈值 0.0425147545、window size 250；CopyKAT R 子进程仍活跃并按分钟写 heartbeat，快照前位于结果保存/热图阶段。主 runner、Python worker 和 R 子进程均存在。
- active v2 的 Stage 02–08 尚未开始，`FINAL_AUDIT.json` 仍不存在；旧 CoVarNet/LIANA/NicheNet 产物继续只作历史 attempt，不作为当前完成证据。
- 本轮结论仍是样本级 raw-count CNV 推断后再合并映射；不得把整合对象上的全局 CNV 解释为样本级恶性证据。
- 下一步：完成剩余 354 个 Stage 01 terminal results，再执行零未映射/零冲突合并、CoVarNet K9、LIANA sample × CM、NicheNet 和最终审计；发布前仍需确认 CopyKAT `genome=hg20`。

## 2026-08-12 10:06 审计快照

- active v2 Stage 01 已顺序完成 461/721 个生物学样本，即 63.9%；尚余 260 个样本没有完成。
- 第 462 个样本为 GSE299340 的 GSM9037562_297：8,734 cells、4,123 candidates、3,647 references。inferCNV 已成功，实际 threshold=0.109578741、window size=250。
- 该样本的 CopyKAT 只形成约 1.106 GB 的 raw gene-by-cell 中间矩阵；未形成 prediction、summary 或 terminal marker，因此不得把它计入 461 个已完成样本。
- 最后可见 heartbeat 停在 2026-08-10 08:02；2026-08-12 快照时主 runner、Python worker 和 CopyKAT R 子进程均不存在。主机之后发生过重启，但现有证据不能证明重启就是这次中断的原始原因。
- 当前应记为 `interrupted_stalled_nonterminal_stage01`，而不是活跃运行或成功完成。Stage 02–08 尚未开始，`FINAL_AUDIT.json` 仍不存在。
- 恢复时应保留第 462 个样本的非终态 attempt，以新的 append-only attempt 从该样本重新开始；在进入 Stage 02 前重新核对 461 个 terminal 结果的 cell-ID 覆盖和 hash，并继续确认 CopyKAT `genome=hg20` 的参数意图。
