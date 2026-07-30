# Project_v3 统一泛癌单细胞工作流

- 状态：按权威运行审计更新
- 最近证据日期：2026-07-31
- 权威运行：`20260725_120105_unified_pancancer_covarnet`
- 正式记录：[运行与参数审计](../runs/20260725_120105_unified_pancancer_covarnet.md)

## 目的

从来源和矩阵语义不完全一致的多数据集对象出发，建立分支明确、逐数据集可审计、能够扩展到全量细胞的 human discovery 对象；随后进行多方法整合、证据化注释、CoVarNet 模块发现和冻结模块的外部验证。

## 1. 库存、矩阵语义与样本键

1. 审计对象来源、物种、矩阵类型、样本键和用途分支。
2. 主计划包含 34 个 discovery datasets、8,887,638 个 QC 前细胞和 806 个 sample units；`GSE162498` 因 6,635,520 个细胞单列 extreme-scale planning-only。
3. 逐对象确认 `.layers["counts"]`。`GSE131907`、`GSE179994`、`GSE188737`、`GSE222315`、`GSE269826`、`GSE278694` 从已判定为 count-like 的 `X` 复制 counts，其余对象使用已有 counts layer。
4. human-only 基因空间检查排除使用小鼠基因标识的 `GSE211602`。

## 2. 逐数据集 QC

每个数据集调用统一 worker，并以已确认样本键在样本内计算阈值和 doublet。

固定参数：

| 参数 | 值 |
| --- | --- |
| floor `min_genes` / `min_counts` | 200 / 500 |
| strict `min_genes` / `min_counts` | 500 / 1,000 |
| MAD 倍数 | 3 |
| 线粒体比例阈值 | `min(25%, max(15%, median + 3 × MAD))` |
| Scrublet expected doublet rate | 0.04 |
| 单样本最大 doublet 移除比例 | 0.08 |
| Scrublet 最小样本细胞数 | 500 |
| Scrublet `min_counts` / `min_cells` | 2 / 3 |
| Scrublet variability percentile | 85 |
| Scrublet PCs | 30 |
| 随机种子 | 0 |

若一个数据集没有可用线粒体基因信号，则跳过 mt% 过滤并记录。代码实现中 strict profile 的选择条件是样本细胞数是否达到 1,000；这与“按测序深度选择 strict”的自然语言理解不同，复用时必须保留这一实现事实。

结果：33/33 个数据集通过，2,252,118 → 2,062,230 个细胞，移除 189,888 个细胞；695 个满足条件的样本全部成功运行 Scrublet。

## 3. 基因统一、归档分片与最终拼接

1. human-only 计划为 32 个数据集、2,040,048 个细胞、794 个样本；39/39 archival shards 可读，共同标准化基因为 11,063。
2. 最终拼接使用 38 个准备完成的 shards，生成 1,963,745 × 11,063 的 sparse Zarr，覆盖 31 个 GSE 和 773 个样本。
3. `GSE278694` 在计划对象到最终拼接之间被过滤。现有 notebook 和日志能证明过滤发生，但未发现科学理由，因此这一范围改变尚未完成治理闭环。
4. 全量拼接对象保留 counts layer、标准化的 cell ID 和必要 provenance columns。

## 4. 全量多方法整合

配置标识：`20260727_full_hvg4000_multimethod_v1`。

| 模块 | 固定参数 |
| --- | --- |
| 归一化 | library size 10,000，随后 `log1p` |
| 技术基因排除 | mt、ribo、Hb、ERCC、`MALAT1` |
| HVG | 4,000；候选 6,500；`seurat_v3` |
| PCA | 50 PCs |
| batch/sample key | `gse_id` / `sampleID` |
| 全局随机种子 | 20260727 |
| Harmony | 最大 20 iterations |
| neighbors | 30 |
| BBKNN | 每 batch 2 neighbors；trim 200；50 PCs |
| scVI | latent 30；2 layers；hidden 128；30 epochs；batch size 2,048 |
| CellTypist | `Immune_All_Low.pkl`；`majority_voting=false` |
| scIB subset | 最多 100,000 个细胞，按 GSE × CellTypist 分层确定性抽样 |

主要输出：

- Harmony graph：91,983,248 个非零连接。
- BBKNN graph：230,498,834 个非零连接。
- scVI graph：87,695,572 个非零连接。
- CellTypist：38 个 shards、97 个参考标签，confidence median 0.4824。

scIB 0.5.6 的 100k 子集缩放总分为 Harmony 0.6156、BBKNN 0.2000、scVI 0.7959。全图 Harmony/BBKNN/scVI 的 iLISI 分别为 0.01557/0.21751/0.00205，cLISI 为 0.99164/0.98741/0.99272，graph connectivity 为 0.42150/0.25459/0.46284。评价维度并不一致，因此不能宣称某一方法在所有方面最佳。证据注释分支选择 BBKNN 图。

## 5. 证据化注释

输入为全量 BBKNN 对象。一级 Leiden resolution 为 0.1；17 个 lineage 的二级 resolution 分别为：

| lineage | resolution | lineage | resolution |
| --- | ---: | --- | ---: |
| Tcell | 1.0 | NK | 0.8 |
| ILC | 0.6 | Bcell | 0.8 |
| Plasma | 0.6 | Monocyte | 0.8 |
| Macrophage | 0.9 | Dendritic | 0.8 |
| Granulocyte | 0.6 | Mast | 0.6 |
| Epithelial | 1.0 | Endothelial | 0.8 |
| Fibroblast | 0.9 | Erythroid | 0.5 |
| Platelet | 0.4 | HSPC | 0.7 |
| Neural | 0.6 |  |  |

种子为 20260729。小于 30 个细胞的 cluster 以 50-PC PCA centroid 合并到稳定 cluster；Neural 因子图缺少稳定 cluster，采用 MiniBatchKMeans fallback。

表达证据从 raw-count concat 流式计算：每块 10,000 个细胞，library size 10,000 后 `log1p`，marker UMAP 最多抽样 180,000 个细胞，marker panel 194 个基因。v2 对 15 个 clusters、404,773 个细胞进行污染标签细化。最终得到 17 个 major lineages、136 个 lineage clusters、106 个 final annotations；交付审计为 `PASS`。仍有 56 个 low-confidence、10 个 unresolved 和 5 个 ambiguous clusters。

## 6. CoVarNet

输入映射为 `sample_uid → sampleID`、`major_celltype → major`、`final_annotation → sub`。

| 参数 | 值 |
| --- | --- |
| 样本最少细胞数 | 100 |
| abundance normalization | min-max |
| correlation | Pearson |
| 目标 K | 9 |
| rank survey | 2–20，`vp` |
| 最终 NMF | nsNMF，`nrun=30`，options `v` |
| 每模块展示 top states | 10；解释报告使用 top 20 |
| 网络相关阈值 | 0.2 |
| FDR 阈值 | 0.05 |

输入 1,963,745 个细胞、773 个样本；过滤 18 个不足 100 cells 的样本后保留 1,962,479 个细胞、755 个样本。最终固定 `W` 为 106 × 9，`H` 为 9 × 755；网络为 101 nodes / 374 edges。

R 端 `sc_cm_recover` 因类型比较错误失败后，以 `coef(nmf_final)` 恢复模块 abundance。R 网络绘图因 `Set3` 调色板最多 12 色而无法覆盖 17 个 major classes，随后用 Python/networkx 重新生成；计算表和模块权重未重新拟合。

## 7. 固定模块外部验证

三个队列均使用同一个冻结的 106-state hierarchical mapper 和 106 × 9 `W`，不在验证队列重新训练。表达处理为 count per 10,000 + `log1p`，cell prediction chunk 为 5,000，confidence threshold 为 0.50。样本层以 soft frequency 为主、hard label 为敏感性分析；发现队列使用 min-max scaling，模块投影使用非负最小二乘。

- `GSE195832`：30,190 个细胞，8 个样本，4 个 nivolumab pre/post 配对；9/9 个模块的 soft/hard median 方向一致，median rho 0.9940。CM06 在 4/4 配对中上升，但 `q=0.641`。
- `GSE123813`：79,040 个细胞，30 个生物样本，15 个完整配对；主分析限于 7 个 sort-matched 配对，15 配对仅作敏感性分析。median soft-hard rho 0.9535，NNLS residual 0.324。
- `GSE169246`：489,490 个细胞，78 个样本、22 位患者；主分析为 20 个 blood pre/post pairs，另有 7 个 blood triples 和 11 个 same-site tumor pairs。mapper gene coverage 3,999/4,000，median confidence 0.6906，median soft-hard rho 0.9362，NNLS residual 0.2992。

三个审计均通过，但 BH 校正后没有显著模块。现阶段只能报告探索性方向，不能写成疗效、生存或复发结论。

## 已知限制

- `GSE162498` 尚未进入本次最终整合。
- `GSE278694` 的移除理由未记录。
- 参考标签和自动 marker 证据不能替代低置信度 cluster 的人工复核。
- 恶性/CNV 分支未在本次审计中确认完成。
- 源阶段状态 manifest 滞后；正式 run record 和最终 audit 文件是本次状态判断的依据。

## 外部证据

- `${RUN_ROOT}/config/full_hvg4000_multimethod_v1.json`
- `${RUN_ROOT}/objects/02_integrate/archival_shards_concat.zarr`
- `${RUN_ROOT}/objects/03_preprocess/20260727_full_hvg4000_multimethod_v1/`
- `${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/FINAL_DELIVERY_AUDIT.json`
- `${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/covarnet_r/Covarnet_R_refined_v2_20260729_072200/`
- `${RUN_ROOT}/09_immunotherapy_validation/*/FINAL_AUDIT.json`
