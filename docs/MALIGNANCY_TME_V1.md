# 恶性细胞判定 v1 与 TME 整合 v1

构建：2026-11-02（UTC）起。代码：`scripts/malignancy_v1/`。输入：[整合数据集 v2](INTEGRATION_V2.md)的72个来源分片（4,720,931个细胞）。组成与一致性表：[inventories/tme_v1/](../inventories/tme_v1/)。

## 为什么要做

后续用 CoVarNet 识别 CM 模块时只用非恶性（TME）细胞。v2 的 UMAP 把恶性细胞和 TME 一起整合，问题如下：
- 上皮、黑色素细胞、胶质细胞区域按来源分成小岛。同来源近邻比例的中位数：上皮0.43、黑色素0.83、胶质0.67。随机混合的期望约0.02–0.04。
- 19%的细胞没有明确的区室。
- 只按来源去批次。
- 没有定量的整合指标。
- 调色板分辨不清，点相互遮挡。

已有的恶性判定口径不统一：core 用 scATOMIC + v7 inferCNV，immune/public 两个 store 只有 infercnvpy，新来源没有判定。因此全部72个来源用同一流程重跑。

## 恶性细胞判定（每个来源一个 Slurm 作业，按单位/样本运行）

| 步骤 | 内容 |
|---|---|
| 快速聚类 | 每个单位用原始计数：normalize 1e4 → log1p → 2,000 HVG → 30 PC → 15近邻 → Leiden 1.0 |
| 谱系 | 按聚类 marker 均值判定。免疫：PTPRC。基质：取成纤维（COL1A1/COL1A2/DCN）、内皮（PECAM1/VWF）、周细胞（RGS5）三项中的最大值。其余为候选。DLBCL 的 B 细胞聚类算候选，不作参考 |
| scATOMIC 2.0.3 | 单位内全部细胞；`run_scATOMIC` + `create_summary_matrix(use_CNVs=FALSE, modify_results=TRUE, min_prop=0.5)`；使用 MAGIC 插补（默认）。scATOMIC 判为肿瘤≥50%的免疫/基质聚类改成候选（`candidate_from_*`），不再作为 CNV 参考 |
| inferCNV 1.18.1 参考 | 参考池为本单位剩下的免疫/基质细胞。不足200个时，从同一来源的其他单位借免疫/基质细胞，补足到400个，记为`ref_borrowed` |
| inferCNV 1.18.1 对照与参数 | 参考池的1/4（最多500个）作为对照细胞参与打分，其余（最多1,500个）作参考；候选细胞为观测。cutoff 0.1，denoise，HMM 关，`analysis_mode="samples"`，基因位置用 hg38 GENCODE v27 |
| CNV 指标 | score：去噪后 CNV 谱偏离1的均方。corr：与本单位 score 前5%观测细胞的平均 CNV 谱的相关（Tirosh 2016 / Neftel 2019） |
| CNV 判定 | 两个阈值都取对照细胞的第99百分位。候选聚类的 score 中位数或 corr 中位数超过阈值即为 CNV 阳性（按聚类判定） |
| 共识 | 免疫/基质谱系 → non_malignant。候选且 scATOMIC 判肿瘤、CNV 阳性 → malignant。候选且 scATOMIC 判正常、CNV 阴性 → non_malignant。两者不一致或 CNV 无法评估 → **uncertain：按用户决定留在 TME 中，并以 `tme_flag` 标记** |

开发过程中修正的问题（均在正式运行前完成；测试结果移到 `malignancy_v1/test_runs/`，未使用）：
1. 被 scATOMIC 判为肿瘤的"基质"聚类曾被当作 CNV 参考。现在先跑 scATOMIC，这些聚类改为候选。
2. 以肿瘤为主的单位几乎没有免疫/基质细胞，整个单位都会变成 uncertain，例如 GSE161277 Patient0 和葡萄膜黑色素瘤。现在从同一来源借参考细胞。
3. 葡萄膜黑色素瘤（GSE138665）中，scATOMIC 判为肿瘤的7,129个细胞大多因均方 score 不足被判为 CNV 阴性。加入 corr 指标后，6,639个判为恶性，536个为 uncertain。

## TME 整合 v1

细胞：v2 中 `malignancy_call != malignant` 的细胞（包括带标记的 uncertain）。配方沿用 v2，并做以下改动：
- HVG 只用 TME 细胞重新计算。
- Harmony 用两个批次变量：来源（theta 2）和单位/样本（theta 1）。
- Leiden 1.0 聚类后按 marker 集 z 分给出粗区室。
- 在固定的20万细胞子样本上，计算 Harmony 前后的 iLISI（来源、单位）、cLISI（区室）和同来源近邻比例。
- UMAP 设 min_dist 0.3。图件包括分面图（每个癌种/来源单独高亮）、密度图，并随机打乱绘制顺序。

结果见下文（运行完成后填写）。
