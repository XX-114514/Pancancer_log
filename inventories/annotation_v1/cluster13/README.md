# 全局簇13子聚类：人工复核草案（job 30667416）

**状态：诊断性草案，不是最终注释，不进入 v2 流程。** [v2 方案](../../../methods/tme-annotation-v2-lodi-workflow-20261105.md)已规定簇13的细胞回到各癌种分析，本运行只作诊断；CoVarNet v1 细胞范围以 v2 方案中 2026-11-06 的用户决定为准。

运行：Slurm job 30667416，COMPLETED，exit 0，41分47秒（2026-11-05）。全局 Leiden 簇13共290,546个细胞，分辨率1下74个子簇，其中非恶性68,374个。

| 文件 | 内容 |
|---|---|
| `manual_labels_draft.tsv` | 74个子簇逐个人工标签、分类、置信度、证据 |
| `residual_tumour_suspects_for_user.tsv` | 自动标记 `residual_tumour_suspect` 的64个子簇及人工分层，待用户决定 |
| `subclusters.tsv`、`panel_scores.tsv`、`de_top50.tsv`、`gene_dot.tsv`、`gene_measured_by_source.tsv.gz` | 运行输出的聚合表（无细胞级数据） |

生成脚本：[cluster13_manual_labels_draft_v1.py](../../../scripts/annotation_v1/cluster13_manual_labels_draft_v1.py)。

## 人工分类（已观察计数）

| 分类 | 子簇 | 细胞 | 非恶性细胞 |
|---|---:|---:|---:|
| program_mixed（应激、管家、缺氧、增殖等程序，无身份） | 13 | 209,571 | 47,852 |
| organ_epithelium | 11 | 25,615 | 12,756 |
| tumour_like | 39 | 47,783 | 6,080 |
| ambient_or_doublet | 5 | 2,742 | 796 |
| low_quality | 4 | 4,618 | 686 |
| non_epithelial_parenchyma（间皮、RPE） | 2 | 217 | 204 |

## 与 v2 诊断文字的差异（解释；第1、2项已于2026-11-08按用户要求更正到 v2 方案）

1. 子簇33、18（v2 文字中的"中性粒和鳞状上皮"双细胞）：CSF3R、FCGR3B、CXCR2 的检出率 ≤2%，S100A8/A9 来自鳞状上皮本身。不支持中性粒双细胞。
2. 子簇10（v2 文字中的"乳腺基底上皮"）：SAA1、CCL28、LTF、MMP7 上调，TP63 检出 1%，ACTA2 15%，更符合管腔祖细胞。
3. 代码与 v1 方案的偏差：`subcluster_v1.py` 的双细胞判定用 z 分加差异基因，没有实现方案中"超过30%细胞同时表达两组 marker"的条件；`residual_tumour_suspect` 只用 uncertain 占比和 CNV/scATOMIC 阈值，64/74 被标记，特异性低。

## 未验证

- 未使用细胞级共表达检验双细胞；子簇65、18 的双细胞解释未确认。
- 恶性判定未更改；"tumour_like"是基于聚合证据的解释，不是恶性判定。
- 未与 v2 逐细胞标签做子簇交叉表。

## 与 v2 标签的交叉表（2026-11-09）

68,374个非恶性细胞按 cell_id 与 `cell_labels_r3.tsv.gz` 一一连接，无缺失。文件：`c13_subcluster_by_v2_fate.tsv`、`c13_subcluster_by_v2_label.tsv`、`c13_retained_lineage_evidence.tsv`。

| 人工分层 | 实质 | 剔除或搁置 | 保留：标记 | 保留：低深度 | 保留：主分析样 | 合计 |
|---|---:|---:|---:|---:|---:|---:|
| likely_residual_tumour | 5,734 | 322 | 5 | 38 | 25 | 6,124 |
| identity_with_flag | 3,464 | 11 | 2 | 2 | 3 | 3,482 |
| program_no_identity | 33,484 | 6,654 | 69 | 4,551 | 3,094 | 47,852 |
| technical_not_a_tumour_call | 410 | 18 | 0 | 4 | 0 | 432 |
| not_flagged | 9,880 | 386 | 0 | 164 | 54 | 10,484 |
| 合计 | 52,972 | 7,391 | 76 | 4,759 | 3,176 | 68,374 |

"保留"指 v2 中仍在 T/NK、髓系、B/浆、基质、内皮或小胶质池（`r3_*`），共8,011个。"主分析样"只按 r1/r2 action 近似，不等于旧服务器工作版的 `include_main_subtype`。

- 已观察：8,011个保留细胞中，likely_residual_tumour 只有68个（主分析样25个）。145个（1.8%）至少检出 EPCAM、KRT8、KRT18、KRT19 中两个、又没有对应谱系 marker，其中137个来自 program_no_identity，没有一个来自 likely_residual_tumour。
- 解释：簇13经疑似残留肿瘤子簇进入 CoVarNet 候选的细胞量可忽略。
- 未验证：低深度细胞本来 marker 检出就少；CNV 中位数在各分层差异不大，未作单细胞恶性复判。

未改动任何标签。
