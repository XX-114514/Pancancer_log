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

## 与 v2 诊断文字的差异（解释，待用户确认）

1. 子簇33、18（v2 文字中的"中性粒和鳞状上皮"双细胞）：CSF3R、FCGR3B、CXCR2 的检出率 ≤2%，S100A8/A9 来自鳞状上皮本身。不支持中性粒双细胞。
2. 子簇10（v2 文字中的"乳腺基底上皮"）：SAA1、CCL28、LTF、MMP7 上调，TP63 检出 1%，ACTA2 15%，更符合管腔祖细胞。
3. 代码与 v1 方案的偏差：`subcluster_v1.py` 的双细胞判定用 z 分加差异基因，没有实现方案中"超过30%细胞同时表达两组 marker"的条件；`residual_tumour_suspect` 只用 uncertain 占比和 CNV/scATOMIC 阈值，64/74 被标记，特异性低。

## 未验证

- 未使用细胞级共表达检验双细胞；子簇65、18 的双细胞解释未确认。
- 恶性判定未更改；"tumour_like"是基于聚合证据的解释，不是恶性判定。
- 未与 v2 逐细胞标签做子簇交叉表。
