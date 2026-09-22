# 投稿用注释上游冻结决定

## 冻结结论

从 2026-09-22 起，本投稿准备分支使用以下不可混用的版本轴：

| 轴 | 投稿权威版本 | 状态 | 可用于何处 |
| --- | --- | --- | --- |
| 全局细胞身份 | V7 | `FROZEN_PROJECT_REFERENCE` | 主文分母、Methods、主/补充图表 |
| CNV 恶性证据 | V7 | `FROZEN_PROJECT_REFERENCE` | 主文和补充分析，skip 必须保持不可评估 |
| V8 谱系分类学 | 五个 lineage freeze | `FROZEN_TAXONOMY_ONLY` | 方法发展、Extended Data 或明确标注的开发证据 |
| V8 全量 metadata candidate | V8 candidate on V7 | `NOT_FROZEN` | 软件合同/审计说明，不作新生物学真值 |
| V9 全量注释候选 | V9 S0–S7 candidate | `FROZEN_CANDIDATE` | 候选/方法开发，不进入正式分母或下游 |
| 通讯/空间/临床证据 | V5 CoVarNet/LIANA snapshot | completed snapshot | 仅以 V5-derived 身份报告 |

## 冻结的全局 V7 合同

- 4,676,787 cells、1,322 samples、54 GSE、43 cancers。
- Cell-index SHA-256：
  `3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16`。
- 主要字段：`major_celltype_v7`、`final_annotation_v7`、
  `annotation_confidence_v7`、`malignancy_call_v7`、
  `malignancy_tier_v7`、`final_annotation_with_malignancy_v7`。
- 严格分析必须使用 release 中定义的 major、subtype 和 malignant eligibility
  字段；不得把 skipped CNV 样本改记为 negative。

## 冻结的 V8 lineage-only 决定

- Myeloid/DC：接受 FCN1 monocyte、Macrophage、cDC2 identity；接受 C1QC 和
  inflammatory state；其余列明状态退回父类。
- T/NK：冻结保守 lineage taxonomy 和 tissue-resident-like state 证据；没有精细
  identity 通过跨 GSE recurrence。
- Fibroblast：只冻结 Fibroblast parent 与 Fibroblast_myCAF；其他状态为父类/overlay。
- Endothelial：冻结 Endothelial parent、arterial 和 venous；capillary、
  angiogenic-tip、lymphatic 退回父类/overlay。
- B/Plasma：冻结 Bcell、Plasma 和 Bcell_naive；Plasma-IgA 退回 Plasma；
  IFN/cycling 为 overlay；held-out predictions 仍封存。

这些决定不得被描述为 `FROZEN_V8_BIOLOGICAL_REFERENCE`，也不得拼接成未经审计的
全局细胞级对象。

## V9 候选边界

V9 保持与 V7 相同的细胞宇宙，但其 L3 人工复核、低稳定性门槛、CNV 角色变化和
632 样本重跑建议尚未闭环。S6 主要继承 V7 CNV，不是 V9 身份下重新推断的恶性轴。
因此本次仅冻结其“候选状态和证据哈希”，不冻结其生物学标签为投稿真值。

## 变更控制

允许在不改变科学结论的情况下继续：语言润色、引用核验、图版布局、颜色与字体统一、
补充材料编排。下列变化必须创建新的 append-only release/run，不能覆盖本快照：

- 改变细胞宇宙、基因映射、marker/接受规则或 candidate/reference；
- 改变 CNV 参数、恶性判定、分母、统计模型或排除标准；
- 解封 held-out 数据后再调参；
- 将 V5 downstream 重新标记为 V7/V8 derived；
- 将 lineage-only freeze 晋升为全局 V8。

## 尚未关闭的科学门槛

1. 决定 V8 sealed benchmark 与 V9 候选中哪一条成为唯一升级路径。
2. 若采用 V9，完成 15 个 L3 cluster 人工复核、稳定性处置和 632 样本 CNV 决策。
3. 冻结统一 identity，核验 Copykat_python `genome='hg20'` 后只重建一次 malignancy。
4. 完成一次最终 robustness pass，再决定是否建立新的全局 biological reference。
5. 在最终身份版本上重跑 CoVarNet/LIANA，或把 V5-derived 通讯证据降为独立历史分析。
6. 补齐 GSE278694 排除依据和低置信/争议注释的发布处置。
