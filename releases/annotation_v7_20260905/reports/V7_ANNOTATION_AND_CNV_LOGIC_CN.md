# V7 注释与 CNV 冻结逻辑

- 细胞宇宙：4,676,787 个细胞，索引哈希 `3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16`。
- 身份主轴：V2 保守大类与亚型；V5 全部字段保留为历史与审计父版本。
- CNV 输入：5M 合并对象中的逐样本原始非负整数 counts；不使用 Harmony/scVI/整合表达矩阵。
- 候选：按癌种指定肿瘤谱系；Unknown/Ambiguous/边界细胞仅在 operational 候选落入该癌种肿瘤谱系时纳入。
- 参考：V2 strict major，且排除候选、QC、doublet_likely 与 ambient_or_doublet_review。
- 严格恶性：inferCNVpy 阳性与固定版本 Copykat_python 非整倍体的交集。
- 主恶性：逐样本 CNV cluster 规则；并集只作为敏感证据层。
- `cnv_no_malignancy_detected` 仅表示未检出 CNV 恶性证据，不等于证明良性。
- CopyKAT 身份：本次使用 `Copykat_python 1.0.0 @ c2390e8a46e2`，不是 R CopyKAT。
- 跳过样本：见 `tables/skipped_samples_with_explanations.tsv`，绝不把跳过强制改为正常。
