# PanCancer 5M V7 冻结参考与注释逻辑

## 1. 版本结论

V7 是本项目当前冻结的细胞身份与恶性证据参考版本。它覆盖 4,676,787 个细胞、
1,322 个原始样本、54 个 GSE 和 43 个癌种。V7 不覆盖或删除任何 V5/V6 字段，
而是在相同且唯一的细胞索引上完成两项更新：

1. 身份轴采用 2026-09-05 不确定性修复后的保守 V2 大类和亚型；
2. 恶性轴依据该新版大类，按原始样本从 raw integer counts 重新建立
   candidate/reference，并重跑 inferCNVpy 与固定版本 Copykat_python。

冻结标识是 `FROZEN_V7.json`，schema 为
`pancancer_annotation_v7_cnv_rerun_20260905`。cell index SHA-256 为
`3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16`。

从 2026-09-06 起，新启动的项目分析应默认引用 V7。旧 V5/V6、旧 CNV、旧
CoVarNet/LIANA 产物仍可用于历史比较，但必须明确版本，不能无说明地称为 V7。

## 2. 冻结对象与不可变性

权威数据产物：

| 产物 | 用途 | SHA-256 |
|---|---|---|
| `checkpoints/final_annotation_v7.pkl.gz` | 4,676,787 × 128 的完整细胞级注释和历史字段 | `e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e` |
| `objects/final_annotation_v7_obs_only.h5ad` | 4,676,787 × 0 的轻量 AnnData；obs 含 V7 关键字段 | `e5b68d7e5f127bcdcdadc3aaa94a9bce027ec032e6b82e68e5847ea5ccf7b38e` |
| `FROZEN_V7.json` | 细胞宇宙、统计、路径和哈希的冻结清单 | 以文件内容为准 |

已重新计算两个主要对象的 SHA-256，和冻结清单一致。V7 采用 append-only
规则：不得覆盖上述文件；任何会改变细胞身份、candidate/reference、CNV 参数、
统计口径或恶性规则的更新都必须生成新 run，并命名为 V8 或更高版本。

## 3. 版本继承关系

```text
5M merge raw counts（4,676,787 × 10,432）
  ├─ V5/V6 历史最终表
  │    └─ 保留原精细注释、旧恶性结果、状态/QC 和下游字段
  └─ uncertainty-resolution V2
       ├─ 修复大类 Unknown/Ambiguous、来源标签重复投票和负 marker 逻辑
       ├─ 保留连续边界、provisional、doublet/ambient review
       └─ 扩充 panel + cluster FindMarkers 处理 subtype unresolved
            ↓
       V7 身份主轴

V7 身份主轴 + 按癌种/样本重建的 candidate/reference + raw counts
  ├─ per-sample inferCNVpy
  ├─ per-sample Copykat_python
  └─ 精确映射、cluster 主判定和交集/并集证据分层
       ↓
       frozen V7
```

身份和恶性是两个独立轴。`final_annotation_v7` 不能替代 CNV 字段；CNV 阳性也
不能自动证明组织病理身份。

## 4. 大类注释逻辑

### 4.1 证据来源

V7 的 `major_celltype_v7` 直接冻结自
`major_celltype_uncertainty_v2`。其父流程基于三类证据：

1. 来源数据真正提供的逐细胞标签；
2. 使用真实全转录组 library size 计算的 marker panel 分数；
3. CellTypist 预测，仅在其适用范围和足够置信时作为支持证据。

`majority_voting`、`predicted_labels`、`categorized_labels` 和
`over_clustering` 等旧自动结果归类为 `previous_automatic`，不再伪装成独立的
文献真值重复投票。只有源文件中明确的逐细胞列，如 `Cell_type.refined`、
`cell_type`、`celltype`、`gse116256_CellType` 和 `celltype_simplified3`，归类为
`source_provided_per_cell`；它们仍然是证据而不是不可推翻的金标准。

### 4.2 负 marker 与冲突

旧规则中单个负 marker 即可硬否决候选，容易把环境 RNA、共享状态基因、EMT
或轻度混合当成大类错误。V2/V7 改为：

- 仅 1 个负 marker 超阈值：`provisional_negative_review`，不再硬否决；
- 至少 2 个方向一致的负 marker：视为 coherent competing program；
- 跨 compartment 的一致冲突保留为 `Ambiguous`；
- 跨 compartment 冲突同时位于样本 doublet score 高分位时，设置
  `exclude_from_primary_major_v2=True`，但不删除该细胞。

本轮阈值为 marker top score ≥ 0.25、top1–top2 margin ≥ 0.08、CellTypist 支持
≥ 0.5、CellTypist 单独提名 ≥ 0.8、负基因表达阈值 0.5、coherent negative 至少
2 个基因、doublet 样本内分位数 0.95。

### 4.3 margin 与连续状态

top1–top2 接近不被一律解释成失败。直接生物邻居之间的低 margin 使用层级边界：

- `T_NK_ILC_boundary`
- `B_Plasma_boundary`
- `Myeloid_boundary`
- `Stromal_boundary`
- `Related_lineage_boundary`
- `Hematopoietic_differentiation_boundary`

这些细胞可用于连续谱、轨迹或联合重聚类，但不进入严格离散大类比较。跨
compartment 的近邻冲突继续保留为 `Ambiguous`，优先复核 doublet、ambient RNA、
吞噬程序和肿瘤状态。

### 4.4 大类统计口径

| 口径 | 细胞数 | 占全部细胞 |
|---|---:|---:|
| 严格大类可用 | 4,057,641 | 86.76% |
| 大类非严格：boundary、provisional、exact unresolved 等 | 619,146 | 13.24% |
| exact `Unknown` | 53,771 | 1.15% |
| exact `Ambiguous` | 76,072 | 1.63% |
| exact Unknown + Ambiguous | 129,843 | 2.78% |

因此，“大类不确定”至少有两个合法定义：exact Unknown/Ambiguous 是 2.78%；若
问题要求严格离散大类，则非严格比例是 13.24%。两者不能混用。

严格大类分析建议同时满足：

```text
strict_major_eligible_v7 == True
exclude_from_primary_major_v2 == False
```

## 5. 亚型与 unresolved 逻辑

V7 的 `final_annotation_v7` 冻结自
`final_annotation_uncertainty_v2`，不是重新强迫所有细胞得到成熟亚型。

父流程先在每个大类内部使用冻结的 lineage-specific Leiden 结构和扩展 marker
panel。panel 评分要求完整 panel 支持、canonical marker 锚定、top–second margin
以及跨样本复现；不能只依赖百万细胞 pooled p 值。

对已知大类内原有 1,073,743 个 subtype unresolved 细胞：

- 91,805 个通过扩展 panel、canonical marker 和跨样本复现，被严格救援；
- 609,412 个由同大类 cluster-vs-rest 全基因 FindMarkers 获得描述性
  `_unresolved_GENE` 后缀；
- 372,526 个仍无跨样本稳健 marker。

FindMarkers 使用整数 counts 流式计算 `log1p(CP10K)` 效应量；技术、线粒体、
核糖体、免疫球蛋白/TCR、应激和性别相关基因不得作为自动成熟亚型名称。默认
de novo 门槛包括 `pct_in ≥ 0.20`、`delta_pct ≥ 0.10`、`log2FC ≥ 0.50` 和至少
3 个样本复现。

`Fibroblast_unresolved_FUNG` 这类名称的定位是“cluster 的可追踪描述”，不是已
验证的新细胞类型。它可以用于绘图、再聚类和后续人工命名，但不能进入严格亚型
假设检验。

| 亚型口径 | 细胞数 | 占全部细胞 |
|---|---:|---:|
| 严格亚型可用 | 2,455,424 | 52.50% |
| 非严格亚型/描述性 unresolved | 2,221,363 | 47.50% |

严格亚型分析只使用：

```text
strict_downstream_subtype_eligible_v7 == True
```

## 6. 基因名与可检测性处理

身份父版本使用 `wide_marker_panel_full_library_masked_v2.h5ad`。关键原则是：

1. CP10K 分母为真实全转录组 library size，不是 marker 子矩阵总和；
2. 每个 GSE × gene 的结构性缺失由 availability mask 显式表示，不能把“未测得”
   当作“生物学零表达”；
3. `KIM1` 统一解析到 HGNC symbol `HAVCR1`；
4. 旧数据集的 legacy symbol 采用逐 GSE fallback，包括
   `SELENOP/SEPP1`、`ACKR1/DARC`、`JCHAIN/IGJ` 和 `CXCL8/IL8`；
5. CNV 使用修复后的 hg38 GENCODE v27 HGNC symbol 坐标文件，10,432/10,432 个
   counts 基因均能精确映射。

这些修复避免因基因别名或面板缺测把细胞错误标成 marker 阴性。

## 7. V7 candidate/reference 重建

### 7.1 总体规则

candidate/reference 按原始样本独立重建，不在 Harmony、scVI、BBKNN 或全局
归一化矩阵上做 CNV。

- candidate：该癌种预期肿瘤谱系、允许的 context 谱系，以及仅当
  `operational_major_candidate_v2` 落入该癌种允许谱系的 Unknown/Ambiguous/边界
  细胞；`doublet_likely` 硬排除。
- reference：必须是 V2 strict major、不是 candidate、`qc_annotation_final` 为空，
  且不属于 `doublet_likely` 或 `ambient_or_doublet_review`。
- 若样本的 primary reference 少于 50 个，才启用癌种策略定义的 strict fallback
  reference。

最终角色：1,156,120 个 candidate、2,465,967 个 reference、1,054,700 个其他
非候选细胞。39,573 个 candidate 来自不确定细胞的 operational 谱系；10,205 个
潜在 candidate 因 `doublet_likely` 被排除。

### 7.2 癌种特异策略

| 策略 | 癌种匹配 | candidate expected/context | primary reference | fallback |
|---|---|---|---|---|
| `epithelial_solid_default` | 其他实体瘤 | Epithelial / 无 | 免疫细胞 | Endothelial、Fibroblast |
| `cns_malignancy` | Glioblastoma | Neural / Epithelial | 免疫细胞 | Endothelial、Fibroblast |
| `melanocytic_malignancy` | melanoma，包括 Uveal/Acral | Melanocytic、Neural / Epithelial | 免疫细胞 | Endothelial、Fibroblast |
| `neuroendocrine_malignancy` | neuroendocrine、SCLC | Epithelial、Neural / 无 | 免疫细胞 | Endothelial、Fibroblast |
| `mesenchymal_malignancy` | Osteosarcoma | Fibroblast / Epithelial | 免疫细胞 | Endothelial |
| `mixed_mesenchymal_malignancy` | Mesothelioma、Carcinosarcoma | Epithelial、Fibroblast / 无 | 免疫细胞 | Endothelial |
| `myeloid_malignancy` | AML/leukemia | HSPC、Monocyte、Macrophage、Dendritic、Granulocyte、Erythroid / 无 | T/B/NK/ILC/Plasma/Mast | Endothelial、Fibroblast |
| `germ_cell_malignancy` | Seminoma | Epithelial、HSPC / 无 | 免疫细胞 | Endothelial、Fibroblast |

“免疫细胞”在上述实体瘤策略中指 Tcell、Bcell、NK、ILC、Plasma、Monocyte、
Macrophage、Dendritic、Mast 和 Granulocyte。

## 8. inferCNVpy 与 Copykat_python

### 8.1 输入和执行边界

- 输入：`5M_merge_concat.zarr` 中逐样本 raw non-negative integer counts；
- 每个样本单独运行，样本间不共享 CNV 阈值；
- inferCNVpy 至少需要 30 个 candidate 和 50 个 reference；reference 超过 5,000
  时用 sample UID 派生的固定随机种子确定性抽样；
- Copykat_python 至少需要 200 个 candidate，且其中至少 200 个细胞检测到超过
  200 个基因；
- 所有依赖、本地基因坐标和 Copykat_python 仓库均在预检中确认不依赖远程连接。

### 8.2 inferCNVpy

实际参数与判定：

- `layer='counts'`；window size = 250；
- 每样本 reference CNV burden 的 99% 分位数作为该样本阳性阈值；
- candidate 的 `cnv_burden >= sample threshold` 定义为单细胞 inferCNV 阳性；
- 在 candidate CNV space 内做 PCA、最多 30 邻居、Leiden resolution 0.15；
- 小于 100 个细胞的小 CNV cluster 在 profile correlation ≥ 0.95 时合并；
- random state 固定为 0，单个样本内部 inferCNV 线程固定为 1。

### 8.3 Copykat_python

本次实际使用：

```text
Copykat_python 1.0.0
commit c2390e8a46e21fb42aa06f05839900c0f12564c5
```

它不是 R 包 `copykat`，文档和图表中必须写作 `Copykat_python`。实际引擎参数包括
`id_type='S'`、`cell_line='no'`、`ngene_chr=5`、`min_gene_per_cell=200`、
`LOW_DR=0.05`、`UP_DR=0.1`、`win_size=25`、`KS_cut=0.1`、
`distance='euclidean'`、`genome='hg20'`、每样本 1 core。

### 8.4 运行完成度

| 方法 | success | 跳过 | failed |
|---|---:|---:|---:|
| inferCNVpy | 938 | 384：candidate/reference 不足 | 0 |
| Copykat_python | 668 | 653：细胞不足；1：`all_cells_filtered` | 0 |

跳过是显式不可评估，不是阴性。唯一算法数据限制样本为
`GSE161529|GSM4909319_mER-PM0178`；该样本 inferCNVpy 成功，Copykat_python
因 `all_cells_filtered` 跳过。

CNS 策略 44 个样本中 inferCNVpy 成功 38、Copykat_python 成功 33；黑色素谱系
策略 62 个样本中分别成功 27 和 18。其余主要因为单样本 candidate/reference
未达到方法门槛。

## 9. 恶性判定与统计口径

### 9.1 cluster 主判定

每个样本、每个合并后的 CNV cluster 独立判断。cluster 至少含 30 个细胞，并满足
以下任一条件时为 `malignant`：

```text
(inferCNV-positive fraction >= 0.60
 AND mean CNV burden >= 1.05 × sample reference threshold)
OR
(Copykat_python evaluable AND aneuploid fraction >= 0.50)
```

未达到 malignant 时，若 inferCNV 阳性比例 ≥ 0.20，或 Copykat_python 中存在任意
aneuploid，则为 `uncertain`；否则为 `cnv_no_malignancy_detected`。未运行成功的
candidate 为 `not_evaluable`；其他细胞为 `not_candidate`。

`malignancy_call_v7` 的全量计数：

| 值 | 细胞数 | 解释 |
|---|---:|---|
| `not_candidate` | 3,520,667 | 未进入该癌种肿瘤谱系候选池 |
| `malignant` | 515,031 | 逐样本 CNV cluster 主规则支持恶性 |
| `uncertain` | 484,919 | 有部分 CNV 证据，但 cluster 主规则不足 |
| `cnv_no_malignancy_detected` | 110,886 | 本方法未检出恶性 CNV 证据，不等于证明良性 |
| `not_evaluable` | 45,284 | candidate，但方法条件不满足 |

### 9.2 单细胞两方法证据层

只有 inferCNVpy 和 Copykat_python 均可评估的 914,868 个细胞才计算交集/并集：

- inferCNVpy 可评估 1,110,836，阳性 352,741；
- Copykat_python 可评估 946,882，aneuploid 427,707；
- 两方法交集阳性 238,418：`high_confidence_malignant`；
- 两方法并集阳性 513,509：敏感证据；
- 两方法一致阴性 401,359；
- 两方法不一致 275,091；
- 其余 241,252 个 candidate 缺少至少一种可评估方法。

### 9.3 `malignancy_tier_v7` 不是主判定的同义词

`malignancy_tier_v7` 按以下优先级把证据层压成互斥类别：

```text
high-confidence intersection
  > cluster-supported malignant
  > sensitive union only
  > uncertain / no detection / not evaluable / not candidate
```

因此不能把 `high_confidence_malignant` 与 `cluster_supported_malignant` 简单相加后
称为 `malignancy_call_v7 == malignant`。单细胞交集阳性可能位于一个未通过
cluster 主规则的 cluster；这是两条证据轴的定义差异，不是重复细胞或映射错误。
实际有 30,080 个两方法交集阳性细胞因所在 cluster 未通过主规则而在
`malignancy_call_v7` 中保持 `uncertain`。相应地，V7 的“严格交集恶性集合”和
“cluster 主恶性集合”应作为两个并行的预定义分析口径，不能假定前者完全包含在
后者中。

实际 tier 计数：

| tier | 细胞数 |
|---|---:|
| `high_confidence_malignant` | 238,418 |
| `cluster_supported_malignant` | 306,693 |
| `sensitive_union_only` | 74,390 |
| `uncertain` | 380,870 |
| `cnv_no_malignancy_detected` | 110,465 |
| `not_evaluable` | 45,284 |
| `not_candidate` | 3,520,667 |

## 10. 字段选择与下游用法

| 目的 | 首选字段 | 推荐条件 |
|---|---|---|
| 大类展示 | `major_celltype_v7` | 可显示全部；严谨比较另加 strict 过滤 |
| 严格大类分析 | `strict_major_eligible_v7` | `True`，并排除 `exclude_from_primary_major_v2` |
| 亚型展示 | `final_annotation_v7` | `_unresolved_GENE` 仅作描述 |
| 严格亚型分析 | `strict_downstream_subtype_eligible_v7` | `True` |
| cluster 主恶性分析 | `is_malignant_primary_v7` | `True` 为恶性；`False` 仅代表 CNV 未检出；NA 不填补 |
| 严格交集恶性集合 | `strict_malignant_downstream_eligible_v7` | `True`，即两方法交集阳性；其中 30,080 个不属于 cluster 主恶性集合 |
| 敏感恶性探索 | `is_malignant_sensitive_v7` | `True`，即两方法并集阳性 |
| 方法一致性 | `malignant_consensus_v7` | 区分交集阳性、一致阴性和方法冲突 |
| CNV 可评估性 | `infercnv_run_status_v7`、`copykat_run_status_v7`、`cnv_skip_reason_v7` | 必须保留 skipped/not evaluable |
| 作图显示 | `final_annotation_with_malignancy_v7` | 仅作标签展示，不代替上述布尔/证据字段 |

`is_malignant_primary_v7` 是 nullable Boolean：只有 cluster malignant 为 True、
`cnv_no_malignancy_detected` 为 False；uncertain、not_evaluable 和 not_candidate 都是
NA。禁止用 `fillna(False)` 把它们变成非恶性。

条件比较、差异表达、通讯和网络统计必须以原始样本/患者为生物学重复，优先采用
sample-level pseudobulk、分层模型或逐样本汇总；不能把 pooled cells 当成独立重复。

## 11. 映射与质量闸门

V7 按 `cell_id + sample_uid + global_position` 进行逐样本结果回填，不使用裸 barcode
跨样本匹配。最终审计结果：

- 4,676,787 个细胞全部且仅映射一次；
- `mapping_conflicts = 0`；
- `unmapped_cells = 0`；
- 10,432/10,432 个 counts 基因可映射到固定坐标；
- raw counts 抽检为非负整数；
- inferCNVpy/Copykat_python 最终 `failed = 0`；
- 两个冻结对象的 SHA-256 与 `FROZEN_V7.json` 一致。

689 个样本至少有一种方法被跳过或部分评估。这是方法覆盖信息，不是映射失败。
逐样本解释见 `tables/skipped_samples_with_explanations.tsv`。

## 12. 解释边界

1. inferCNVpy 和 Copykat_python 提供转录组 CNV 证据，不是病理或 DNA CNV 金标准。
2. `cnv_no_malignancy_detected` 表示当前方法未检出，并不证明细胞生物学良性；
   CNV-low 肿瘤、低深度细胞或参考选择限制仍可能造成假阴性。
3. `not_candidate` 仅表示没有进入癌种特异候选谱系，不能理解为正常。
4. `uncertain`、hierarchical boundary 和 `_unresolved_GENE` 应保留；它们分别代表
   CNV 证据不足、连续谱系边界和描述性 cluster marker，而不是同一种失败。
5. UMAP 只用于局部结构审阅，不能单独确定谱系、亚型或恶性。
6. V7 没有自动重跑旧 V5 分支的 CoVarNet/LIANA。继续使用旧结果时必须标明其
   注释/恶性父版本；若要声称为 V7 下游结果，应以 V7 字段重新构建输入并重跑。

## 13. 主要审计文件

- `config/run_config.json`：实际输入、阈值、癌种 candidate/reference 策略和引擎版本。
- `reports/00_local_preflight.json`：本地依赖、counts、基因坐标和固定提交预检。
- `attempts/manifest_v7_relaxed_candidate/manifests/manifest_audit.json`：角色重建审计。
- `tables/per_sample_cnv_targets.tsv`：逐样本细胞数、角色数和方法 eligibility。
- `manifests/infercnv_run_summary.tsv`、`copykat_run_summary.tsv`：逐样本运行状态。
- `tables/malignant_cluster_summary.tsv`：cluster 级恶性规则输入与结论。
- `tables/integrated_cnv_mapping_summary.tsv`：回填映射审计。
- `tables/skipped_samples_with_explanations.tsv`：所有跳过/部分评估原因。
- `tables/major_v5_to_v7_transition.tsv`：V5 到 V7 大类差异。
- `tables/malignancy_v3_to_v7_transition.tsv`：旧恶性 v3 到 V7 的差异。
- `reports/03_v7_integration_summary.json`：整合质量闸门。
- `reports/04_v7_freeze_summary.json` 和 `FROZEN_V7.json`：最终冻结统计。

## 14. 参考的历史文档

本说明整合并更新了以下先前文档的核心口径：

1. `07_fine_annotation/20260823_corrected_evidence_v1/ANNOTATION_SYSTEM_SUMMARY_CN.md`
2. `07_fine_annotation/20260824_refined_subtypes_v2/REFINED_LOGIC_CN.md`
3. `07_fine_annotation/20260905_uncertainty_resolution_v1/reports/UNCERTAINTY_RESOLUTION_V2_CN.md`
4. `08_downstream_analysis/20260905_uncertainty_annotation_visualization_v1/reports/ANNOTATION_STATUS_INTERPRETATION_V1_CN.md`
5. `10_final_v5_malignancy_covarnet_liana/20260825_v5_final_v1/reports/ANNOTATION_MALIGNANCY_SCHEMA_CN.md`
6. `10_final_v5_malignancy_covarnet_liana/20260825_v5_final_v1/reports/FINAL_V5_MALIGNANCY_COVARNET_LIANA_REPORT_CN.md`

若历史文档与本说明在“当前默认版本”或恶性证据数字上不一致，以 V7 的
`FROZEN_V7.json`、实际字段和本说明为准；历史文档继续用于解释方法演进。

