# TME 注释面板 v2（合并版）

2026-11-05上皮逐基因审阅补充：[185模块兼容版本](../v7_epithelial_extension_v1/TME_V2_COMPATIBILITY.md)保留Lodi模块，拆分本版本中的胃主/壁、前列腺基底/club和肾脏混合身份。下文及原JSON保留为179模块历史版本；后续上皮审阅应显式使用兼容文件。

由 `scripts/annotation_v2/build_panel_registry_v2.py` 生成，输入是 V7 registry 和 Lodi 作者仓库 @35f054b。注释逻辑见 [methods/tme-annotation-v2-lodi-workflow-20261105.md](../../../methods/tme-annotation-v2-lodi-workflow-20261105.md)。这是证据注册表，不是细胞分类器：`auto_final_label` 均为 false，不选赢家。

| 文件 | 内容 |
|---|---|
| [registry.json](registry.json) | 179个模块，schema 与 V7 相同，另加 `source_layer` 和 `source_locator`；可由 `scripts/v7_panel_repair_v1/evidence_engine.py` 直接读取（validate_registry 通过） |
| [ALL_PANELS.md](ALL_PANELS.md) / [all_panels.tsv](all_panels.tsv) | 全部模块与基因，按来源分层 |
| [references.tsv](references.tsv) | 49条文献：V7 原有21条，v2 新增28条（PubMed/Crossref 已核对；1条 Lodi 引文未能确认，标为 unresolved） |
| [gene_reference_trace_v2_additions.tsv](gene_reference_trace_v2_additions.tsv) | 新增模块的逐基因依据状态。V7 原有记录见 `../v7_panel_repair_v1/gene_reference_trace.tsv` |
| [lodi_subtype_vocabulary.tsv](lodi_subtype_vocabulary.tsv) | Lodi 的 Intermediate / Minor / Minor_2.0 命名词汇 |
| BROWNE_INTERFERON_RESPONSIVE_GENES.tsv、SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN.tsv | Lodi 使用的 MSigDB 基因集，通过 msigdbr 25.1.1 取得（MSigDB 2025.1.Hs），附版本 |

## 来源层

- **V7_panel_repair_v1（112）**：原样保留。
- **Lodi2025_code（20）**：
  - 29基因大类检查，9个模块；
  - 5个癌种特有群体；
  - hypoxia（15基因）和上皮分化（21基因）；
  - 4个外部命名集（IFN、EMT、stress、IG）。stress 列表未公开。
- **Lodi2025_marker_genes（21）**：作者 FindAllMarkers 结果中每个亚型 avg_log2FC 前15、pct.1 ≥0.25 的基因，只作参照谱。
- **TME_de_novo_v1（26）**：V7 中没有的器官上皮、骨/软骨、破骨样细胞、脂肪、骨骼肌、apCAF、PI16 成纤维、CXCL9 巨噬。
