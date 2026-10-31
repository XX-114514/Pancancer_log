# 癌种名称统一 v1（2026-10-31）

## 问题

当前TME合格范围（770个工作单位、43个来源）中的`cancer_type`沿用各来源的原始写法，共35个字符串，同一癌种存在多种写法（`CRC` / `Colorectal cancer` / `Left/right CRC`；`Lung cancer` / `Non-small cell lung cancer` / `LUAD` 等），并有错误标签（GSE182109标为`Glioblastoma`，但含grade II/III胶质瘤）。按癌种分层或比较前必须先统一。

## 方法

- 映射粒度：**来源级**。770单位表中每个来源只有一个原标签（已核：0个来源内冲突）。
- 本体：OncoTree `oncotree_latest_stable`（2025-10-03发布），从官方API下载（`tumorTypes.json`，SHA256 `132b7b67…60ca`，存服务器`${PC_ROOT}/migration/public_refetch/oncotree_20261031/`）。表中每个代码都在该版本中存在，名称、mainType、tissue直接取自API。
- 证据：每个GSE的GEO系列标题/摘要，以及在需要时的GSM特征字段（`histology`、`tumor subtype`、`tumor grade`、`cancer type`、`tumor type`），于2026-10-31下载，存服务器`geo_series_text_20261031/`。不凭记忆判断组织学。
- 两层标签：
  - `analysis_main`：默认分析分层（23类），用于"按癌种"分析；
  - `analysis_fine` + `oncotree_code`：在证据允许的最细层级，混合来源取共同父节点。
- `resolution`：`exact`（标题或特征明确）、`assumed_dominant_histology`（只写到器官癌，按该器官的主要组织学赋码，例如胃癌→STAD）、`parent_of_mixed`（来源内有多个组织学亚型）、`parent_unresolved`（组织学没有公开字段）、`study_level`（依据研究或合集描述）。

表格：[cancer_type_map_v1.tsv](../inventories/cancer_type/cancer_type_map_v1.tsv)。

## 结果（已观察）

- 770/770单位均有映射，35个原写法归为23个`analysis_main`。有≥2个独立来源的分层只有9个（NSCLC 8、结直肠5、肾细胞癌3、胃3，乳腺/前列腺/卵巢/宫颈/葡萄膜黑色素瘤各2），覆盖595个单位；其余14个癌种都只有单一来源，**按癌种的结论无法与来源效应分开**。
- 按`resolution`统计单位数：exact 272，assumed 253，parent_of_mixed 130，parent_unresolved 106，study_level 5。
- 纠正：GSE182109 `Glioblastoma` → `Diffuse glioma`（DIFG；GSM中grade IV 35、II 8、III 1）。GSE297041的`CESC`没有依据（GEO只写primary cervical cancer）→ CERVIX父节点。

## 使用规则

1. 按癌种分层或比较时使用`analysis_main`；只有证据为`exact`的来源可以直接用`analysis_fine`。
2. `parent_of_mixed`的6个来源（GSE127465、GSE159115、GSE161529、GSE182109、GSE197461、GSE273127）需要把GSM级亚型映射到单位，之后才能做细亚型分析（后续任务）。
3. 本表不修改冻结存储和原`cancer_type`字段；是侧表，按`source`关联。
4. 单来源癌种不能单独支撑"癌种特异"结论。

## 未解决

- GSE148071、GSE179994（NSCLC）、GSE297041（宫颈）的组织学需要从原文补充材料获取。
- CXG SPECTRUM是否全部为HGSOC，需要核合集元数据。
- 单位级亚型映射（上面第2条）。
