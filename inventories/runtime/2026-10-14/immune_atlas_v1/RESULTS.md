# Project_v5免疫专用单细胞档案 v1 正式交付

完成时间：2026-10-14 05:08:43 UTC /13:08:43 UTC+8。

## 正式对象

`${PROJECT_ROOT}/Project_v5/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr`

9来源、400原始物理文库、2,670,057细胞、49,563原生特征；实际分配29.816GB，24,431文件复制前后SHA逐一一致。成品独立包含数组，临时卷不作为权威入口。

已完成PCA/Harmony、全量邻接图、全细胞UMAP坐标、18个major类别（含未决/边界）、270描述分群、T谱系16群与5个连续状态轴。264群有本群原生表达对比marker、5群仅非特异代表marker、1群没有有效marker。细身份候选保留证据等级。

## 来源与验证资格

| 来源 | 细胞 | 原库 | 来源内终点资格细胞 | CNV支持恶性 |
|---|---:|---:|---:|---:|
| CAR产品（CXG来源） | 396,593 | 59 | 224,628 | 0 |
| GSE169246 | 489,360 | 78 | 484,326 | 0 |
| GSE176021 | 868,494 | 110 | 0 | 0 |
| GSE256326 | 22,452 | 9 | 0 | 0 |
| GSE286827 | 375,019 | 57 | 0 | 11,867 |
| GSE292046 | 52,525 | 8 | 0 | 14,933 |
| GSE294109 | 130,243 | 25 | 0 | 1,984 |
| GSE305078 | 252,822 | 32 | 0 | 23,873 |
| GSE316195 | 82,549 | 22 | 0 | 0 |

终点资格总计708,954细胞，按来源各自使用；细胞数不是患者数。0个CNV支持恶性不能解释为该队列无肿瘤；捕获不适合、参考不足及单方法阴性候选均有独立状态。临床角色及研究别名最终AC仍需独立版本冻结。

## 恶性预测

同原库完整原生RNA分母的log1p CP10K输入，149原库success、251skip，无方法失败。52,657个CNV支持恶性、57,958不确定、151,735阴性候选、642,560谱系参考、1,760,348 skip和4,799其他未评估。单infercnvpy支持，gold=false；旧raw整数直接输入推断已隔离。

## 新CM发现

15个可计算范围全部完成，10有支持网络成员、5无支持成员；另6个准备范围不足。频率按身份/捕获定义，CNV独立保存，当前发现没有按CNV状态拆分节点。完成CNV归档后从正式对象重新计算全部频率，最大误差约1.11e-16。范围重叠，rank不能相加当独立模块总数。

| 范围 | NMF rank | 支持成员的CM数 | 网络边 |
|---|---:|---:|---:|
| whole_TME | 2 | 2 | 70 |
| immune_enriched | 9 | 9 | 51 |
| CAR_product | 12 | 12 | 36 |
| T_lineage_biopsy | 12 | 12 | 10 |
| whole_TME_primary_tumor | 12 | 12 | 95 |
| T_whole_TME_primary_tumor | 4 | 4 | 7 |
| whole_TME_metastatic_tumor | 2 | 0 | 0 |
| immune_enriched_primary_tumor | 3 | 3 | 13 |
| T_immune_enriched_primary_tumor | 3 | 3 | 2 |
| immune_enriched_blood | 5 | 0 | 0 |
| T_immune_enriched_blood | 5 | 0 | 0 |
| immune_enriched_normal_adjacent | 3 | 3 | 15 |
| T_immune_enriched_normal_adjacent | 3 | 2 | 3 |
| immune_enriched_unknown_material | 2 | 0 | 1 |
| T_immune_enriched_unknown_material | 2 | 0 | 0 |

新发现的模块保留rank调查边界与来源留出不可估计记录。旧Project_v3固定9CM的多基因/原W外部复现独立推进，当前未完成。

## 完整验收证据

- [全行计数与污染/X验收](independent_immune_fullgene_audit_public.json)：全部原生count精确、最大X绝对误差0、mask/ID顺序及全部坐标通过。
- [正式容量门槛](immune_final_storage_budget_public.json)：500GB十进制预留+1GB复制余量，实测fits。
- [实际正式发布](immune_authority_release_public.json)：data4读回通过，完成后有效余量约514.3GB。
- [全24,431文件复制指纹](immune_final_copy_fingerprint_public.json)：私有逐文件清单SHA 25f6117eee3ea6b4428e60df106047e5bd9907ad3e58ac49f1d960f95919ec65。
- [最终接口](immune_final_interface_audit_v1_public.json)：270节点、128,599,006连接稀疏条目与82,771,767距离条目合法、原库CNV映射一致、无外部链接。
- [15范围成品CM频率重计](CM_final_archive_recount_v1_public.json)：source/patient-block精确、全部频率误差<=1e-12。
- [成品来源和资格](immune_final_cohort_summary_public.json)：与冻结原来源、患者键资格/原生终点/离体刺激计数一致。
- [流程实际完成](authority_metadata_recovery_v1_pipeline_terminal_public.json)及[审计exit0](final_interface_audit_v1_terminal_public.json)。

## 使用与边界

该版本为计数/来源/候选注释与新CM结果的技术档案。`independent_validation_frozen=false`、`source_final_AC=false`、`malignancy_gold_standard=false`保留。GSE169246/GSE176021/GSE316195沿用历史QC，本轮未重新估计污染；新来源采用的扣除量、残留及拟合排除字段保留。

[全部270群表](subcluster_catalog_aggregate_public.tsv) · [下游使用协议](../../../../docs/IMMUNE_ATLAS_V1_HANDOFF.md) · [方法](../../../../methods/immune-specific-atlas-v1.md) · [运行记录](../../../../runs/20261013_immunotherapy_atlas_v1.md)
