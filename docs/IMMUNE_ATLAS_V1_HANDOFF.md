# 免疫治疗单细胞档案 v1：下游使用与验收

## 成品与发布门槛

运行根为`${PROJECT_ROOT}/Project_v5/runs/20261013_immunotherapy_atlas_v1`，正式对象为其`outputs/immune_validation_authority_v1.zarr`。临时卷上的candidate和历史CNV目录均不是正式入口。

固定范围9来源、400原始文库、2,670,057细胞、49,563原生特征。正式使用前，必须同时存在并通过：

1. `manifests/immune_authority_release_public.json`：全行计数/污染扣除/X验收及复制后的所有文件校验。
2. `manifests/immune_final_interface_audit_v1_public.json`：270节点、全量图、原库CNV映射和成品接口。
3. `manifests/CM_final_archive_recount_v1_public.json`：从正式对象重新计算15范围的原始CM输入频率。

只有实际成品通过才可使用；本协议本身不证明发布完成。公开汇总与当前阶段见[运行记录](../runs/20261013_immunotherapy_atlas_v1.md)。

## 数据与主要字段

| 槽位/字段 | 含义与用途 |
|---|---|
| `layers/counts` | 全部原生整数计数；原始RNA及CNV复核入口 |
| `layers/ambient_removed_counts` | 实际采用的污染扣除量；没有采用污染估计的细胞不伪称去污完成 |
| `X` | 选用计数的log1p CP10K；分母为本细胞完整原生特征总量 |
| `varm/measured_in_source`、`uns/source_order` | 哪些特征在各来源可测；未测与表达为零分别处理 |
| `major`、`majorCluster` | 多基因证据的大类身份与CoVarNet兼容字段 |
| `subcluster`、`subCluster` | `major_marker_number`展示名；一个marker是群描述，不能单独证明身份 |
| `node_id` | 独立稳定节点键；数字不是旧Project_v3群的对应序号 |
| 节点catalog的`subtype_candidate`及证据字段 | 保存在`uns/subcluster_catalog_json`的细身份候选；尚未升级为权威成熟细胞亚型 |
| `state_score_*` | naive/memory、cytotoxic、dysfunction、IFN、cycling独立连续读出 |
| `obsm/X_pca`、`X_pca_harmony`、`X_umap` | 所有细胞的坐标；UMAP有明确的拟合子集和全量转换范围 |
| `obsp` | 全量邻接图；整合不覆盖原始RNA计数 |
| `dataset`、`sampleID`、`capture_scope`、组织/时间字段 | 原研究、物理文库和捕获材料；不能仅凭barcode关联跨研究患者 |
| `patient_block`、原生临床终点及资格字段 | 只按已经核验的来源映射；未知保持缺失 |
| `malignant_status`、`cnv_evaluable`、`CNV_method_status` | 恶性候选结果、可评估范围与实际方法状态 |
| `CNV_native_total_counts` | 同原库CNV归一化所用完整原生RNA分母，数值型；未评估可缺失 |
| `uns/CM_scope_reports_json`、`CM_node_membership_json` | 新CM发现范围、实际网络支持与节点归属 |

## 分析选择

- 原生RNA分析使用计数层并保留来源可测mask；Harmony是坐标去批次，不是校正后的计数矩阵。
- 细身份候选和状态分开报告。270群中264使用本群原生表达对比、5仅有非特异代表marker、1没有有效marker。XIST、DNAJ等描述不自动命名新细胞类型。
- `CNV_supported_malignant`是同原库单infercnvpy支持。阴性候选、谱系参考、未评估均不是双方法正常金标准。CAR产品、免疫富集和离体刺激不用于完整肿瘤阴性推断。旧raw整数直接输入的CNV结果已隔离，不得混用。
- CM分别使用完整TME、免疫富集、CAR产品、T谱系及材料分层。15个可计算范围有10个支持网络成员、5个无支持成员，另6个准备范围不足；范围重叠，不能相加成独立模块数。
- CM保存的是本轮新发现。旧Project_v3固定9CM的外部复现需冻结多基因/节点映射、原W与捕获覆盖，不可按相同数字或单marker对接。
- 临床验证先限定来源、治疗组、基线/随访和终点。PR/SD/PD、病理退缩、CAR产品3月R/NR及疫苗时点不能合并为同一个疗效标签。原生终点资格不等于独立患者数量。
- 免疫富集、血液、正常邻近、离体刺激、核分选和产品保留材料区别。少患者或捕获不完整的队列可做限定表达参照，不能补足完整TME组成证据。

## 权威范围

该版本可作为可追溯的原生计数、来源、分群与候选结果档案。`independent_validation_frozen=false`、来源最终AC状态及单方法CNV等级保留。技术发布不代表所有患者结局、细身份或旧CM已完成独立临床验证。

本轮未新增公共矩阵下载，也未修改原权威Zarr。完整方法、QC复用/污染评估范围见[方法合同](../methods/immune-specific-atlas-v1.md)。
