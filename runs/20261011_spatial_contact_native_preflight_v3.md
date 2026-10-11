# 固定候选的原生空间表达预检

最新取证：2026-10-11 19:27:55.750429 UTC+8 / 11:27:55.750429 UTC。v3正式接受：19:27:38.038378 UTC+8，证据v14，0状态冲突，本里程碑作业终止。原注释v20261010、172节点、9CM、全发现BH和Zarr保持。

## 问题、范围与结论

检查13条固定接触记录（9种LR）在5份CosMx文件/254,458细胞中是否有真实可测端点。全部65个记录—文件组合保留，缺失/零/覆盖不足均报告；另2条程序记录未混入。没有新增通讯p值或正式状态赋值。

关键结果：

| 口径 | 有至少1个FOV双方各≥20阳性细胞的记录—文件组合 | 原候选记录数 |
| --- | ---: | ---: |
| 全细胞，原生计数>0 | 21/65 | 10 |
| 通用QC，计数>0 | 21/65 | 10 |
| 通用QC，计数≥2 | 0/65 | 0 |
| 通用QC，超过负探针背景，RNA谱系 | 5/65 | 5 |
| 上一口径再排除混合RNA谱系 | 0/65 | 0 |
| 通用QC，超过背景，作者粗谱系 | 13/65 | 7 |

“0符合门槛”不等于0表达；≥20/FOV是预声明的描述性覆盖条件，不是通讯成立的定义或p值。作者粗谱系口径19组合不可估计（13个未获取作者标签、6个ICON没有独立B类别）；两组合缺HLA-DRA，缺失未填阴性。余31组合覆盖不足。作者粗谱系口径的13组合去重复后为4种LR，在LUAD TMA1/2和MESO TMA2有覆盖：LGALS9–PTPRC、LGALS9–CD44、HLA-DRA–CD4、TNFSF10–TNFRSF10B。共享程序/发现背景记录不重复计成独立空间证据。

以上是可测性。没有完成物理结合、邻近、通路响应或CM驱动验证。严格RNA背景口径的5记录只集中于MESO1同一个FOV，多数与竞争谱系共表达；正确作者标签改变可测端点位置，不能把原结果简单解释为T细胞机制。作者注释来自本空间处理，非正交独立真值。

## 方法与分母

- 全原生X计算RNA总量/检测基因；精确唯一基因别名，缺失/歧义分别保留。复合物需要全部组件在同一细胞表达；HLA-DRA是MHC-II代理，不能仅凭alpha转录本证明完整结合。
- 按原FOV+cell_ID找回5份CSV中的每份10负探针；每FOV所有细胞×负探针值的95百分位，目标基因严格超过该阈值并>0。这是背景诊断，不是逐基因FDR。
- 通用QC为RNA>30、Area有限正且小于5倍文件几何均值；是敏感性，未声称复现GSE287472作者QC。保留全部254,458缓存细胞，QC口径209,157细胞。
- RNA身份需要至少2个独立谱系marker+CCL5/IL7R/MS4A1描述基因，排除当前LR/命名基因；混合谱系只称RNA共表达代理，不诊断分割错误。不把表达程序或作者Treg/proliferating标签正式转成新状态。
- CCL5/IL7R重叠记录单独计数，不能把同一细胞当接触。固定programme按原model_key映射，覆盖不足不隐藏端点结果。
- 作者粗谱系敏感性用原T/B类别替代双RNA谱系marker，保留同一QC、背景、命名、LR条件；非人工指定新端点。患者/core映射和独立物理标定未完成，FOV/文件/平台不能当独立患者。

[预声明v1](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_plan_v1.json)、[v2类型修复](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_plan_v2.json)、[v3作者修正及敏感性](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_plan_v3.json)。原始大对象/RDS/细胞键仅留服务器。

## 执行、失败与审阅更正

| UTC+8 | 事件 | 结果 |
| --- | --- | --- |
| 19:09:41.873332 | 原计划 | 全65检查预声明；19:10:18.733176普通push |
| 19:12:23.632543—19:13:23.073664 | v1生产 | exit0 |
| 19:14:28.110588—19:14:32.676149 | v1审计 | exit1，细胞键object数组不允许pickle读回；保留失败缓存 |
| 19:15:30.129747 | v2计划 | 显式Unicode；合成键回归通过；19:15:47.045144普通push |
| 19:16:03.222726—19:17:00.981835 | v2生产 | exit0 |
| 19:17:32.509656—19:18:56.559578 | v2审计 | 计算重放PASS；未独立确认作者键语义 |
| 19:19:27.417489 | 证据v13 | 接受计算预检；随后发现作者方向假设错误，作者结果被v3明确替代 |
| 19:22:54.765780 | v3计划 | 修正方向并预声明作者粗谱系敏感性；19:23:46.076113普通push |
| 19:23:52.834498—19:24:55.036235 | v3生产 | exit0 |
| 19:25:24.493557—19:26:55.545403 | v3审计 | 独立原计数、背景、全部作者键及390口径复算PASS |
| 19:27:38.038378 | 证据v14 | 接受原生预检，0冲突；没有空间/机制晋升 |

作者原RDS条码实际为TMA_cell_ID_FOV，我先前误当TMA_FOV_cell_ID；正确方向181,056条全部匹配、0孤儿，错误方向仅5,570偶合。v2作者一致性计数无生物学意义，不以审计PASS掩盖错误。官方[Seurat读取器](https://github.com/satijalab/seurat/blob/master/R/preprocessing.R)与[作者调用代码](https://github.com/KChen-lab/ST_Comparison/blob/main/Seurat/CosMx_ICON.R)支持方向；[全量方向证据](../inventories/analysis/2026-10-11/spatial_author_barcode_orientation_v1.json)保存源SHA/数量。原作者[注释说明](https://github.com/KChen-lab/ST_Comparison/blob/main/Annotations/READ.me)说明QC过滤的细胞未包括，缺标签不当阴性。

GSE287472作者raw路径404/API403限流，未获得作者QC代码；已记录而非暗示已验证。v1/v2/v13均保留。关键节点普通push收据和HTTPS核验另附，不伪造事件时点。

## 复算与产物

scanpy Python3.10.19、NumPy1.26.4、pandas2.3.3、SciPy1.11.4及pyreadr；生产/审计代码SHA在[终态](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/checkpoint.json)。调用：

```bash
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" PANCANCER_AUTHOR_ROOT="${AUTHOR_INPUT_ROOT}" python scripts/spatial_contact_native_preflight_v3.py
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" PANCANCER_AUTHOR_ROOT="${AUTHOR_INPUT_ROOT}" python scripts/audit_spatial_contact_native_preflight_v3.py
```

目标新目录的predeclared_plan.json必须先存在，脚本拒绝重复产出。独立审计未import生产者：HDF5 CSR直接读、整表CSV MultiIndex对齐、排序统计背景阈值、作者标签独立reindex；原生选定基因14,856,530值、负探针2,544,580值、65组合×6口径、30,654 FOV口径核对。新作者全键零孤儿，原文件/分母/坐标不变。计算独立不等于独立生物学验证。

- [全部65组合](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/all_contact_record_testability.tsv)
- [QC/坐标/单位](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/coordinate_donor_and_segmentation_review.tsv)
- [全部选定基因](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/all_selected_gene_cell_QC.tsv)
- [正确作者谱系比较](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/author_lineage_agreement.tsv)
- [独立审计](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/independent_audit.json)、[登记](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/registration_receipt.json)、[UTC/UTC+8时间线](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/timeline.tsv)
- [覆盖图](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/figures/all9_interfaces_native_spatial_coverage.png)、[图源](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/figures/figure_source.tsv)、[像素审阅](../inventories/analysis/2026-10-11/spatial_contact_native_preflight_v3/figures/visual_review.json)

## 下一步

先核实固定9LR的结合/复合物证据、物理坐标及组织/FOV/患者core边界；预先固定去自接触、CCL5/IL7R重叠、RNA深度/局部密度条件化的邻近零模型。全部固定接口报告，重复13记录映射一次，不按当前结果更换端点。先对可测关系作条件空间检验，稀疏关系明确记录功效；癌种限定关系保留价值，跨癌种转移与同癌种复现分开。患者未知时只作条件FOV描述性检验，不声称患者独立验证。目标仍active，尚未完成通讯轴验证。
