# 原CM全数据库通讯提取：输入验收完成

观察UTC 2026-10-12T05:20:43.440767+00:00 / UTC+8 2026-10-12T13:20:43.440767+08:00。producer v2 04:52:14.680248—04:57:34.877336 UTC、audit04:57:34.882807—04:59:08.314816 UTC均exit0；driver379438/children379442、17914已终止，注册活跃空/状态冲突0。[独立审计](independent_audit.json)与[执行/公开投影](public_checkpoint.json)。

## 全范围结果

- 332个不重复文库、8GSE、10个当前审阅情境、34条原边×情境记录，展开68方向；原9CM/172节点不变。332不能称332独立患者。
- 4624个固定consensus LR ×68=314432完整候选覆盖行；未输出不是无通讯。全量扫描6996788原LIANA行、提取381170条目标库×LR观测，不做双rank/p筛选。
- 独立原CSV重放全6996788行、全部观测与40321280数组值（含缺失），全部314432候选/native缺亚基/输出/患者块计数，2180次NNLS KKT通过；这是输入/数值验收，不是生物验证。
- 5377条source/context/direction/LR案例通过20文库、10分析块、≥80%输出覆盖，全部expr_prod有变异。覆盖不保证完整模型可识别，也不意味着显著。
- 去重复后是3118个方向性细胞对×LR接口；598个至少两GSE有覆盖、595个至少两当前癌种标签有覆盖。

[全314432计数/缺失](all314432_public_case_counts.tsv.gz)、[全3118接口覆盖](all3118_directional_LR_interface_coverage.tsv)、[68方向情境](all68_directed_group_summary.tsv)。公开版本删去4个表达/score均值，所有候选键、counts/缺失/flags仍完整；原分数、库/患者keys与数组保留server，SHA可追溯。

## 来源分布与限制

| 来源 | 达数据覆盖门槛案例数 |
|---|---:|
| GSE178341 |4131|
| GSE161529 |426|
| GSE182109 |207|
| GSE183904 |192|
| GSE148071 |171|
| GSE179994 |111|
| GSE270680 |84|
| GSE131907 |55|

4131/5377约77%来自同一CRC来源，案例还包含相同源的多个技术/临床情境，并非4131独立发现。多CM映射共享同一底层LR向量，也不增加独立证据。598多来源仅表示覆盖机会，尚未检验效应一致/留一来源。癌种字段保留旧层级，肺癌/NSCLC表述差异不当不同癌种证据。既有发现来源不是完全独立验证。

旧LIANA固定expr_prop=.1、return_all_lrs=False；官方[参数说明](https://liana.readthedocs.io/en/stable/generated/liana.mt.rank_aggregate.__call__.html)明确过滤条件。不能补零/当阴性；observed numerical zero可保留。native缺亚基与measured-but-no-output分别统计。所有数据库入口尚未经逐项结合/定位审阅，ACTR2警告与C1QA–CD93争议保持；全基因可测不等于同细胞复合物、蛋白结合或接收通路响应。

## 先前分子轴已有可测入口（不是选后显著性）

[预先讨论轴的本次明确端点](preexisting_axes_current_endpoint_coverage.tsv)：这是全库提取后按既有问题举例，不改变全5377案例范围或宣称复现原不同端点。

- COL1A1–DDR1，Fibroblast_SFRP2_17→Epithelial_cell_KRT18_52：乳腺GSE161529观测34/36文库（分析块34，患者未审阅）、CRC GSE178341 21/21（20来源内患者块）、胃癌GSE183904 20/21（20患者块）；三个来源均通过数据覆盖。
- CCL5–CCR1，CD8_T_cell_NKG7_9→Macrophage_IFI6_18：CRC两情境36/36和42/42（29、38分析块）、GBM29/29（13患者块），不能把CRC两个情境当两个来源。
- GRN–TNFRSF1B和SPP1–CD44，Macrophage_IFI6_18→Regulatory_T_cell_CTLA4_15：当前达标记录主要在CRC两个情境，属于来源限定入口，同样保留。

当前只是RNA输出/覆盖，不推断关联方向、跨癌种效应、CM驱动或免疫抑制。原编号、身份置信度和没有正式状态的限制保持；NKG7_9与NKG7_20不能凭共享marker合并。

## 下一项具体工作

对全5377案例及其原CM映射，检查原6列主设计（去两端活动、双方major内比例、双方log细胞数、截距）、实际可识别/共线/影响点/患者重复，再冻结完整检验家族并拟合。结构常数与CM曝光本身无独立变异分开报告；不能把数学不可估计当不存在通讯。主/敏感性、整体BH、患者块方法和来源汇总必须预声明，原H对照保留。源内关联/方向/不确定性完成后，再从留出来源与接收响应验证可行轴。

## 失败历史

v1 04:43:01—04:43:09 UTC exit1，旧audit摘要不存input_signature字符串而保留independently_reconstructed flag；读分数前失败，无科学score。旧目录/脚本/151byte表头/log保留。v2仅按真实schema检查已审flag、固定全局status与逐库JSON签名/SHA及父manifest/contract；范围、指标、门槛、NNLS不变。没有重跑LIANA、下载、删除、身份/canonical修改。
