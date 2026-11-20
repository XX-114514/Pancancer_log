# R4/K15 下游验证组织与证据登记

本次核查范围为2026-11-15 R4/K15的全部15模块。不是新一轮NMF/CPDB或空间拟合，不将其他版本同号CM混用。

- [下一阶段工作计划](../../../../methods/R4-CM-downstream-workplan-20261120.md)：六个工作包、依赖和验收条件。
- [脱敏技术回执](technical_receipt.json)：15模块、24空间模型、336项哈希检查及独立复核。
- [本轮任务轨迹](../../../../logs/codex/2026/2026-11-20_R4-downstream-organization.md)。
- [执行脚本](../../../../scripts/build_R4_downstream_register_20261120.py)。

服务器逻辑位置为`${RESULTS_ROOT}/R4_downstream_20261120_v1/evidence_register_v1/`：

|产物|用途|
|---|---|
|module_evidence.tsv|15CM网络成员、CPDB事件及合格单元级表达比较|
|spatial_edge_setting_coverage.tsv|每CM/后验丰度摘要/N下正负零缺失设置分母|
|spatial_N_sensitivity_existing_summary.tsv|已有全24参数比较汇总的可追溯链接副本|
|spatial_qualification_register.json|成像身份、测量、患者与图像核查的未解决门槛|
|spatial_qualification_register_v2.json|并行进展补充：作者图像完整性已通过，配准/病理/患者门槛仍开放；保留v1截点|
|receipt.json / independent_audit.json|输入/输出SHA及独立聚合复核范围|

详细科学表、输入、原生数据和个体记录继续留服务器。本目录只有技术状态、方法与哈希。事件数/设置数不能作为独立发现数或用于按数量挑模块。

目前优先继续CM12血管核心的组织内空间验证；CM01/CM04作为组织结构对照；CM13/LIHC另作外部预后验证资格支线。完整三方hub、差异通讯和稳定治疗响应指标仍未建立。

作者HE/spot图像已下载和完整性验收，下一动作是定量网格—像素配准与病理区室核对。不能把完整性通过写成配准通过；直接marker稀疏检出也不能由解卷积结果补作实测证据。
