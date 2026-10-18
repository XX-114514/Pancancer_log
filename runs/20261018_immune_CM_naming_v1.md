# 2026-10-14免疫图谱CM（v1）全量暂定命名与临床解释

## 任务与已完成范围

用户要求命名本次CM、解释生物学及进一步临床价值。2026-10-18只读重核15范围全部79因子实例、62有网络成员/17无成员，命名不修改底层编号/成员/W/H/Zarr；不与v20261010九CM或后续MCP混称。

[完整解释与临床路线](../methods/immune-CM-naming-and-clinical-interpretation-20261018.md) · [79名称中文目录](../inventories/analysis/2026-10-18/immune_CM_naming_v1/CM_NAMES_CN.md) · [对应表](../inventories/analysis/2026-10-18/immune_CM_naming_v1/CM_name_catalog_v1.tsv) · [输入SHA及验收](../inventories/analysis/2026-10-18/immune_CM_naming_v1/naming_receipt_public.json)。

## 执行与验证

服务器已有scanpy Python/pandas/numpy，只读输入`${PROJECT_ROOT}/Project_v5/runs/20261013_immunotherapy_atlas_v1`的最终catalog、CM索引及各scope最终attempt的W/H/成员/边/元数据。脚本`summarize_immune_CM_naming_evidence_v1.py`从stdin执行，无新科学模型/临床检验。后续加入原始W集中度与实质载荷边的描述诊断，第二次完整提取成功；旧版证据保留。

`build_immune_CM_name_catalog_v1.py`生成79个scope-qualified别名，重点18个人工依据核心拟定，其余按原节点描述，不强行赋予机制。79键唯一、每行原成员精确一致、全部功能/临床晋级标志false、目录SHA一致、源码compile通过。权重归一化/节点覆盖及W-H轴、来源块映射在真实全输入中断言通过。完整数值载荷只留本地工作证据并由.gitignore保护；公共仅名称、节点别名、方法、源SHA与日志，不公开患者/细胞行、H或学习权重。

已读服务器pan-cancer-covarnet-translational-analysis及feature guardrails，按本次已验收scope文件替代skill历史h5ad布局；未为满足旧路径而重建/重拟。文献核对CoVarNet、TCF7与TLS原研究，以及GEO官方的TNBC/NSCLC/HNSCC设计。网页Nature/PMC部分受跳转或验证码限制，使用可读原论文/官方登记及既有精确本地合同，不把不可读网页当新证据。

## 关键决定

原W出现近单节点及载荷/网络不一致诊断，不能按成员表行数就赋复杂机制。该诊断只针对原W，nsNMF实际F=W S H和原临床A=W S合同明确保留；1%/90%为描述门槛，非原CM接受门槛。whole_TME来源分布偏移、患者未核单位、XIST及B/Plasma候选身份都限制功能命名。

已核2026-10-16头颈癌四范围30组件基线临床全家族及独立审核完成、无可晋级稳健CM。本轮不重开旧家族或搬用其他表达程序结论。建议下一解决身份/测量，再进行真正独立同癌种病例、预先冻结问题、功能/空间及临床增量；仅有时间/事件才进入生存。新名称并未证明疗效、耐药或直接通讯。

## 时间与交接

取证/目录生成时间见naming_receipt_public.json；用户时间采用Australia/Melbourne，日志同时记录UTC/UTC+8。完整原目录保持；命名作为独立解释sidecar。日志见[本任务](../logs/codex/2026/2026-10-18_immune-CM-naming-and-clinical-interpretation.md)。普通commit/push及远端SHA另核后再声明同步成功。
