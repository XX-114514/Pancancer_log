# 全18大类有限原生RNA输入提取

[冻结计划](../inventories/analysis/2026-10-14/all_lineage_state_input_plan_v1.json)、[实测容量](../inventories/analysis/2026-10-14/all_lineage_state_input_budget_v1.json)。选择100766细胞，原标签不变；18major按实际选中来源共同轴，一份稀疏int32 counts，保留原生全RNA分母和检测基因数，source缺基因不填0。

预算3GiB/500GiB空闲预留，行指针实测176509672非零上限，rawCSR最大1412883648bytes；取证free507.374GiB。按256选中行/8192core行跨度小块读取，每major复核配额。新矩阵全留服务器，公开只汇总。正式state/全数据投影/9CM/状态后DIALOGUE/空间未完成。

入口extract_all_lineage_state_inputs_v1.py；PANCANCER_PREFLIGHT_ROOT/PANCANCER_ANNOTATION_ROOT/PANCANCER_ALL_LINEAGE_INPUT_PLAN。独立all100766原始count/分母/axis readback后才能拟合。模型自身预算另核。

## v1故障/v2恢复

2026-10-14T04:31:52.815985+00:00 /UTC+8 2026-10-14T12:31:52.817708+08:00：v1 04:30:11—04:30:13 UTC exit1/CSR索引顺序断言，尚无H5AD；首128行sort-only canonical、0重复、RNA总量exact。v2只在读入RAM副本sort_indices，仍断言无重复，权威counts不变。[v2计划](../inventories/analysis/2026-10-14/all_lineage_state_input_plan_v2.json)。所有100766独立COO重组验收后进入模型。
