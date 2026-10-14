# 全18大类cNMF全部K稳定性检查

原18major K10/15/20×100仍训练；仅完成major的producer completed记录和完整300iteration/3consensus可进入诊断。

[事前合同](../inventories/analysis/2026-10-14/all_lineage_cnmf_diagnostic_plan_v1.json)。逐45programme/4500topic匹配/650crossK关系，全部K和重复保留。旧repeat/crossK门槛保持；零profile余弦NA/恒定分数rankNA，惩罚成本完成匹配清单但不把undefined写0，>=90%finite重复匹配才具备repeat候选资格。主top50只取正系数，不用0权重基因补足；实际签名size列出。

源/node及RNA复杂度、原生同cell支持一起诊断，source-high不默认artifact；跨K重复不当独立states。统计调用前后原训练输出全部SHA不变。入口`diagnose_all_lineage_cnmf_v1.py --major <frozen_original_major> --plan <frozen_plan>`，PANCANCER_PREFLIGHT_ROOT。尚未运行实际诊断，不做cell级p、临床/LR或正式状态。
