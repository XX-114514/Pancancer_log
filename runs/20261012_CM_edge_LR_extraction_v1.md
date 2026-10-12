# Run：原CM全数据库LR提取与可比性

预声明UTC 2026-10-12T04:40:29.183334+00:00 / UTC+8 2026-10-12T12:40:29.183334+08:00，PREDECLARED_NOT_STARTED；继承已审60边/34合格pair×context，不按旧显著LR或top3来源选择。LIANA1.7.1 consensus4624，68双向group×4624=314432候选覆盖记录。

[合同、所有代码/输入SHA](../inventories/analysis/2026-10-12/CM_edge_LR_extraction_execution_plan_v1.json)。完整每样本LIANA不是未过滤RNA全空间：旧运行expr_prop=.1/return_all_lrs=False，官方[参数说明](https://liana.readthedocs.io/en/stable/generated/liana.mt.rank_aggregate.__call__.html)说明输出表达门槛；本轮引用当前已安装1.7.1及历史代码核实值，不套最新默认值。所有未输出保留missing，native缺亚基分别统计，不填零或改别名。

主后续指标log1p(expr_prod)，lr_means/两rank保留但本轮不测试p；rank不是绝对强度/FDR。原身份/块、W/H/F保持固定。NNLS去两端活动与H均输出，仅点估计就绪，不当CM驱动。源/情境/候选覆盖可公开，私有库/患者数组只留服务器。

提供PANCANCER_PREFLIGHT_ROOT、PANCANCER_ANNOTATION_ROOT、PANCANCER_CM_LR_PLAN，依次运行producer/auditor；普通push后派发，真实首尾/exit和所有模型接受另记录。

## 2026-10-12T12:47:00.606453+08:00 v1失败与v2修复

v1真实04:43:01—04:43:09 UTC exit1，旧审计摘要不存在input_signature字符串字段，读LR分数前失败，audit未启动；旧目录/代码/log保留，两个PID不存在。新[修复合同](../inventories/analysis/2026-10-12/CM_edge_LR_extraction_execution_plan_v2.json)/v2独立输出，以已审flag、逐库JSON签名/SHA和父合同验收；生物范围/score/门槛不变，v2未运行。

## 2026-10-12T13:20:43.440767+08:00 v2实际完整终态

producer04:52:14—04:57:34、audit04:57:34—04:59:08 UTC均exit0；332库/7M原行/381170观测和全314432cases/40M值及2180KKT独立验收。[结果](../inventories/analysis/2026-10-12/CM_edge_LR_extraction_v2/RESULTS.md)。当前COMPLETE_AUDITED_FULL_LR_INPUT_NO_ASSOCIATION，private保留原分数/块，公开仅counts投影；原v1失败不改，新p尚未计算。
