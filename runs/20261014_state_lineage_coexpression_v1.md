# 原生同细胞谱系RNA诊断 v1

2026-10-14T02:12:29.693688+00:00 /UTC+8 2026-10-14T10:12:29.693688+08:00冻结，尚未执行。 [计划](../inventories/analysis/2026-10-14/state_lineage_coexpression_plan_v1.json)。240533/5sources原H5 SHA固定，11panel/4foreign同细胞及全90module tail。>=3标记/.0/.01/.05完整RNA比例、owncore>=2是描述门槛，不是正式celltype、doublet/ambient/phagocytosis判定；无p或自动排除。`python check_state_lineage_coexpression_v1.py`，PANCANCER_PREFLIGHT_ROOT/PANCANCER_LINEAGE_DIAGNOSTIC_PLAN。逐cell统计私有、source-major/原node聚合单独review。

## v1真实失败及v2恢复

2026-10-14T02:16:32.712041+00:00 /UTC+8 2026-10-14T10:16:32.712041+08:00：v1实际02:14:32—02:14:34 exit1，worker288108/child288127已退出；groupby对categorical默认observed=False产生未观察空组，assert len(ii)>=1触发。保留日志、execution及首source私有panel量，新[v2](../inventories/analysis/2026-10-14/state_lineage_coexpression_plan_v2.json)仅observed=True，不改cell/genes/门槛，输出新v2。将v1/v2首source全部quantity直接比对。还没有完整诊断结果，不能将失败当生物问题。

## v2真实完成与独立验收

2026-10-14T02:43:07.957687+00:00 /UTC+8 2026-10-14T10:43:07.957687+08:00补录：02:17:12—02:17:42 exit0，02:20:03全7937589值/110/212/990记录复算PASS，v1首source量完全一致；全部55panel完整。强foreignRNA罕见，未判doublet/ambient/phagocytosis，不能作为大范围身份错误证据。SCP/大载荷观察超时，短SSH正常、alias不使用ProxyCommand；IPv4/QoS未证实修复，9014byte tar恢复缺的2file，全SHA通过，没有重跑science。 [结果](../inventories/analysis/2026-10-14/state_lineage_coexpression_v2/RESULTS.md)。

02:40:10.955108 UTC /10:40:10 UTC+8：numeric/cNMF/pyUCell/native RNA合并handoff登记完成，control135目录/0冲突；本任务已无active scientific job，原CM证据不动，formal state/UCell biological/all9 mapping均false。最新goal继续，下一字典/投影及缺谱系输入。 [登记](../inventories/analysis/2026-10-14/state_lineage_coexpression_v2/registration_receipt.json)。
