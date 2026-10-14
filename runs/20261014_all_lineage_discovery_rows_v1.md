# 全18major状态发现均衡行选择 v1

2026-10-14T04:12:18.067948+00:00 /UTC+8 2026-10-14T12:12:18.067948+08:00冻结未运行。[计划](../inventories/analysis/2026-10-14/all_lineage_discovery_rows_plan_v1.json)。全部18原major，各<=6000cells/<=6source，按继承块/库/细胞数metadata排名、source/block均衡fixed seed；single-source不删，标来源限定。共轴>=5000/cells>=1000才标model可训练，未满足保留并另报，不能强填缺基因0或强命名state。所有4532933cells仍保留后续流式评分，当前只选<=108000发现行，无RNA拷贝，500GiB预留和training预算后再提取。

入口`select_all_lineage_state_discovery_rows_v1.py`，PANCANCER_PREFLIGHT_ROOT/PANCANCER_ANNOTATION_ROOT/PANCANCER_STATE_DISCOVERY_ROW_PLAN。全9CM/172原identity不变，不以当前2major完成goal；不按LR/临床p选择。

## 实际终态和独立审计

2026-10-14T04:27:54.083445+00:00 /UTC+8 2026-10-14T12:27:54.083445+08:00取证：真实100766行/18大类，04:15 terminal exit0、04:26独立PASS，未提取新RNA。详见[完整结果](../inventories/analysis/2026-10-14/all_lineage_discovery_rows_v1/RESULTS.md)。
