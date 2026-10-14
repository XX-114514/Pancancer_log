# Spectra全173因子稳定性与来源/身份/复杂度诊断 v1

2026-10-14T01:44:33.011327+00:00 /UTC+8 2026-10-14T09:44:33.011327+08:00，尚未运行。冻结[计划](../inventories/analysis/2026-10-14/Spectra_multiseed_diagnostics_plan_v1.json)。输入3seed12000cells/6512features，首seed固定anchor，所有scope内3对Hungarian matching/cycle；权重cosine≥.9、top50Jaccard≥.5、每适用major的去来源score秩相关≥.7为探索稳定门槛。不是FDR、生物验证或状态命名。

所有173候选保留；按source+node后RNA/检测基因残差关联识别复杂度影响，来源高解释度保留队列/癌种解释，亚型高解释度不自动当新状态。不使用任何临床或LR p值、最优seed选择。`PANCANCER_PREFLIGHT_ROOT`与`PANCANCER_SPECTRA_DIAGNOSTIC_PLAN`，scanpy Python、2线程；新目录拒绝覆盖。标量cosine独立复算及正/负/常量秩相关边界检查，完整SHA/axes/finite门槛。

运行入口：`python diagnose_spectra_multiseed_v1.py`。仅汇总表公开；cell/patient/score向量仍私有。cNMF全grid运行，不重新启动；固定UCell/all9CM/state-DIALOGUE/空间仍待。

## 实际结果及独立审计

2026-10-14T01:49:26.197281+00:00 /UTC+8 2026-10-14T09:49:26.197281+08:00取证补录：01:45真实exit0，01:47:41全519/957/638/173规则复算通过；8候选稳定、2复杂度标记；3零prior交集描述在新安全表纠正，原表不改。仅计算候选，无正式state/UCell。 [结果](../inventories/analysis/2026-10-14/Spectra_multiseed_diagnostics_v1/RESULTS.md)。

01:50:11 UTC /09:50 UTC+8：自己的numeric诊断handoff及state evidence register已登记，原CM v17不改；control refresh 135目录/0冲突。cNMF原worker/child命令核实仍运行，406/600factorization文件；T全K已完成、Mac尚待。 [交接](../inventories/analysis/2026-10-14/Spectra_multiseed_diagnostics_v1/registration_receipt.json)。
