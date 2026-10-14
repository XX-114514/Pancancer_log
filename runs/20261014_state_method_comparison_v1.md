# cNMF完整K稳定性与Spectra全因子对照 v1

2026-10-14T01:55:32.953664+00:00 /UTC+8 2026-10-14T09:55:32.953664+08:00冻结，未开始分析指标。 [计划](../inventories/analysis/2026-10-14/state_method_comparison_plan_v1.json)。等待原cNMF真实terminal、全600+6consensus/SHA及worker absent后，官方cnmf1.7隔离环境入口`python compare_state_cnmf_spectra_v1.py`，PANCANCER_PREFLIGHT_ROOT/PANCANCER_STATE_COMPARISON_PLAN；2线程，新目录不覆盖。

全部90cNMF候选/9000repeat配对/1300crossK/14355Spectra对照。官方Kstats明确无密度过滤，不与density0.5共识混写；两种gene-spectrum和absolute/normalized usage并列，原输入/outputs不得改。冻结探索repeat/crossK门槛只是程序候选稳定性，不是状态生物学/跨患者验证；不按clinical/LR p选K。所有K保留，其他谱系/原9CM/UCell/state-DIALOGUE/空间保持完整目标。

2026-10-14T02:02:46.571061+00:00 /UTC+8 2026-10-14T10:02:46.571061+08:00代码审阅修正、尚未执行v1/v2：恒定score/零gene-vector数学不可估计不能显示0相关。新[v2计划](../inventories/analysis/2026-10-14/state_method_comparison_plan_v2.json)/v2新输出明确NaN及defined字段，全部范围/门槛保持，v1保留未运行；入口改`compare_state_cnmf_spectra_v2.py`。
