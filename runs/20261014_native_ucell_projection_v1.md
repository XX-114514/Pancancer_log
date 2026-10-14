# 固定程序字典及原生轴UCell评分 v1

2026-10-14T03:26:23.757554+00:00 /UTC+8 2026-10-14T11:26:23.757554+08:00预声明尚未执行。 [字典计划](../inventories/analysis/2026-10-14/state_programme_dictionary_plan_v1.json)、[评分计划](../inventories/analysis/2026-10-14/native_ucell_projection_plan_v1.json)。全90cNMF/8严格稳定Spectra98sig，main稳定65含跨K重复，不是65state；SP三个matched seed L2权重中位数top50，CN正回归coef top50，TPM替代保留。不按LR/临床p选K/seed/genes。

`build_state_programme_dictionary_v1.py`以PANCANCER_PREFLIGHT_ROOT/PANCANCER_STATE_DICTIONARY_PLAN先生成不覆盖字典；`project_state_ucell_native_v1.py`在官方隔离pyucell0.8环境以PANCANCER_UCELL_PROJECTION_PLAN执行。完整native轴/CPU平均并列秩/主maxrank1500及1000/3000敏感性，缺失coverage至少80%且10genes，之后skip，scope外NaN，实测零保持0。每cell检测>=5genes仅信息支持，不当状态阳性门槛；不按每库百分位强造状态/不改单一原type。源cell IDs/score/sample-node矩阵私有，只公开source-major/node汇总。

500GiB预留+1GiB本轮磁盘预算逐source检查；配额不通过保留partial/execution，不静默继续。原9CM仅0—26.98%core weight当前输入、CM06未覆盖，所有missing谱系/状态LR/接收/DIALOGUE/空间完整目标仍需完成。

## 字典v1展示失败与v2恢复

2026-10-14T03:35:24.953478+00:00 /UTC+8 2026-10-14T11:35:24.953478+08:00：v1完整98记录JSON已写成，但展示TSV join matched_factor_indices的integer列表TypeError，producer exit1；保留JSON和代码，未启动任何评分。新[v2字典](../inventories/analysis/2026-10-14/state_programme_dictionary_plan_v2.json)仅str展示修正，要求与v1全部98record精确相同；[v2评分](../inventories/analysis/2026-10-14/native_ucell_projection_plan_v2.json)仅新字典路径/producer SHA断言/输出版本，所有参数不变。

2026-10-14T03:59:00.501325+00:00 /UTC+8 2026-10-14T11:59:00.501325+08:00真实v2字典03:38:11成功，98原记录精确一致；评分224310/22431103:42:30—03:52:16exit0/退出，5sources全部98/三rank与TPM完成。独立audit派发；formal states仍待。 [字典](../inventories/analysis/2026-10-14/frozen_state_programme_dictionary_v2/RESULTS.md)。

2026-10-14T04:10:17.089296+00:00 /UTC+8 2026-10-14T12:10:17.089296+08:00：03:59:43独立真实UCell audit PASS（非全score再算）；330main中92明显RNA相关、TPM signature rho中位.512。gene轴50/50 complete，保留low gene支持及scope NA；正式state命名/原9CM扩展仍待。[结果](../inventories/analysis/2026-10-14/native_ucell_projection_v2/RESULTS.md)。
