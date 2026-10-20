# 真实文库参考背景与来源迁移诊断

承接2026-10-14免疫图谱CM（v1，T_whole_TME_primary_tumor），[效果前方法](../methods/CM-natural-library-reference-diagnostic-20261021.md)、[实际启动技术证据](../inventories/analysis/2026-10-21/CM_natural_library_reference_v1/technical_start.json)。

实际v1启动2026-10-20T22:03:38.199665+00:00 UTC，在原RNA首段后发生小整数聚合索引溢出，尚无模型效果；原代码/失败输出保留。v2仅将块编码转int64，增加90块×19类索引合成例子，配置SHA911954f5bd044ce160a69fe308ac1b209dd9dd2804e0a8a0c198f39b1e85286c保持；实际启动2026-10-20T22:07:28.410790+00:00 UTC，源码SHAff340330f9f27a81ed25d4cda0a33e1ad5b80e1baef31fad1cc1933cacf1131d。scanpy Python3.10/单BLAS、原项目只读，原生RNA/学习profiles/基因尺度/文库真值/预测全RAM。

```text
python -B scripts/run_CM_natural_library_reference_diagnostic_v2.py < PRIVATE_FROZEN_PLAN
```

环境合同：PANCANCER_PREFLIGHT_ROOT、PANCANCER_PRIOR_REFERENCE_CONTROL、PANCANCER_PRIOR_REFERENCE_OUTPUT、PANCANCER_MIXTURE_SOLVER与PANCANCER_RAM_OUTPUT。新任务tmpfs仅小型来源汇总/技术及私有配置，未新增原始数据下载。

完整24类进入文库真值分母，支持的19类进入参考，未覆盖成分单报。64问题/1376目标位置效果前固定；不可估计、零及未覆盖不当阴性，不按结果改变门槛。完整原作业终态、独立来源指标复算及实际科学/图审待单独接受，启动和合成成功不等于科学通过。旧RNA临床家族、正式身份与两版CM保持。


## 2026-10-20T22:14:29.324059+00:00 UTC /2026-10-21T09:14:29.324059+11:00 Australia/Melbourne /2026-10-21T06:14:29.324059+08:00 UTC+8 全部自然文库问题、来源科学与实际四图接受

本固定诊断科学生产于2026-10-20T22:10:38.250141+00:00 UTC实际完成。独立来源聚合审核于2026-10-20T22:11:30.040476+00:00 UTC结束：完整64位置/1376目标、全部来源分母/缺失原因、误差/相关/校准的充分汇总量逐项重放，以及相同文库架构差值代数核对通过。原完整目标KKT与预定直接全基因求解均接受，原生T缓存计数/深度/检出、24类RNA质量和全基因聚合和/精确选择SHA复原；没有声称第二次原生RNA/逐个体独立refit。完整96来源类别元数据及原68非T计数/总分区另实现对齐。

本机仅接收十份1903918字节来源聚合与技术，显式排除私有计划/基因定义/学习数组/个体行，逐SHA/schema核对；完整两尺度/两系数口径/全部目标及源内失覆盖问题私有科学审阅结束。来源科学数量、效应和结论细项保持私有；不由单一改善/相关替代细状态或细胞量资格。旧临床RNA程序证据、正式身份和两版CM保持。

图v1因旧colorbar接口位置参数失败，v2修接口后完整四图实际查看，v3仅去除拥挤低刻度；原版本保留，最终四图全部实际查看接受，科学表不改。当前本有限诊断已收束；下一先明确连续RNA程序/来源与测量对象，再另定确有资格的新临床/空间问题，不反复放宽原阈值。整体细胞互作目标尚未满足，工具goal仍paused而用户已明确人工继续；未虚报active或complete。新随访CSV和公开push特定许可仍待，无重试/绕过。
