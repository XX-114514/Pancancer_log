# 全细胞背景图谱：分层注释与整合比较实施

## 完整目标与当前阶段

按用户明确实施指令执行[原完整方案](../methods/full-context-atlas-hierarchical-annotation-plan-20261023.md)。八类要求见[验收登记](../inventories/analysis/2026-10-23/full_context_annotation_v1/objective_acceptance_register_v1.json)，不因某个读取/合成通过缩小成功定义。当前为stage1：既有全来源原生计数的可读索引、来源内标签继承和结构审计；新多marker身份、A/B/C、全量模型和最终交付仍待。

## 输入、资源与约束

输入为既有冻结全范围合同、主分区集合和旧原生库/精确细胞roster；具体来源行数、原始路径、资源/组配额及模型资格仅私有。现场核实际科学环境、GPU和现有软件，未安装scVI/scArches或更改共享环境。现行空闲预留下大写入未资格，先复用矩阵并限制技术目录大小，不删除或覆盖历史。

父项目AGENTS禁止删除；旧annotation_rebuild规则仅在其范围适用，本次只读其原生输入，新目录位于Project_v5。代码、索引、配置和来源汇总只写新版本目录；不把来源聚合或患者字段返回为个体行。

## 已执行

四份实际脚本SHA绑定见[启动技术收据](../inventories/analysis/2026-10-23/full_context_annotation_v1/technical_stage1_launch_v1.json)。新FullContextAtlas桥接原生旧库与现有分区；核精确roster、样本范围、原生H5AD基因/细胞顺序和source mask，保留原计数全基因分母及未决作者类别。此读取器不提名新生物身份，也不使用临床结果训练。

科学环境六个数值反例actual exit0，涵盖乱序/重复/排除样本、未测与真实零、完整RNA分母、别名冲突、输入SHA/计数改变及原生细胞顺序错位。首次真实执行已在原生读取前冻结，启动UTC 2026-10-23T11:30:33.896217+00:00；观察UTC 2026-10-23T11:33:12.503288+00:00 / UTC+8 2026-10-23T19:33:12.503288+08:00，精确PID/startticks同一进程live，来源检查点推进。

执行入口为已冻结scanpy Python，OMP/BLAS/MKL/NumExpr各4线程，墙时5400秒，目录技术产物上限16MiB。真实软件版本见私有preflight，实际输出/源聚合留私有。进程存活、完整结构审计、独立复算、身份/科学发布和最终交付分别验收。

## 下一步

等待同一handle完整终态，不能因观察超时重启；逐来源计数/单位/可测和原生CSR重放完成后独立验收，再审来源/分选/病例桥与多marker粗身份。类型内RNA状态、稀有/未知和A/B/C的资格/范围/种子须先冻结。旧CM、正式注释、权威计数及并行TDLN任务保持。

本阶段运行中，不称单体Zarr完成或新身份通过。public push此前自动审批拒绝保持，普通本地归档后云端待同步。


## 2026-10-23T11:51:00.455728+00:00 UTC / 2026-10-23T19:51:00.455728+08:00 UTC+8: full structure accepted; native evidence running

Stage1 complete outputs and original process absence were independently checked.
The auditor returned actual exit0; every original legacy source major-count table
matched, the complete code/input/output SHA chain and directory cap passed, and
no expression arrays were written. This accepts structure and inherited labels,
not new biological identity or cross-study independence. See the technical
acceptance in inventories/analysis/2026-10-23/full_context_annotation_v1/.

The first driver did not capture its OS exit status: its terminal exit0 field is
a producer claim, not OS-exit proof. Acceptance uses the complete outputs,
independent full-source count/CSR checks and the auditor's actual exit0. A new
supervisor now captures the real scientific child exit for subsequent stages.

Stage2 was frozen and launched after stage1 acceptance. Three native-evidence
counterexamples passed in the actual scientific environment. The exact supervisor
and child were live and checkpoints advanced. Full native count blocks are read
with 2048 selected rows, four threads, a 10800-second wall budget and a 16MiB
aggregate-output cap. Source-unmeasured genes remain null; zero-native-RNA cells
stay in the archive but not expression denominators. Only inherited-group marker
evidence is produced, with no identity, state, malignancy, clinical P or A/B/C
promotion. The original count matrices, old CM and parallel TDLN task remain.

All frozen source keys were joined to the existing design registry using the
actual legacy columns (including cancerType). An earlier auxiliary schema probe
guessed a different cancer column; it did not use missing output as a scientific
conclusion. Source aliases and patient/phase/endpoint qualification remain open.

An initial local archival heredoc failed UTF-8 parsing and the subsequent checker
used the parent working directory; neither performed the intended update/check.
This explicit UTF-8 file and repository cwd repair the documentation step without
changing science. Full repository validation follows. Goal remains active;
this turn is progress, not a completion or a blocked audit. Public push retains
its prior auto-review restriction; local archive is pending cloud synchronization.
