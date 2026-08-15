# 重要变更记录

本文件只记录会影响项目理解、方法、范围或复现的重要变化。日常细节见 `logs/`。

## 2026-08-16

- 确认 CM Stage01 recovery 在第 581/721 个样本因 `/data4` 的 `USER002` 组 hard quota 耗尽而终止；严格完成合同通过前 580 个样本，Stage 02–08 未开始。
- 用户释放空间后复核仅得到约 2.50 GiB 组配额余量；按 119 个真实 recovery attempts 的磁盘块占用外推，剩余 Stage01 约需 129–157 GiB，因此未提交可预见会再次 EDQUOT 的任务，并把重启门槛设为至少 200 GiB 组配额余量。
- 第 581 个样本的旧 attempt 和原生 CopyKAT prediction 均保留；后续只允许在 cell-ID 一对一覆盖审计通过后于新 append-only attempt 复用。

## 2026-08-12

- 非破坏修补 CM Stage01 completion 合同与 attempt provenance；封存第 462 个旧非终态 attempt，并以 `20260812_124600_cm_stage01_from_0462` 从零基索引 461 重启。新 inferCNV 已成功，CopyKAT heartbeat 活跃。
- 确认免疫治疗四队列在 2026-08-08 后没有新 Stage 13 attempt；19/19、71,398 spots、40-factor merge 仍通过，但最终验收策略冲突和缺失 final audit 未解除。
- 将 CM active v2 Stage 01 更新为 461/721；第 462 个样本 inferCNV 成功但 CopyKAT 仅形成非终态 raw matrix，heartbeat 停于 2026-08-10 08:02，当前无进程，因此状态改为 interrupted/stalled，而非运行中。
- 记录下载 retry v3 于 2026-08-09 07:12 完成 141/141 terminal tasks：127 verified、7 skipped_verified、6 unavailable_upstream、1 failed；GSE201347 partial 继续作为独立 recovery 任务。
- 同步三份 run record、路线图、待办、runs/artifacts inventories 和本日 daily/Codex 日志；保留所有历史快照。

## 2026-08-08

- 记录免疫治疗 cell2location 已完成 19/19 slices、71,398 spots 和 40-factor 严格合并；因 Stage 13 对 11 个参考不足 inferCNV skip 的验收语义与上游 0 method failures 审计冲突，整体改记为最终验收阻断。
- 更新 CM active v2 Stage 01 为 367/721，第 368 个样本 CopyKAT 活跃；Stage 02–08 仍未开始。
- 更新下载 retry v3 为 90/141 terminal：76 verified、7 skipped_verified、6 unavailable_upstream、1 failed；单列 GSE201347 partial RDS 的本地重试失败。
- 同步 STATUS、三份正式 run record、TODO、runs/artifacts inventories 和本日 daily/Codex 日志。

## 2026-08-06

- 将当前状态从“无运行中任务”纠正为三个 active extensions：四队列免疫治疗整合、恶性/CM 细胞通讯 v2、24-GSE 外部下载重试。
- 新增三份正式运行记录，固定各任务的输入范围、分析流程、方法、实际参数、失败恢复、当前分母和完成门槛。
- 记录免疫治疗分支完成 6/19 张 cell2location 切片，第 7 张已开始 spatial training；最终 19-slice merge 和 acceptance audit 待完成。
- 记录 CM active manifest 因 GSE166555 从 2 个技术组重组为 25 个生物学样本而从 698 变为 721 samples，细胞数保持 1,849,413；Stage 01 当前完成 316/721。
- 记录 CM 旧 attempt 阶段历史不代表 active v2 完成，当前 Stage 02–08 canonical artifacts 不存在，将在 Stage 01 后重算；增加 CopyKAT `genome=hg20` 参数核对任务。
- 记录下载初始 run 的 4 verified/118 failed 及 active 141-task retry v3 的 74 terminal tasks；将 `unavailable_upstream` 与 verified 明确分开。
- 更新状态、路线图、待办和 runs/artifacts inventories；只保存逻辑路径与汇总证据，没有复制大型对象、完整日志或样本级敏感数据。

## 2026-07-31

- 将权威运行 `20260725_120105_unified_pancancer_covarnet` 的状态从“Phase 07 整合完成、注释待进行”更正为“全量多方法整合、证据注释、CoVarNet 和三个固定模块外部验证完成，带审计警告”。
- 新增正式运行记录，固定逐数据集 QC、全量 Harmony/BBKNN/scVI、CellTypist、scIB、证据注释、CoVarNet 和外部验证的参数及证据路径。
- 记录最终发现对象从 human-only 计划的 32 个数据集、2,040,048 个细胞变为 31 个 GSE、1,963,745 个细胞；`GSE278694` 的移除理由尚未在源记录中发现。
- 记录源 `manifests/stage_status.tsv` 与最终产物状态不一致，并建立事件记录；未改写源分析目录的 manifest。
- 更新路线图、任务、数据集/运行/产物清单和统一工作流方法摘要。
- 按数据政策仅同步轻量日志摘要和原日志路径，没有复制大型对象、完整调度日志或细胞级数据。

## 2026-07-26

- 初始化独立的 PanCancer 项目进展与工作日志仓库。
- 确定以 `Project_v3` 为当前主要记录范围，`Project_v1/2` 作为历史背景与证据来源。
- 建立状态、路线图、任务、方法、运行、决策、日志和结构化索引的职责分离。
- 根据当时可见的外部运行证据，将当前状态记录为 Phase 07 统一整合完成、细胞注释待进行。
- 建立只报告、不自动修改的仓库验证脚本和日志生成脚本。
