# 事件：统一运行的记录状态漂移与未说明范围变化

- 日期：2026-07-31
- 影响运行：`20260725_120105_unified_pancancer_covarnet`
- 类型：记录治理 / 复现边界
- 严重度：中
- 当前状态：开放，已完成记录层纠偏

## 摘要

运行目录的最终产物和审计已经覆盖全量整合、证据注释、CoVarNet 和三个外部验证，但源 `manifests/stage_status.tsv` 仍把后续阶段标为 `pending`。此外，human-only 计划中的 `GSE278694` 在最终 concat 前被移除，现有代码和日志能证明移除发生，但没有发现科学理由或正式决定记录。

## 已观察事实

- 2026-07-30 11:18 前，最终 integration、annotation、CoVarNet 和 validation 产物已写出。
- annotation `FINAL_DELIVERY_AUDIT.json` 为 `PASS`。
- 三个 validation `FINAL_AUDIT.json` 均通过。
- 最终 concat 为 31 个 GSE、1,963,745 个细胞和 773 个样本。
- 之前的 human-only plan 为 32 个数据集、2,040,048 个细胞和 794 个样本。
- 过滤 notebook/代码和 “shard intentionally deleted” 日志指向 `GSE278694` 被移除。
- 源 stage manifest 和旧 run index 仍停留在 annotation pending。

## 解释

这是“记录状态漂移”和“范围决定未文档化”，不是已证实的数值计算错误。当前最终对象、注释、模块和投影审计仍可读取，但最终 discovery cohort 的纳入标准不能在缺少 `GSE278694` 决定依据时被称为完全冻结。

## 影响

- 自动读取 stage manifest 的报告会低估项目进度。
- 不知道 `GSE278694` 的移除理由，会阻碍严格复现和敏感性解释。
- 旧 capped branch 的 32-dataset 统计与最终 31-GSE 全量对象不能混写。
- source commit 未冻结进一步降低精确复现能力。

## 已采取措施

- 以实际产物、脚本、日志和最终 audit 为准，新增正式运行记录。
- 更新状态、路线图和结构化清单，明确区分 32-dataset planning snapshot 与 31-GSE final concat。
- 保留源 manifest 和原运行目录不变，没有删除或覆盖历史证据。
- 完整日志仍留在 `${RUN_ROOT}`，记录仓库只保留路径和必要错误摘要。

## 待解决

1. 由原执行者或项目负责人说明 `GSE278694` 的排除理由、影响和永久/临时状态，并形成 ADR。
2. 在保留历史失败/恢复链的前提下更新源 `stage_status.tsv`。
3. 为未来 production run 强制冻结 source commit、config checksum 和 active-output manifest。

## 相关记录

- [正式运行记录](../../runs/20260725_120105_unified_pancancer_covarnet.md)
- [当前状态](../../STATUS.md)
- [统一工作流](../../methods/project-v3-unified-workflow.md)
