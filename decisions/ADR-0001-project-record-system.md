# ADR-0001: 建立独立的项目权威记录系统

- Date: 2026-07-26
- Status: Accepted
- Related stage: 全项目

## Context

PanCancer 工作区已有 `Project_v1`、`Project_v2` 和 `Project_v3`，包含代码、数据、运行、历史文档和大型分析产物。上级目录不是 Git 仓库，内部另有两个嵌套 Git 仓库。需要长期记录当前状态、任务、方法、正式运行和技术决策，同时避免误跟踪原始数据或覆盖现有文档。

## Problem

零散 Markdown 无法清楚区分当前事实、时间历史、方法、运行和决策，也难以持续审计 Codex 的工作。直接把整个分析工作区初始化为单一 Git 仓库则会增加大型文件、敏感信息和嵌套仓库误入版本控制的风险。

## Considered options

1. 在 PanCancer 根目录直接初始化 Git 并重构现有文档。
2. 继续在各 Project 目录中分别增加零散日志。
3. 建立独立 `project-records/` Git 仓库，只保存轻量级记录和索引。

## Decision

采用选项 3：

- `project-records/` 为独立 Git 仓库；
- 当前主要范围为 `Project_v3`，`Project_v1/2` 作为历史背景与证据来源；
- 状态、日志、方法、运行、决策和 inventories 分离；
- Codex 必须为每次有意义任务记录自身检查、修改和验证；
- 数据和大型产物留在分析工作区，只用逻辑路径和稳定 ID 引用；
- GitHub 目标为 private。

常规 Git commit/push 的一般安全默认是需要用户确认。本仓库用户已于 2026-07-26 给出长期明确授权：Codex 可对本独立记录仓库自动 commit，并在 private 远程已经明确配置后普通 push。该授权不包括强制推送、删除、改写历史或操作上级/嵌套仓库。

## Rationale

独立仓库能隔离大型分析内容和嵌套 Git 历史；职责分离使当前状态不会与时间流水账混淆；强制 Codex 日志让自动化工作可追溯；运行记录与 inventories 连接实际证据而不复制数据。

## Consequences

- 需要维护逻辑路径与真实工作区的映射。
- 同一事实必须确定权威位置，减少重复。
- 正式运行需要额外创建轻量级 run record 和索引。
- 推送前仍需敏感信息扫描和 diff 审查。
- 尚未配置远程时只能本地提交，不能猜测 GitHub 地址。

## Validation

- 初始化后运行 `scripts/validate_repository.sh`。
- 检查 Git 只跟踪 `project-records/` 内的轻量级文件。
- 检查上级 `README.md`、`AGENTS.md` 与分析目录未被修改。

## Conditions for revisiting

- 记录范围扩展为多个独立大项目；
- GitHub 可见性改变；
- 需要对结构化元数据使用数据库而非 TSV；
- 现有记录规模导致 Markdown/TSV 无法可靠维护；
- 安全策略或 commit/push 授权发生变化。
