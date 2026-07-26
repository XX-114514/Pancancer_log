# 记录维护工作流

## 主流程

```text
实际项目产生证据
→ Codex 做最小范围检查
→ 创建或更新 Codex 任务日志
→ 正式运行时创建 run record
→ 更新状态、方法、决策或索引中的权威位置
→ 执行验证
→ 人工检查 git diff
→ commit
→ 向已确认的 private 远程 push
```

## 一次有意义任务

1. 阅读 `AGENTS.md` 要求的上下文并执行 `git status --short`。
2. 创建 `logs/codex/YYYY/YYYY-MM-DD_task-slug.md`。
3. 在日志中写用户要求、初始状态和计划。
4. 检查实际项目证据，记录命令与文件，但不复制大型输出。
5. 完成任务，并在日志中区分事实、解释、假设和建议。
6. 只更新真正发生变化的权威文档：
   - 当前状态变化 → `STATUS.md`
   - 阶段/依赖变化 → `ROADMAP.md`
   - 可执行动作变化 → `TODO.md`
   - 方法变化 → `methods/`
   - 重要选择 → `decisions/`
   - 正式运行 → `runs/` 与 `inventories/runs.tsv`
7. 执行验证，更新任务日志中的结果。
8. 查看完整 diff 后提交；远程明确且为 private 时可普通 push。

## 正式运行的记录时点

- **运行前**：分配 run ID，冻结输入、代码 commit、环境、参数、命令和完成标准。
- **运行中**：记录开始时间、scheduler job ID、关键事件和失败，不复制完整 stdout/stderr。
- **运行后**：记录结束时间、退出状态、输出标识、验证结果、错误和解释。
- **重试**：使用新的 run ID 或显式 attempt ID，保留原失败。

## 状态更新原则

`STATUS.md` 只回答“现在是什么状态”。任务过程写在日志，正式运行细节写在 run record，方案理由写在 ADR。不得把旧的阶段性结果留在当前状态中造成误导。
