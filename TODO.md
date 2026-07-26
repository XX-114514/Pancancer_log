# 可执行任务

| ID | Priority | Status | Task | Dependency | Completion criterion | Related files |
| --- | --- | --- | --- | --- | --- | --- |
| REC-001 | P0 | in_progress | 配置 GitHub private 远程并核验目标 | GitHub 仓库地址 | `git remote -v` 指向确认的 private 仓库，普通 push 成功 | `docs/DATA_POLICY.md` |
| ANN-001 | P0 | todo | 冻结 Phase 07 annotation 输入与注释策略 | P3 integrated model | run record 写明对象标识、commit、环境、marker/reference 和完成标准 | `STATUS.md`, `runs/` |
| ANN-002 | P0 | todo | 执行并验证主类群注释 | ANN-001 | 标签列、证据表、UMAP、计数表及人工复核结果齐全 | `methods/project-v3-unified-workflow.md` |
| EXT-001 | P1 | todo | 决定 `GSE162498` 的 extreme-scale 分支处理 | 资源与分块策略复核 | 形成 ADR 或正式 run，并明确是否进入主整合 | `ROADMAP.md` |
| INV-001 | P1 | todo | 补全 32 个 human integration 数据集的逐项记录 | 冻结 manifest | 每行含 dataset ID、角色、状态和权威证据 | `inventories/datasets.tsv` |
| REP-001 | P2 | todo | 建立每次有意义任务的日志习惯 | 无 | 新任务同时具有 Codex log；正式运行另有 run record | `AGENTS.md`, `templates/` |
| SEC-001 | P1 | todo | 首次 push 前进行人工敏感信息审查 | REC-001 | 验证脚本通过且人工确认 diff 不含敏感路径/数据 | `docs/DATA_POLICY.md` |

状态建议使用：`todo`、`in_progress`、`blocked`、`done`。任务完成后保留该行并链接完成证据；定期将历史性已完成任务汇总到 `CHANGELOG.md`。
