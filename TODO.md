# 可执行任务

| ID | Priority | Status | Task | Dependency | Completion criterion | Related files |
| --- | --- | --- | --- | --- | --- | --- |
| REC-001 | P0 | blocked | 核验 GitHub private 远程并完成普通 push | 服务器 GitHub 认证 | `git remote -v` 指向确认的 private 仓库，普通 push 成功 | `docs/DATA_POLICY.md` |
| REC-002 | P0 | done | 审计并同步统一 CoVarNet 运行的处理、参数和日志索引 | 权威运行目录可读 | 正式 run record、状态、方法、清单、Codex log 和 incident 均通过仓库验证 | `runs/20260725_120105_unified_pancancer_covarnet.md` |
| GOV-001 | P0 | todo | 对齐源 `stage_status.tsv` 与最终产物状态 | 最终审计文件 | manifest 不再把已完成阶段标记为 `pending`，且保留历史失败/恢复链 | `logs/incidents/2026-07-31_unified-run-record-drift.md` |
| GOV-002 | P0 | blocked | 补齐 `GSE278694` 从最终整合移除的科学理由 | 原执行者决定或原始证据 | ADR/run record 写明理由、影响和是否永久排除 | `logs/incidents/2026-07-31_unified-run-record-drift.md` |
| ANN-001 | P0 | done | 冻结 Phase 08 annotation 输入与注释策略 | P3 integrated model | run record 写明对象、环境、marker/reference 和完成标准 | `runs/20260725_120105_unified_pancancer_covarnet.md` |
| ANN-002 | P0 | done | 执行并验证主类群和亚型注释 | ANN-001 | 标签列、证据表、UMAP、计数表和最终审计齐全 | `methods/project-v3-unified-workflow.md` |
| ANN-003 | P1 | todo | 人工复核低置信度和争议注释 | ANN-002 | 56 low-confidence、10 unresolved、5 ambiguous clusters 有复核结论 | `STATUS.md` |
| EXT-001 | P1 | todo | 决定 `GSE162498` 的 extreme-scale 分支处理 | 资源与分块策略复核 | 形成 ADR 或正式 run，并明确是否进入主整合 | `ROADMAP.md` |
| CNV-001 | P1 | todo | 核验恶性/CNV 分支是否已执行及其证据 | P4 annotations | 形成独立 run record，或明确本次范围不含该分支 | `ROADMAP.md` |
| VAL-001 | P1 | done | 对三个免疫治疗队列执行固定 9 模块投影 | 冻结 mapper 和 W | 三个队列 `FINAL_AUDIT.json` 均通过 | `runs/20260725_120105_unified_pancancer_covarnet.md` |
| VAL-002 | P1 | todo | 选择含可靠临床 endpoint 的独立验证队列 | 数据许可与正式 metadata | 预先冻结 endpoint、配对规则和统计模型并完成验证 | `STATUS.md` |
| INV-001 | P1 | todo | 补全最终 31 个 discovery GSE 的逐项记录 | 冻结 final concat manifest | 每行含 dataset ID、角色、状态和权威证据 | `inventories/datasets.tsv` |
| REP-001 | P2 | in_progress | 持续为每次有意义任务写入 Codex log | 无 | 新任务均有日志；正式运行另有 run record | `AGENTS.md`, `templates/` |
| SEC-001 | P1 | todo | push 前进行人工敏感信息审查 | REC-001 | 验证脚本通过且人工确认 diff 不含敏感路径/数据 | `docs/DATA_POLICY.md` |

状态使用：`todo`、`in_progress`、`blocked`、`done`。已完成任务保留原行并链接证据。
