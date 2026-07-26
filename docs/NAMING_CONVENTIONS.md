# 命名与运行编号规范

## 通用规则

- 日期：`YYYY-MM-DD`
- 时间戳：`YYYYMMDD_HHMMSS`，默认时区为 `Asia/Shanghai`
- 稳定标识符使用小写英文、数字和连字符，便于搜索和跨平台处理。
- 正式 accession、dataset ID、sample ID、run ID、output ID 和阶段名不得仅为美观而修改。
- 已有名称含义不清时，在本文或 `GLOSSARY.md` 解释，不擅自重命名历史文件。
- 路径文档化时优先用逻辑变量，避免真实用户名和服务器绝对路径。

## 文件命名

| 类型 | 格式 | 示例 |
| --- | --- | --- |
| 每日日志 | `logs/daily/YYYY/YYYY-MM-DD.md` | `logs/daily/2026/2026-07-26.md` |
| Codex 日志 | `logs/codex/YYYY/YYYY-MM-DD_task-slug.md` | `logs/codex/2026/2026-07-26_review-phase07.md` |
| 运行记录 | `runs/YYYY/YYYYMMDD_HHMMSS_run-name.md` | `runs/2026/20260726_103000_phase07-annotation.md` |
| 技术决策 | `decisions/ADR-NNNN-short-title.md` | `decisions/ADR-0001-project-record-system.md` |
| 故障记录 | `logs/incidents/YYYY/YYYY-MM-DD_short-title.md` | `logs/incidents/2026/2026-07-26_annotation-failure.md` |

## Slug 规则

- 转为小写；
- 空格和下划线转为单个连字符；
- 移除除 `a-z`、`0-9`、连字符之外的字符；
- 避免含患者名、用户名、内部主机名或秘密；
- 名称无法安全转换时使用通用 slug，并在正文写原任务名称。

## Run ID

默认格式：

```text
YYYYMMDD_HHMMSS_run-name
```

Run ID 一旦用于正式记录即保持稳定。重试推荐追加新时间戳或 `attemptNN`，不要覆盖旧 run。dry-run、smoke、正式运行应在 Metadata 中显式区分。

## 阶段名

当前统一主线使用 `P0`–`P7` 作为路线图层级；外部 `Project_v3` 运行中的 `Phase 06`、`Phase 07` 保留原正式名称。两者是不同层次，不应强行重命名：

- 路线图 `P2` 对应当前主线的逐数据集 QC；
- 路线图 `P3` 对应外部 `Phase 07` 的统一整合；
- 外部阶段编号仍以原 run record 为准。
