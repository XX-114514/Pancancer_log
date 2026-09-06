# PanCancer 项目进展与工作日志

本仓库是 PanCancer 项目的轻量级、可审计记录系统，当前以 `Project_v3` 的泛癌单细胞转录组主线为主要范围。它记录当前状态、任务历史、方法、正式运行、技术决策及证据索引，是项目进展的权威入口，但不是分析数据仓库。
## 2026-09-06 权威入口

- 当前细胞身份/CNV 恶性参考：[V7 冻结 release](releases/annotation_v7_20260905/README.md)。
- 当前论文通讯证据与 Prism 成图：[20260906 evidence chain](releases/communication_evidence_chain_20260906/README.md)。
- Git 维护：[维护工作流](docs/GIT_MAINTENANCE.md)。
- Prism 与 supplementary：[交付指南](docs/PRISM_SUPPLEMENTARY_GUIDE.md)。
- 版本边界：[版本政策](docs/VERSION_POLICY.md)。


## 不存储的内容

本仓库不保存原始测序数据、患者身份信息、受控访问数据、H5AD/RDS/矩阵对象、完整作业日志、缓存、压缩归档或大型分析产物。真实文件留在分析工作区，本仓库只记录稳定标识符、逻辑路径、校验信息和验证结论。详见 [数据政策](docs/DATA_POLICY.md)。

## 快速查看当前状态

1. 阅读 [STATUS.md](STATUS.md) 获取当前阶段、阻塞和下一步。
2. 阅读 [ROADMAP.md](ROADMAP.md) 获取阶段与里程碑。
3. 阅读 [TODO.md](TODO.md) 获取可执行任务。
4. 通过 [运行索引](inventories/runs.tsv) 找到正式运行记录及外部证据。
5. 通过 [术语表](docs/GLOSSARY.md) 和 [命名规范](docs/NAMING_CONVENTIONS.md) 理解标识符。

## 日常使用

```bash
# 创建当天日志
bash scripts/new_daily_log.sh

# 创建 Codex 任务日志
bash scripts/new_codex_log.sh "review phase 07 annotation"

# 创建正式运行记录
bash scripts/new_run_record.sh "phase07 annotation"

# 输出轻量级项目快照；结果目录参数可选
bash scripts/project_snapshot.sh "${RESULTS_ROOT}"

# 只检查、不修改
bash scripts/validate_repository.sh
```

新日志创建后必须补充实际内容，不能保留未替换占位符。正式运行还应同步更新 `inventories/runs.tsv`；新数据集或重要产物分别更新 `inventories/datasets.tsv` 和 `inventories/artifacts.tsv`。

## Git 工作流

```bash
git status --short
git diff --stat
git diff
bash scripts/validate_repository.sh
git add <reviewed-files>
git commit -m "type: concise description"
git push
```

本仓库当前获得了用户对常规 `commit` 和向已确认 private 远程执行普通 `push` 的长期授权；尚未配置远程时不得猜测地址。强制推送、改写历史、删除文件和破坏性 Git 操作始终禁止，详见 [AGENTS.md](AGENTS.md)。

## 文档职责

| 权威位置 | 只负责 |
| --- | --- |
| `STATUS.md` | 当前真实状态 |
| `ROADMAP.md` | 阶段、里程碑、依赖和完成标准 |
| `TODO.md` | 可执行任务 |
| `CHANGELOG.md` | 重要变化 |
| `logs/daily/` | 当天发生的事情 |
| `logs/codex/` | Codex 每次任务的检查、修改和验证 |
| `runs/` | 正式运行的输入、参数、命令、输出和结果 |
| `methods/` | 当前认可的方法 |
| `decisions/` | 重要方案选择及理由 |
| `inventories/` | 数据集、运行和产物的结构化索引 |

相同事实只在一个权威位置详细维护，其他文档通过相对链接引用。完整关系见 [目录结构](docs/DIRECTORY_STRUCTURE.md) 和 [维护工作流](docs/WORKFLOW.md)。
