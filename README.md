# PanCancer 项目进展与工作日志

## 最新公共下载检查

2026-10-11 18:56—19:02 UTC+8：[下载/QC/验收及下一步报告](inventories/runtime/2026-10-11/DOWNLOAD_READINESS.md)。5,341,204源细胞；4份QC分片993,312细胞；正式新增验收0。下载进程已退出，第二批8/14文件校验，未完成文件不计矩阵到库。逐来源QC/验收可继续，全量整合未就绪。下方11:22入口和动态描述保留为历史快照。


## 最新CM探索节点与时间线

2026-10-11 19:27:55 UTC+8：[全部65空间端点预检完成](runs/20261011_spatial_contact_native_preflight_v3.md)。正确作者181,056标签全匹配，6口径复算通过；4种LR在作者粗谱系口径有描述性覆盖，尚待邻近检验。失败和方向错误保留、显式更正。

2026-10-11 18:52:59 UTC+8：[空间输入修复完成](runs/20261011_spatial_native_input_repair_v1.md)；保留失败及旧对象，新对象原CSV全矩阵独立核对通过。下一阶段全部13×5候选端点预检；通讯/CM机制仍未验证。

2026-10-11 18:36:08 UTC+8：[空间输入检查与CSV索引故障](logs/incidents/2026-10-11_spatial-csv-index-column.md)。5份文件只读汇总完成；TMA2暂停归一化/空间推断，先修复并审计；患者映射/背景/邻近尚待验证。详见[持续任务日志](logs/codex/2026/2026-10-11_milestone-log-and-cloud-sync.md)。

2026-10-11 18:20:36 UTC+8：[关键节点证据](inventories/analysis/2026-10-11/README.md)、[完整时间线](inventories/analysis/2026-10-11/task_timeline_v1.tsv)、[正式run](runs/20261011_family_candidate_validation_v1.md)。336条参考支持记录完成检查，自动15条后续验证记录；身份内部程序检查完成，空间接触和机制尚待验证。关键/探索节点的记录与普通push规则见[RUNTIME_SYNC](docs/RUNTIME_SYNC.md)。

## 2026-10-11 当前入口：服务器运行与Codex云端开发

- **先看下载与Project_v3动态**：[本次实测报告](inventories/runtime/2026-10-11/REPORT.md)，采集于2026-10-11 11:22:53（UTC+8）。GitHub保存快照，不是实时作业面板。
- **当前状态/下一步**：[STATUS](STATUS.md)、[任务板](TODO.md)、[阶段路线](ROADMAP.md)。
- **开始云端开发**：[Codex Cloud配置与工作流](docs/CODEX_CLOUD_DEVELOPMENT.md)，[可直接提交的首个任务](prompts/cloud_next_task.md)。
- **检查代码**：Python ≥3.10环境执行 `bash scripts/cloud_check.sh`；状态采集与同步见[操作说明](docs/RUNTIME_SYNC.md)。
- 新公共源数据5,341,204细胞；608,684细胞已有QC分片；最终新增验收尚未登记。Project_v3表达标签及新CoVarNet/LIANA计算已完成，最终科学细身份仍待GSE299623审阅。
- 目标仓库已核实为**public**。2026-10-11用户明确授权向此仓库提交与普通push；仅上传已审查的代码、汇总与脱敏记录。下文此前版本入口保留为历史依据，不覆盖本次状态。

本仓库是 PanCancer 项目的轻量级、可审计记录系统，当前以 `Project_v3` 的泛癌单细胞转录组主线为主要范围。它记录当前状态、任务历史、方法、正式运行、技术决策及证据索引，是项目进展的权威入口，但不是分析数据仓库。

## 2026-09-22 投稿冻结入口

- Nature 投稿准备与注释上游锁：[`manuscript_upstream_freeze_nature_20260922`](releases/manuscript_upstream_freeze_nature_20260922/README.md)。
- V8 最新谱系审计限定：[`annotation_v8_lineage_audit_20260917`](releases/annotation_v8_lineage_audit_20260917/README.md)。
- V9 最新候选报告与开放复核表：[`annotation_v9_candidate_20260919`](releases/annotation_v9_candidate_20260919/README.md)。
- 全局投稿参考仍是 V7；V8 仅为 lineage taxonomy，V9 为 `FROZEN_CANDIDATE`，均未晋升为全局参考。
- 当前 V5-derived 通讯/空间/临床图不能无说明地称为 V7/V9-derived。

## 2026-09-12 权威入口

- 当前细胞身份/CNV 恶性参考：[V7 冻结 release](releases/annotation_v7_20260905/README.md)。
- [V8 development software candidate](releases/annotation_v8_development_20260912/README.md)
  是已索引的 `NOT_FROZEN` 全量元数据软件候选；它不替代 V7，也不授权正式 downstream
  变更或生物学完成声明。
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
