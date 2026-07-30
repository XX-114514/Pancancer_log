# Codex 日志：审计并同步统一 CoVarNet 运行

- 日期：2026-07-31
- 任务：报告当前项目进展，提取此前运行代码的处理和参数，并同步轻量记录到项目记录仓库
- 关联运行：`20260725_120105_unified_pancancer_covarnet`

## 用户请求

说明当前项目进展，重点记录已运行代码做了哪些处理、使用了哪些参数，并把处理和日志同步到项目记录仓库。

## 执行计划

1. 盘点运行目录、阶段状态、进程和产物完整性。
2. 从脚本、配置、命令记录和日志提取处理流程与参数。
3. 生成符合记录仓库政策的轻量审计包。
4. 同步文件并运行仓库验证。

## 检查范围

### 源运行

- `RUN_RECORD.md`、`INPUT_AND_CODE_AUDIT.md`
- `manifests/stage_status.tsv`
- Phase 06/07 状态表、QC worker、concat 和 preprocessing scripts
- `config/full_hvg4000_multimethod_v1.json`
- full integration 的 stage logs、SLURM stdout/stderr 和 direct recovery logs
- annotation v1/v2 records、tables、reports 和 final audit
- CoVarNet parameters、NMF/network tables、module reports 和 recovery records
- 三个 validation run 的记录与 final audits
- Python/R session information

### 记录仓库

- `AGENTS.md`
- `STATUS.md`、`ROADMAP.md`、`TODO.md`、`CHANGELOG.md`
- data policy、workflow、reproducibility 文档
- methods、inventories、recent daily/Codex logs
- Git branch、worktree 和 remote 状态

## 主要观察事实

- 当前没有相关运行中进程。
- 源 stage manifest 滞后于最终产物。
- QC 33/33 完成，2,062,230 个细胞通过。
- human-only 计划为 32 datasets / 2,040,048 cells；最终 concat 为 31 GSE / 1,963,745 cells。
- `GSE278694` 的过滤操作可见，但科学理由未记录。
- Harmony、BBKNN、scVI 和 CellTypist 全量运行完成；scIB finalization 成功。
- annotation v2 最终审计 `PASS`：17 major lineages、136 clusters、106 annotations。
- CoVarNet K9 完成；R abundance 和 network figure 步骤均有已记录 fallback。
- 三个固定模块外部验证审计均通过，但 BH 后无显著模块。

## 关键参数摘要

- QC：genes 200/500、counts 500/1000、MAD 3、mt 15–25%、Scrublet rate 0.04、max removal 0.08、30 PCs、seed 0。
- integration：10k/log1p、4,000 HVG、50 PCs、Harmony 20 iterations、30 neighbors、BBKNN 2/batch + trim 200、scVI latent 30 / 2 layers / hidden 128 / 30 epochs / batch 2048、seed 20260727。
- annotation：L1 resolution 0.1、lineage resolution 0.4–1.0、small cluster threshold 30、chunk 10,000、marker panel 194、seed 20260729。
- CoVarNet：min cells 100、min-max、Pearson、rank 2–20、K=9、nrun 30、corr 0.2、FDR 0.05。
- validation：frozen mapper/W、10k/log1p、chunk 5,000、confidence 0.50、soft primary、hard sensitivity、NNLS projection、paired Wilcoxon + BH。

## 失败与恢复

- 初次 sparse-Zarr validator、33-dataset harmonization、capped HVG/Harmony 和 full HVG job 均失败后由修正版/retry 成功。
- scVI 训练完成后因对 dict 类型 history 调用 `.to_csv` 失败；r2 direct resume 复用模型成功。
- CoVarNet `sc_cm_recover` 失败后从 `coef(nmf_final)` 恢复 abundance。
- RColorBrewer 调色板不足导致 R 网络图失败；Python/networkx 基于既有网络表重绘成功。
- locale warning 无计算影响。

## 工具与环境问题

- 默认 sandbox 多次出现 `bwrap: Creating new namespace failed`；对必要的只读审计和用户明确要求的同步操作使用受控权限执行。
- 基础环境没有可用 `jq`，并且 base Python 不含 pandas；改用已有文本/结构化记录和已验证的项目环境信息。
- 部分初次文件匹配路径不准确或输出过长，随后以限定深度的 `find` 和精确路径重查。
- 本地 Git 版本不支持某些较新的 branch 查询选项，改用兼容命令。

## 解释与判断

- “当前进展”采用最终产物与 audit，而不是滞后的 source manifest。
- `GSE278694` 的排除原因没有证据，故不推断理由，只记录为治理缺口。
- scIB 指标存在维度权衡，不宣称单一整合方法全面最佳。
- 外部验证缺少 BH 后显著性和直接临床 endpoint，不上升为疗效或生存结论。

## 本次记录变更

- 新增正式 run record。
- 更新 `STATUS.md`、`ROADMAP.md`、`TODO.md`、`CHANGELOG.md`。
- 更新统一方法摘要和 runs/artifacts/datasets inventories。
- 新增 daily log、本 Codex log 和记录漂移 incident。
- 未复制完整日志、大型 H5AD/Zarr、细胞级表或完整样本 metadata。
- 未改写源运行目录的 stage manifest。

## 验证

同步后执行：

- `git status --short`
- `git diff --stat`
- `git diff --check`
- 仓库内置验证脚本
- TSV 列数一致性与 Markdown 相对链接检查
- 敏感绝对路径和超大文件检查

验证结果：仓库内置验证脚本报告 `0 error(s), 0 warning(s)`；`git diff --check` 通过，三个 TSV 列数一致，没有敏感绝对路径、超 5 MiB 文件、禁止数据格式或断裂的本地 Markdown 链接。远程 push 因既有认证阻塞不在本次范围内强行执行。

## 未决事项

- `GSE278694` 排除理由。
- source stage manifest 对齐。
- low-confidence / unresolved / ambiguous annotation 人工复核。
- `GSE162498` extreme-scale 分支决定。
- 恶性/CNV 分支完成状态。
- 带可靠临床 endpoint 的独立验证。
