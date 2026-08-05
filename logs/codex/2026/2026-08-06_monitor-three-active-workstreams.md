# Codex 日志：监控三个 active workstreams

- 日期：2026-08-06
- 任务：监控免疫治疗整合、CM 细胞通讯和新数据下载，并把进度、流程、方法与参数写入记录仓库
- 关联运行：`immunotherapy_4cohort_20260803_v1`、`20260731_copykat_infercnv_intersection_v1`、`20260806_025419_slurm_retry_v3`

## 用户要求

监控三个正在进行的任务，记录当前进度，并按既有要求整理具体分析流程、方法和参数到 `project-records`。

## 开始前项目状态

- 记录仓库 branch `main`，工作树 clean，最近提交 `ff51d0e`。
- `STATUS.md` 最近更新为 2026-07-31，仍写“当前无运行中任务”。
- 主统一运行记录完整，但三个后续扩展尚无独立正式运行记录。

## 初始计划

1. 读取记录仓库规则、状态、路线图、待办、方法、最近日志与结构化清单。
2. 用进程、Slurm、stage markers、append-only logs 和文件双时点增长核验三个任务的实时状态。
3. 从实际入口脚本、审计 JSON/Markdown 和命令行提取分析流程、方法和参数。
4. 新增三份 run records，并更新权威状态和 inventories。
5. 检查 diff、TSV 列数、链接、敏感路径和仓库验证脚本。

## 使用的项目技能规范

- `pan-cancer-covarnet-cell2location-spatial-projection`：用于检查 cell2location 必须保留完整 spot × factor abundance、q05 数值合同和跨癌种解释边界。
- `pan-cancer-covarnet-translational-analysis`：用于区分直接临床字段、derived/proxy interpretation 和不能支持的疗效/生存结论。
- `pan-cancer-geo-supplement-audit`：用于要求 GEO/SRA 上游状态、官方 metadata、失败与 `unavailable_upstream` 分开记录，不自动改本地注释。

这些规范影响了记录中的验证门槛与解释边界；未用它们修改正在运行的分析。

## 检查过的记录仓库文件

- `AGENTS.md`、`README.md`、`STATUS.md`、`ROADMAP.md`、`TODO.md`、`CHANGELOG.md`。
- `docs/WORKFLOW.md`、`docs/REPRODUCIBILITY.md`、`docs/DATA_POLICY.md`。
- `methods/project-v3-unified-workflow.md`。
- 最近的 daily、Codex、incident 日志。
- `runs/20260725_120105_unified_pancancer_covarnet.md`。
- `inventories/runs.tsv`、`inventories/artifacts.tsv`、`inventories/datasets.tsv`。

## 检查过的外部证据

### 免疫治疗四队列

- 运行审计、stage events、Phase 5/6 runners、cell2location v1/v2/v3 实现。
- sc/snRNA、annotation、malignancy、spatial、reference 和 regression audits。
- per-slide `COMPLETED.json`、实时 epoch 输出和进程资源状态。

### CM 细胞通讯

- 后台编排、active manifest、GSE exclusion/regroup audits、stage history 和实时 sample CNV log。
- inferCNV、CopyKAT、scMalignantFinder、CoVarNet、LIANA、NicheNet 与 final audit scripts/parameters。
- 当前 Stage 02–08 expected artifacts 是否存在。

### 新数据下载

- decision/file/recovery manifests、initial summary、retry progress、events 和 Slurm stdout/stderr。
- 下载与 retry runners、进程/Slurm 资源状态。
- GSE274934 RAW tar 的 20 秒双时点大小和 heartbeat。

## 执行过的命令类型

```bash
git status --short
git symbolic-ref --short HEAD
git log -5 --oneline
pgrep -af <task patterns>
ps -p <task pids> -o <status fields>
squeue -u <user>
scontrol show job <job_id>
find <scoped directories> -type f
rg -n <method and parameter patterns> <scoped scripts/logs>
head / tail / sed / wc / awk / stat
sleep 20  # only for one read-only file-growth confirmation
bash scripts/validate_repository.sh
```

命令记录使用逻辑结构；内部绝对路径和完整大型日志不进入记录仓库。

## 关键发现

### 免疫治疗

- 1,518,479 total observations；sc/snRNA singlet integration 1,440,403；Visium 71,398 spots。
- 43 final sc/snRNA labels；GSE316195 strict malignant intersection 5,648 cells。
- cell2location reference 28,357 cells、16,544 signature genes、40 factors；RegressionModel 250 epochs。
- spatial model 1,000 epochs、N_cells_per_location 30、detection_alpha 200、posterior 1,000、12 torch threads；当前 6/19 completed；第 7 张已开始训练。

### CM 通讯

- active v2 721 samples / 1,849,413 cells；316/721 sequence completed。
- GSE166555 从 2 个技术组变为 25 个生物学样本，细胞覆盖保持 48,819。
- inferCNV window 250/reference q99；CopyKAT `ngene.chr=5`、`win.size=25`、`KS.cut=0.1`、Euclidean、2 cores；scMalignant 10k/no-log1p。
- CoVarNet K9、rank 2–20、nrun 30、corr 0.2、FDR 0.05；LIANA top-members 10、min group 20、expr-prop 0.10；NicheNet 使用本地 v2 model。

### 下载

- initial 122 tasks：4 verified、118 failed，主要为 connection refused。
- retry v3 141 tasks：74 terminal；61 verified、7 skipped_verified、6 unavailable_upstream。
- 单 worker、12 秒最小请求间隔、300 秒 cooldown/timeout、8 MiB chunk；GSE274934 当前文件持续增长。

## 修改过的文件及原因

- `runs/20260803_immunotherapy_4cohort_analysis.md`：正式记录四队列流程、参数、进度和解释边界。
- `runs/20260731_malignancy_communication.md`：正式记录 active v2 分母、逐样本恶性、CoVarNet/LIANA/NicheNet 和当前进度。
- `runs/20260805_external_geo_download.md`：正式记录 initial/retry 下载范围、参数、状态和完整性门槛。
- `STATUS.md`：纠正“当前无任务”，列出三个 active extensions。
- `ROADMAP.md`、`TODO.md`、`CHANGELOG.md`：同步阶段、完成标准、可执行动作与重要变化。
- `inventories/runs.tsv`、`inventories/artifacts.tsv`：增加结构化运行与产物索引。
- `logs/daily/2026/2026-08-06.md` 和本日志：记录本次工作与证据链。

## 失败尝试与恢复

- 默认 sandbox 出现 bwrap namespace failure；必要的只读审计和目标仓库写入改用受控权限。
- 本机旧 Git 不支持 `git branch --show-current`；改用兼容命令。
- 首轮宽范围输出被截断；改为精确文件、terminal markers 和结构化计数。
- 补丁工具读取 /tmp 既有文件时受 sandbox helper 影响；使用补丁生成完整 replacement，再在 /tmp 暂存区复制为待同步文件。没有删除任何文件。

## 未验证内容

- 三项任务的最终结果，因为它们仍在运行。
- 免疫治疗 19-slice merge 与 final acceptance。
- CM active v2 下游重算和 CopyKAT `hg20` 意图。
- 下载 retry final summary 与全部文件完整性。

## 尚未解决的问题

- 三项 active tasks 的 terminal status。
- 数据盘容量是否足以完成所有大型下载。
- GitHub 认证阻塞。
- 既有 `GSE278694` 决定、source stage drift 和注释人工复核缺口。

## 建议的下一步

- 在三个任务出现关键阶段变化、失败或完成时更新同一组 run records 和新快照日志。
- 只依据 final audit/summary 更新 completed 状态，不依据 PID 消失或日志停止单独判断。

## 建议的 Git commit message

```text
docs: record three active pancancer workstreams
```

## 验证方式和结果

同步后执行：

- `git status --short`、`git diff --stat`、`git diff --check`。
- `bash scripts/validate_repository.sh`。
- runs/artifacts/datasets TSV 列数一致性检查。
- 新增文件的凭据模式、敏感绝对路径、大文件和相对 Markdown 链接检查。
- 人工查看本次 tracked diff，并逐份读取新 run records 的关键状态和参数段。

结果：内置验证脚本报告 `0 error(s), 0 warning(s)`；`git diff --check` 通过；三个 inventories TSV 列数一致；未发现敏感绝对路径、明显凭据、大于 5 MiB 的候选文件、禁止数据格式或断裂的本地相对链接。locale warning 来自服务器环境，不影响验证结论。
