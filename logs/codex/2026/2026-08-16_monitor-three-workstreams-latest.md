# Codex 任务日志：三项任务最新监控、允许边界复核与提交

## 用户要求

重新检查免疫治疗整合、CM 细胞通讯和新数据下载三项任务；回顾此前允许情况，输出最新状态和中断/卡住原因，更新项目记录并 Git 提交。

## 开始前状态

- `project-records` 当前 HEAD 为 `ce56fc8`，此前已记录 CM Stage01 在第 581 个样本因 group hard quota 耗尽而中断。
- `STATUS.md` 存在用户所有的未提交首行空白改动；本任务保持不修改、不暂存、不提交。
- 三项任务的最近完整联合快照为 2026-08-12；CM 另有 2026-08-16 03:37 容量审计。

## 初始计划

1. 回顾三项任务的允许边界、最近 run/incident 和 Git 状态。
2. 同时检查进程、调度器、日志、terminal marker、final audit 和关键文件大小。
3. 用任务自身的完成合同重新验证 CM 断点，并重新查询文件系统与 group quota。
4. 仅以追加方式更新权威 run、incident、daily、roadmap/TODO/changelog 和本日志。
5. 验证完整 diff，显式排除 `STATUS.md` 后普通 commit。

## 回顾的允许边界

- 免疫治疗：19-slice cell2location 只能按严格 abundance/spot 合同保留；四队列完成仍要求 final audit。尚未发现将 reference-insufficient inferCNV skip 认定为允许终态的书面决定。
- CM：只允许按生物学样本、raw counts 运行 CNV；整合对象只作结果载体。恢复必须用新 recovery ID 和 append-only attempt，保留旧 partial；第 581 个原生 CopyKAT prediction 只有通过 cell-ID 一对一与来源 hash 审计后才可复用。发布前仍需确认 `genome=hg20`。
- 下载：`verified`、`skipped_verified`、`unavailable_upstream` 和 `failed` 必须分开；partial 不得记为 verified。GSE201347 只能用新的版本化恢复尝试，旧 partial 和 retry summary 必须保留。
- 工作区 AGENTS 规则禁止删除；本任务没有执行任何删除、覆盖历史或破坏性 Git 操作。

## 检查的证据

- 免疫：Stage 13 attempt、失败 marker/log、19 个 cell2location completion markers、merge audit、malignancy audit、final audit 路径和新文件时间窗。
- CM：active v2 manifest、recovery log、第 581 个 attempt、runner `validate_completed_attempt()`、下游 canonical audit、进程/tmux/Slurm、XFS group quota、filesystem/inode。
- 下载：retry v3 final summary、stdout/stderr、GSE201347 partial 的 mtime/size、下载树大小、新文件时间窗和进程/Slurm。
- 记录仓库：AGENTS、README、STATUS、ROADMAP、TODO、CHANGELOG、工作流、数据政策、三份 run、相关 Codex log 和 CM incidents。

## 执行过的检查

- `git status --short`、最近提交和 `STATUS.md` 单文件 diff。
- `pgrep -af`、`squeue -u`、`tmux list-sessions`。
- `find`/`rg --files`/`stat`/`tail`/JSON 格式化读取关键状态、日志和产物。
- 在 `scanpy` 环境中只读导入 CM runner，并对 721 个样本调用其 completion contract。
- `df -h`、`df -i` 与 XFS group quota 精确查询；只做数值换算，不修改运行目录。

## 关键发现

### 已观察事实

- 冻结时间：2026-08-16 07:00 +08:00。未发现三个任务的匹配活动进程、Slurm 作业或 tmux 会话；任务目录在各自最近快照后无新写入。
- 免疫：19/19 cell2location completion markers；71,398 spots、40 factors 严格 merge 继续通过。Stage 13 仍只有失败的 attempt 01，final audit 不存在。
- CM：strict valid total/prefix 均为 580/721；首个未完成仍是第 581 个样本，剩余 141 个样本、501,318 cells。Stage 02–08 和 final audit 均未开始。
- CM 最新 quota 为 1,466,340,820/1,572,864,000 blocks，余 106,523,180 blocks，约 101.588 GiB；文件系统约有 871 GiB available、inode 使用约 4%。
- 下载：retry v3 保持 141/141 terminal，127 verified、7 skipped_verified、6 unavailable_upstream、1 failed。GSE201347 partial 保持 14,413,317,040/19,071,156,087 bytes（75.577%），尚缺约 4.338 GiB。

### 中断或卡住原因

- 免疫是验收策略阻断：11 个 reference-insufficient inferCNV skip 被 finalizer 当作方法失败，而上游审计记为 0 method failures；在书面决定缺失时不能擅自放宽。
- CM 是已证实的 group hard-quota 中断：CopyKAT 写表明确返回 `Disk quota exceeded`。余量虽已改善至 101.588 GiB，仍低于 129.2–156.9 GiB 外推并低于 200 GiB 门槛，故保持 blocked。
- 下载 retry v3 已正常终止但带 1 个失败；GSE201347 的直接失败链是 NCBI FTP connection refused/远端断开、指数熔断和 transfer attempt limit。当前卡在尚未建立独立版本化 recovery，而非把 partial 误作完成。

## 解释与决定

- 三任务当前均非 active processing：免疫 blocked、CM blocked、下载 retry complete_with_failures 且单文件 recovery todo。
- CM group quota 改善是重要状态变化，但未越过既定执行门槛，本次不提交计算任务。
- 下载现有 quota 余量理论上大于缺失字节不等于授权或完整恢复条件已经满足；本任务仅监控和记录，不启动外部下载。
- `STATUS.md` 虽为当前状态入口，但因存在用户未提交改动且已明确隔离，本轮用 run/daily/roadmap/TODO 追加纠正，不触碰该文件。

## 失败尝试与限制

- 初次沙箱只读命令因系统不支持 unprivileged user namespace 而失败；改用获准的沙箱外只读命令完成审计。
- 专用 `apply_patch` 及沙箱外同名包装器均因同一 bwrap 基础设施故障无法读取目标仓库；经 dry-run 后使用系统 `patch` 应用完全相同的只追加 unified patch。
- `xfs_quota report -gh` 产生部分 `Operation not permitted`，但精确的 `quota -g -N USER002` 查询成功；本日志只采用后者的数值。
- 本次没有重新校验 GSE201347 的 gzip 流，因为文件大小未达到 expected size，已足以判定未通过；没有访问外部网络。

## 修改文件及原因

- 三份 run record：追加 07:00 最新证据、阻断原因和恢复边界。
- CM quota incident：追加精确余量改善但仍不足的状态。
- 2026-08-16 daily、ROADMAP、TODO、CHANGELOG：追加最新联合结论和执行门槛纠正。
- 本 Codex 日志：记录计划、命令类别、事实/解释、限制和提交范围。
- 未修改 `STATUS.md`、inventories 或任何分析运行目录；三个 run 状态分类未发生改变，因此不制造重复 inventory 行。

## 验证计划

- `git diff --check`。
- `bash scripts/validate_repository.sh`。
- 人工检查未暂存与已暂存 diff、文件大小、敏感绝对路径和暂存文件清单。
- 显式 `git add` 本日志列出的文件，确认 `STATUS.md` 不在 index 后普通 commit。

## 未验证与尚未解决

- 免疫 reference-insufficient skip 的验收决定仍待用户/分析负责人书面确认。
- CM `genome=hg20` 参数意图、额外至少 98.412 GiB group-quota 余量和 Stage 02–08 容量预算仍待解决。
- GSE201347 独立 recovery、完整性验证及 24-GSE core/sidecar guardrail 复核仍待执行。
- 本次用户只要求提交，未执行 push。

## 建议的下一步

1. 书面确认免疫 inferCNV skip 语义，再新建 Stage 13 attempt。
2. 将 CM group-quota 余量提高到至少 200 GiB，再从零基索引 580 新建 recovery。
3. 为 GSE201347 建立版本化可续传 recovery，并在 expected size 与 gzip 验证通过后才改为 verified。

## 建议提交信息

`docs: refresh three workstream status at 07:00`

## 验证结果

- 应用前 unified patch dry-run 通过；应用后所有目标文件均只增加行，零删除。
- `git diff --check` 通过。
- `bash scripts/validate_repository.sh`：0 errors、0 warnings。
- 本轮文件安全扫描未发现真实绝对路径、常见凭据模式或私钥标记。
- 提交前将再次检查已暂存文件列表和 cached diff，并确认 `STATUS.md` 仅保留在工作树。
