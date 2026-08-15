# Codex 任务日志：复核 CM recovery 容量、重提判定与记录提交

## 用户要求

用户说明已释放空间，要求重新确定并提交任务；延续前序要求更新 project-records 并进行 Git 交付。

## 开始前状态

- CM recovery 已于第 581/721 个样本因 `Disk quota exceeded` 退出。
- 严格完成合同通过前 580 个样本。
- project-records `main` 领先 `origin/main` 两个提交。
- `STATUS.md` 存在用户所有的未提交首行空白改动，本任务保持不修改、不暂存。

## 计划

1. 复核 group quota、进程和严格恢复断点。
2. 用已完成 attempts 的真实分配块估算剩余 Stage01 容量。
3. 只有资源门槛通过才建立新 recovery ID 并启动。
4. 更新 run、incident、daily、roadmap、TODO 和 inventories。
5. 验证并仅提交本任务文件；private origin 未明确确认前不绕过 push 安全门槛。

## 检查的文件

- `${CM_RUN_ROOT}/manifests/samples_analysis_v2.tsv`
- `${CM_RUN_ROOT}/logs/recovery/20260812_124600_cm_stage01_from_0462.log`
- `${CM_RUN_ROOT}/scripts/01_run_sample_cnv.py`
- 第 581 个样本的非终态 `attempt_01`
- project-records 的 AGENTS、README、STATUS、ROADMAP、TODO、数据政策、维护工作流、CM run 与 incident。

## 执行的检查

- `df -h` 与 group quota 审计。
- `pgrep` 检查 runner、Python 和 CopyKAT。
- 用 runner 自身 `validate_completed_attempt()` 逐样本验证。
- 用文件 `st_blocks` 统计样本 462–580 最新有效 attempts 的实际分配磁盘空间，并按剩余细胞数外推。
- Git status、diff、仓库验证和显式暂存边界检查。

## 关键发现

### 已观察事实

- 580/721 个样本连续通过严格完成合同；第 581 个样本非终态，剩余 141 个样本、501,318 cells。
- 当前没有重复 runner 或 CopyKAT 进程。
- 组配额仅余 2,623,356 blocks，约 2.50 GiB。
- 119 个 recovery attempts 的实际占用为 92.189 GiB；剩余需求估算约 129–157 GiB。

### 解释与决定

- 文件系统整体 free space 不能代表 `USER002` 组可写余量。
- 2.50 GiB 不足以安全完成 Stage01；本日不提交计算任务，避免确定性重复 EDQUOT。
- 安全重启门槛设为至少 200 GiB 组配额余量。

## 修改

- 更新 CM run record、ROADMAP、TODO、CHANGELOG、runs/artifacts inventories。
- 新增本日志、2026-08-16 daily log 和组配额 incident。
- 未修改 `STATUS.md`，以保留其既有用户改动；当前 CM 权威状态由 run record 和 incident 纠正。

## 验证

- `bash scripts/validate_repository.sh`：0 errors、0 warnings；`git diff --check` 通过；三个 inventory TSV 列结构一致；无明显 secret、敏感绝对路径、超 5 MiB 文件、禁止数据格式或断裂本地链接。暂存后将再次核对 payload，且明确排除 `STATUS.md`。

## 未验证与阻断

- 未提交新 CM recovery，因为容量门槛未通过。
- Stage 02–08 和 final audit 均未开始。
- 普通 push 仍要求用户明确确认已配置 origin 是可信 private 目的地。

## 建议提交信息

`docs: record CM group-quota recovery blocker`
