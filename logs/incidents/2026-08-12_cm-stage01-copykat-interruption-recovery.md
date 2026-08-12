# Incident: CM Stage01 第 462 个样本 CopyKAT 非终态中断与恢复

## Summary

active v2 Stage 01 在第 462/721 个生物学样本停止。inferCNV 已成功，但 CopyKAT 在 Step 7 后失去进程，只留下 raw CNA 中间矩阵，没有 prediction、summary 或 terminal marker。2026-08-12 完成非破坏修补并以新 attempt 重新计算。

## Impact

- 前 461 个顺序样本保持已完成，不重算、不覆盖。
- 第 462 个样本在恢复前不得计入完成数。
- active v2 Stage 02–08 未开始，`FINAL_AUDIT.json` 仍待生成。
- 旧 attempt 的 5 个诊断文件完整保留，其中 raw gene-by-cell 矩阵约 1.106 GB。

## Detection

2026-08-12 审计发现最后 heartbeat 停在 2026-08-10 08:02，主 runner、Python worker 和 CopyKAT R 子进程均不存在；样本目录无 prediction、summary、`COMPLETED.json` 或 `FAILED.json`。

## Timeline

| Time (+08:00) | Event | Evidence |
| --- | --- | --- |
| 2026-08-10 05:51 | 第 462 个样本 inferCNV 完成并启动 CopyKAT | `${CM_RUN_ROOT}/logs/stages/01_sample_cnv.log` |
| 2026-08-10 08:02 | 最后一条旧 heartbeat；Step 7、RSS 约 18.4 GiB、CPU 约 199.5% | 同上 |
| 2026-08-12 10:06 | 确认无活动进程、无 terminal marker | run record 审计快照 |
| 2026-08-12 12:42 | 旧 `attempt_01` 封存为 interrupted；创建 `attempt_02` 并重跑 | `INTERRUPTED.json`、`STARTED.json` |
| 2026-08-12 12:44 | 新 CopyKAT 首条 heartbeat，Step 3、RSS 6,879,140 kB、无 swap | recovery log |

## Evidence

- Recovery log：`${CM_RUN_ROOT}/logs/recovery/20260812_124600_cm_stage01_from_0462.log`
- Recovery history：`${CM_RUN_ROOT}/status/recovery/recovery_history.tsv`
- 旧现场封存：`${CM_RUN_ROOT}/results/per_sample/GSE299340__GSE299340_GSE299340__GSM9037562_297_50a832cca71a/attempt_01/INTERRUPTED.json`
- 新 attempt provenance：`${CM_RUN_ROOT}/results/per_sample/GSE299340__GSE299340_GSE299340__GSM9037562_297_50a832cca71a/attempt_02/STARTED.json`
- Runner：`${CM_RUN_ROOT}/scripts/01_run_sample_cnv.py`
- Recovery launcher：`${CM_RUN_ROOT}/run_stage01_recovery_20260812.sh`

## Root cause

未确认。旧日志没有 Python/R traceback，也没有进程退出码；当前账户无权读取对应内核 journal。主机后来发生重启，但时间证据不足以把重启或 OOM 写成已证实根因。

## Contributing factors

- CopyKAT 1.1.0 baseline/GMM 分支计算时间和内存随样本变化显著。
- 原 attempt 生成时还没有 `STARTED.json`，外部硬中断后无法自动形成终态。
- raw CNA 中间矩阵不是 CopyKAT 支持的可靠断点，不能冒充 prediction 复用。

## Resolution

1. completion skip 改为验证样本 UID/token、细胞数、三种方法状态、`cell_evidence.tsv.gz` 路径和最小列合同。
2. 恢复前检查非终态 attempt；确认无活动 runner 后追加 `INTERRUPTED.json`，不删除或覆盖旧文件。
3. 新 attempt 写 `STARTED.json`，记录 PID、样本、分母、manifest、CopyKAT cores 和 heartbeat。
4. 新增带 `flock`、独立 recovery log/history、默认 `start-index=461` 的 Stage01 恢复启动器。
5. 因普通 `setsid` 会被当前执行器清理，最终使用 tmux 会话 `cm_stage01_20260812_1246` 托管。

## Methods and parameters

- 样本：GSE299340 / GSM9037562_297；8,734 cells，4,123 candidates，3,647 references。
- inferCNV：raw counts、window size 250、reference q99；本样本 threshold 0.10957874104380605。
- CopyKAT：1.1.0；`genome=hg20`、`id.type=S`、`ngene.chr=5`、`win.size=25`、`KS.cut=0.1`、Euclidean、2 cores。
- scMalignantFinder：仅在 CopyKAT 后运行，library-size 10,000、no log1p；不改变 inferCNV ∩ CopyKAT 严格交集。
- Heartbeat：60 秒，记录 PID、elapsed、CPU、RSS、swap、threads、stdout/stderr bytes 和 stdout 尾行。

## Validation

- 目标 scanpy Python `py_compile`：通过。
- recovery launcher `bash -n`：通过。
- CopyKAT/Matrix R preflight：通过，CopyKAT 1.1.0。
- 真实已完成样本 completion 合同：通过；零基索引 461 精确对应第 462 个样本。
- tmux、flock、runner、R 子进程和首条 CopyKAT heartbeat：通过。
- 恢复启动时主机约 949 GiB available memory，`/data4` 约 867 GiB available。

## Prevention

- 只跳过通过 artifact 合同的 `COMPLETED.json`。
- 每个新 attempt 必须有 start + terminal provenance。
- 长 R 子进程必须保留 60 秒 heartbeat。
- 任何可复用 prediction 必须满足 cell ID 行数相等、零 missing/extra/duplicates，并记录 hash。
- Stage 01 全 terminal 前不得发布 active v2 的 CoVarNet、LIANA、NicheNet 或 final audit。

## Follow-up tasks

- `CM-001`：继续监控剩余 260 个样本，Stage 01 全 terminal 后重算 Stage 02–08。
- `CM-002`：发布前书面确认实际 `genome=hg20` 参数意图。
