# Incident: CM Stage01 第 581 个样本因组配额耗尽终止

## Summary

recovery `20260812_124600_cm_stage01_from_0462` 将 active v2 Stage01 从 461 个严格完成样本推进到 580 个。2026-08-14 03:40 +08:00，第 581/721 个样本在 CopyKAT 写表时收到 `Disk quota exceeded`，runner 随后退出。

## Impact

- 前 580 个样本的有效 `COMPLETED.json` 和 `cell_evidence.tsv.gz` 保持不变。
- 第 581 个样本没有 canonical prediction/summary、cell evidence 或 terminal marker，不能计入完成。
- Stage 02–08 未开始，`FINAL_AUDIT.json` 不存在。
- 当前还剩 141 个样本、501,318 个细胞。

## Evidence

- Recovery log：`${CM_RUN_ROOT}/logs/recovery/20260812_124600_cm_stage01_from_0462.log`
- 非终态 attempt：`${CM_RUN_ROOT}/results/per_sample/GSE127465__GSE127465_human_p5t1_26fd90449a28/attempt_01/`
- 原生 CopyKAT prediction、CNA results、clustering RDS 和约 22.1 MB raw CNA matrix 均保留。
- 2026-08-16 严格合同审计：580/721 连续有效，第一个未通过样本为 581。

## Root cause

已确认是 `/data4` 的 `USER002` 组 hard quota，而不是 inode 耗尽或 CopyKAT 算法停滞：

- 文件系统：XFS，启用 user/group quota。
- 故障时组用量达到 1,572,864,000-block hard limit。
- 文件系统整体仍有约 769 GiB available；inode 使用约 4%。
- CopyKAT stderr 明确返回 `Error writing to connection: Disk quota exceeded`。

## Capacity reassessment

用户释放空间后，2026-08-16 组用量为 1,570,240,644/1,572,864,000 blocks，仅余 2,623,356 blocks，约 2.50 GiB。

对样本 462–580 的 119 个有效 recovery attempts 按实际分配磁盘块审计：

- 357,794 cells 共占 92.189 GiB。
- 单样本占用中位数 0.625 GiB；P75 1.102 GiB；P90 1.691 GiB；P95 2.491 GiB；最大 3.383 GiB。
- 剩余 501,318 cells 按加权均值估算 129.2 GiB，按 P75 bytes/cell 估算 156.9 GiB。

因此 2.50 GiB 不满足安全重启门槛；Stage01 建议至少保留 200 GiB 组配额余量，下游 Stage 02–08 另行预算。

## Resolution status

- 2026-08-16 没有提交新 recovery，避免可预见的重复 EDQUOT。
- 没有删除、覆盖或改写任何历史 attempt。
- 配额门槛满足后，从零基索引 580 创建新的 recovery ID 和 append-only attempt。
- 当前原生 CopyKAT prediction 只在 cell IDs 行数相等、零 missing、零 extra、零 duplicate 并记录来源与 hash 后允许复用。

## Prevention

1. 每次恢复前同时核对 filesystem free space 和实际 user/group quota headroom。
2. 根据已完成 attempt 的分配块占用和剩余细胞数估算容量，不以 `df` 可用空间替代 quota 检查。
3. 推荐门槛未满足时不提交长任务；记录为资源阻断。
4. 继续使用 append-only attempt、严格 completion artifact 合同和 60 秒 heartbeat。
5. Stage01 全 terminal 前不得发布 active-v2 CoVarNet、LIANA、NicheNet 或 final audit。
