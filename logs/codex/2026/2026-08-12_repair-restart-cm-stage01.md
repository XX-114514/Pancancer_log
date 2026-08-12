# Codex 任务日志：修复并重启 CM Stage01

## 任务

检查 CM 通讯报错，修复恢复逻辑，并从第 462 个样本重新运行 active v2 Stage 01。

## 结论

这不是 LIANA 通讯统计代码报错，而是其上游样本级 CopyKAT 非终态中断。旧 attempt 没有 prediction、summary 或终态 marker，不能从 raw CNA 中间矩阵可靠断点恢复。

## 执行

1. 审计 active manifest、尾部日志、进程、attempt 文件和系统事件可见性。
2. 保留前 461 个完整样本和第 462 个旧 partial。
3. 增加 completion artifact 合同、孤儿 attempt 封存和 STARTED provenance。
4. 新增独立 Stage01 recovery launcher，固定零基 `start-index=461`、CopyKAT 2 cores、60 秒 heartbeat。
5. 目标 Python/R/bash 预检通过后，用 tmux 会话启动。
6. 第一条新 heartbeat 已确认 CopyKAT 在 Step 3 持续计算。

## 当前状态

- Recovery ID：`20260812_124600_cm_stage01_from_0462`
- tmux：`cm_stage01_20260812_1246`
- 样本：462/721，GSE299340 / GSM9037562_297
- 新 attempt：`attempt_02`
- inferCNV：success，threshold 0.10957874104380605
- CopyKAT：active，参数与主队列保持一致
- Stage 02–08：pending

## 记录

详见 [incident](../../incidents/2026-08-12_cm-stage01-copykat-interruption-recovery.md) 和 [CM run record](../../../runs/20260731_malignancy_communication.md)。
