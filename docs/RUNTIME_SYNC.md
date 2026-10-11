# 服务器状态与GitHub同步

## 权威和新鲜度

最新已提交快照：[2026-10-11报告](../inventories/runtime/2026-10-11/REPORT.md)。时间是取证时间，不保证打开GitHub时服务器状态相同。服务器可读下载日志继续更新；本次没有创建自动GitHub push或定时任务。

证据顺序：最终交接登记 + 独立审计 → 当前状态汇总 → producer初始manifest。早期pending状态应保留，不能覆盖后续成功登记。进程角色分为pipeline、monitor、infrastructure；代理不能算科学任务。

## 安全取证命令

在服务器检出明确代码commit，使用Python ≥3.10，并通过外部环境变量提供真实根目录：

```bash
python scripts/collect_runtime_snapshot.py \
  --project-root "${PANCANCER_ROOT}" \
  --output "${SNAPSHOT_OUTPUT}"
python scripts/validate_runtime_snapshot.py "${SNAPSHOT_OUTPUT}"
```

输出文件必须是新的日期/版本名，采集器拒绝覆盖旧快照。采集器只读取当前用户可读的已知路径和进程、输出白名单汇总，**不重新运行科学分析**。外部逐患者报告只留下文件SHA、大小和逻辑路径。

将审查后的JSON复制到 `inventories/runtime/YYYY-MM-DD/`，在本仓库生成新的可读报告：

```bash
python scripts/render_runtime_report.py "${SNAPSHOT_JSON}" --output "${REPORT_OUTPUT}"
bash scripts/cloud_check.sh
git diff --check
```

随后更新STATUS、任务板、run索引和Codex日志，显式stage文件，普通commit/push。检查远端main SHA与本地提交一致。新schema/新增队列/正式acceptance登记出现时，先更新采集合同及测试，不能让旧白名单默默漏算。

## 状态定义

- download_verified：有下载证明/完整源库存，不是临床结局验收。
- source_cells：源文件数量，可能包含基线重复、非primary或不同技术。
- QC分片：当前写入后读回结果，仍待来源级最终验收。
- final_accepted_new_cells：只有来源、去重、QC和模态门槛均通过的正式acceptance登记可增加。
- Tier A/B/C：结局证据层级，不能凭标题或治疗关键词推断疗效/生存。
- Metadata-only：CRC细胞metadata数量不进入已下载计数矩阵总量。

当前exporter适配固定7个CXG、3个GEO和RMC源库存，第二批按文件校验记录；新增数据需显式扩展合同。当前final acceptance为未登记，不是永久固定科学结论。
