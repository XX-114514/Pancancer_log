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

## 关键和探索节点即时维护

2026-10-11用户再次授权并要求在关键/探索节点持续维护本仓库和普通push。这是进行任务时的节点工作流，没有新建定时任务或宣称GitHub实时镜像服务器。

新计划、阶段变化、方法/结论改变、失败与恢复、审计/交接完成时，尽快更新同任务Codex日志、daily、正式run及状态/索引；长任务的pending检查点也可同步。每个节点记录事件UTC/UTC+8和取证时刻、范围/版本、方法、结果、失败修正、证据SHA、结论边界与下一步。后补日志明确为补录，不伪造早期提交时点。

发布前检查已有修改、脱敏内容、仓库验证与diff；显式stage、普通commit/push到用户指定origin main，随后ls-remote或独立fetch核对完整SHA。保存内容commit/发布时间/核验方法。失败时记录待同步，继续排查，直到远端核验完成才声明云端成功。不得为了push覆盖他人修改或重写历史。

最新科学节点：[2026-10-11 18:20:36 UTC+8证据](../inventories/analysis/2026-10-11/README.md)。下载库存仍使用前述独立取证时间。

## 每次任务同步与最新下载依据（2026-10-11）

用户重申每一次任务都需及时普通push；以§关键和探索节点流程执行，不只完成全部项目时同步。最新下载证据为[18:56—19:02报告](../inventories/runtime/2026-10-11/DOWNLOAD_READINESS.md)；此下载报告优先于旧REPORT中的下载线，科学家族/空间线仍以其专属较晚登记为准。v6监控刷新不能证明下载运行，queue_terminal不能证明全部成功。

最新科学空间节点更新为[19:27:55预检](../runs/20261011_spatial_contact_native_preflight_v3.md)，不改变下载线取证。原生端点/正确作者映射v14已接受，空间邻近及机制仍待；v2作者方向错误已显式替代。事件与发布收据分开，关键节点的普通push仍逐次独立核对。

最新科学节点[19:54条件空间验收](../runs/20261011_spatial_conditional_neighborhood_v1.md)：v15/0冲突/本阶段进程全部终止、主45BH无显著；下载线不随此刷新。关键节点普通push独立SHA逐次保存。

最新科学节点[20:29接收程序验收](../runs/20261011_spatial_receiver_programme_exposure_v1.md)：v16/0冲突，70主模型/22可估计/0名义及主BH显著；生产和审计均終止。此节点不刷新下载库存，当前通用程序阴性与计算通过分开记录。

## 技术QC与正式来源验收补充

21:18 UTC+8四分片已[全量技术验收](../inventories/runtime/2026-10-11/QC_ACCEPTANCE.md)，technical_QC_acceptance_20261011_v1登记，不以accepted_*旧命名覆盖source门槛。patientfix_v2按稳定ID附加头颈癌元数据，原QC不改；collector旧pending仍指正式来源，技术状态以本新登记为准。clinical_sidecar_sha256=null不等于无来源证据，须核preflight.clinical_source、官方metadata与真实registry/obs。正式source acceptance正数登记出现时再扩展export合同。
