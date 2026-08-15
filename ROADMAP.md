# 项目路线图

状态值：`completed`、`in_progress`、`blocked`、`planned`、`待确认`。完成状态必须有可检查证据。

| 阶段 | 目标 | 输入 | 输出 | 依赖 | 完成标准 | 当前状态 |
| --- | --- | --- | --- | --- | --- | --- |
| P0 数据来源与库存 | 明确数据集、对象、矩阵语义和用途分支 | Project_v1/2/3 既有对象与官方元数据 | 数据集库存、对象库存、分支清单 | 可读取的对象元数据 | 数据集来源、物种、样本键和用途有可审计记录 | completed |
| P1 数据恢复与标准化 | 将可用来源转换为可追踪的 raw/counts 对象 | GEO/既有对象 | 每数据集 rawdata/checkpoint | 输入语义可确认 | 输入标识、矩阵语义和转换验证齐全 | in_progress |
| P2 逐数据集 QC | 对主发现队列进行一致且保留数据集差异的 QC | 冻结的 discovery manifest | QC 对象与验证表 | P0、P1 | 当前主线 33 个数据集全部通过验证 | completed |
| P3 统一整合 | 构建全量共同基因对象并完成多方法整合 | QC-pass human datasets | concat、Harmony、BBKNN、scVI、CellTypist、scIB | P2、基因空间审计 | 最终对象和三种整合均通过检查 | completed |
| P4 细胞注释 | 结合图聚类、marker、参考模型和污染审计生成主/亚型标签 | P3 BBKNN 对象 | 注释对象、证据表、复核清单 | P3 | 17 个 lineage 完整且交付审计通过 | completed |
| P5 肿瘤细胞与分支分析 | 进行恶性识别、CoVarNet 及必要敏感性分析 | P4 注释对象 | 恶性标签、CoVarNet 与分支结果 | P4 | active v2 的 721 样本、CoVarNet、LIANA、NicheNet 和最终审计完成 | in_progress |
| P6 跨队列验证 | 处理免疫治疗、snRNA、spatial、processed-only 分支 | 冻结的验证/敏感性清单 | 独立验证结果 | P0、P4/P5 | 分支边界明确，不与 discovery 泄漏，最终交付审计通过 | in_progress |
| P7 报告与发布 | 汇总方法、结果、局限和可复现证据 | 已验证结果 | 报告、图表、共享记录 | P4-P6 | 结论可追溯到 run、输入、参数和验证 | in_progress |

## 当前里程碑

- P1：24-GSE `Download_now` initial run 保留 4 verified/118 failed 的失败证据；141-task retry v3 已结束，127 verified、7 skipped_verified、6 unavailable_upstream、1 failed。GSE201347 partial 仍需独立 recovery attempt。
- P3：31 个 GSE、1,963,745 个细胞的全量对象完成 Harmony、BBKNN、scVI、CellTypist 和 scIB 评价。
- P4：17 major lineages、136 lineage clusters、106 final annotations，最终审计 `PASS`。
- P5：主 CoVarNet 发现分支完成固定 `K=9`；恶性/通讯 v2 Stage 01 严格完成 580/721。recovery `20260812_124600_cm_stage01_from_0462` 在第 581 个样本因组配额耗尽终止；当前约 2.50 GiB 余量远低于剩余 Stage01 约 129–157 GiB 的估算，至少恢复 200 GiB 组配额余量后才从零基索引 580 新建 attempt；Stage 02–08 未开始。
- P6：旧三个固定模块免疫治疗投影完成；新的四队列整合已完成 19/19 cell2location 与 71,398-spot 严格 merge，但截至 2026-08-12 无新 Stage 13 attempt，最终验收仍因 11 个 reference-insufficient inferCNV skip 的状态策略冲突而阻断。
- P7：三个 extension 均已建立正式运行记录；下载 retry v3 已生成 final summary，免疫与 CM 的最终 audit 仍未齐全，远程同步仍受 GitHub 认证阻塞。

## active extension 完成标准

### 四队列免疫治疗整合

1. 已完成：19/19 张 GSE273952 切片各自通过 cell2location attempt 验证。
2. 已完成：71,398 个 spots 精确合并，无 missing/extra、无非有限/负值/zero-sum abundance。
3. 待完成：统一 11 个 inferCNV reference-insufficient skip 的验收语义，`FINAL_DELIVERABLE_AUDIT.json` 通过，并保留 7 位 unresolved response 的不确定性。

### 恶性识别与 CM 通讯

1. 721 个 active samples 全部得到 non-failed terminal Stage 01 结果。
2. GSE 排除、GSE166555 重组、global positions 和恶性合并零冲突。
3. active v2 的 CoVarNet K9、sample × CM LIANA、NicheNet 非空验证和 `FINAL_AUDIT.json` 通过。
4. 发布前确认 CopyKAT 实际 `genome=hg20` 是否符合预期。

### 新数据下载

1. 已完成：141-task retry 终止并生成 machine-readable final summary；终态为 127 verified、7 skipped_verified、6 unavailable_upstream、1 failed。
2. 每个成功文件通过与类型相符的 gzip/tar/zip/size/MD5 检查。
3. `unavailable_upstream` 与网络失败不伪装为 verified。
4. core GEX、sidecar、human/xenograft、baseline/repeated-measure guardrail 经人工复核后才进入预处理。
5. GSE201347 的 14,413,317,040/19,071,156,087-byte partial 必须在新的版本化 recovery attempt 中恢复并重新验证；不得改写 retry v3 的失败历史。

## 关键依赖

- `/data4` 文件系统整体仍有约 771 GiB available，但 `USER002` 组受 1,572,864,000-block hard limit 约束；2026-08-16 仅余约 2.50 GiB。CM Stage01 新 recovery 至少需要 200 GiB 组配额余量，且下游 Stage 02–08 另需容量预算。
- `GSE278694` 的移除原因必须形成可审计决定，否则最终发现队列范围不能算完全冻结。
- `GSE162498` 的 extreme-scale 处理决定影响主队列是否扩展。
- low-confidence、unresolved 和 ambiguous annotations 需要人工复核后才能作为发布版标签。
- 外部验证必须继续使用冻结 mapper 和固定 CoVarNet 权重，避免在验证队列重新拟合造成信息泄漏。
- 每个阶段的权威证据通过 `inventories/runs.tsv` 和正式 run record 链接，不在路线图复制完整日志。
