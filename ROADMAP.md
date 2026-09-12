# 项目路线图

状态值：`completed`、`in_progress`、`blocked`、`planned`、`待确认`。完成状态必须有可检查证据。
## 2026-09-12 V8 开发候选索引

- V8 rerun2 已形成全量元数据软件候选及外部 hash/状态索引；其完整性、schema、round-trip 和
  128 个 V7 字段合同均为 `PASS`，但 release status 为 `NOT_FROZEN`。
- 最新 root audit 记录 `PASS_DEVELOPMENT_SOFTWARE_CANDIDATE_NOT_FROZEN`、V7 default
  与外部记录的 `pytest_passed: 167`；这些是软件候选门控证据，而非 biological freeze 或
  独立准确率结论。
- V7 继续是默认 identity/CNV 参考。该候选不代表独立准确率、classifier 执行、真值验证或
  生物学完成，且 policy 明确禁止自动参考晋升和正式 downstream mutation。
- 独立 scANVI reference 仍为 `BLOCKED_NO_LOCAL_CELLXGENE_CENSUS_CLIENT`。下一步仅在明确
  授权安装或提供可审计 client 后，先进行固定 Census 版本的前台 metadata-only 查询并审查
  source/donor/GEO、ontology 与 license；解除该门控前不得训练、mapping 或主张独立性能。
- 下一里程碑是在新的版本化验证 run 中完成独立生物学证据和下游影响审计；不能以本软件候选
  覆盖 V7 或已有 downstream 历史。

## 2026-09-06 路线更新

- P4/P5 身份与 CNV：V7 已冻结并成为当前默认参考，进入 append-only 维护。
- P7 发布治理：V7 与通讯证据链已形成轻量 Git releases。
- 下一里程碑：在 V7 身份/CNV 上版本化重跑依赖身份的 CoVarNet/LIANA 分支，并与
  当前 V5-derived downstream release 做差异审计。
- 下一里程碑：建立 Prism panel-to-source 索引、补齐软件/环境摘要和期刊要求的
  supplementary README。
- 数据发布继续采用 Git 元数据/小表/图 + 外部重型对象哈希的双层结构。

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
- P4/P7：V8 full-metadata rerun2/root audit 为
  `PASS_DEVELOPMENT_SOFTWARE_CANDIDATE_NOT_FROZEN`，外部记录 167 项 pytest 通过；它保留
  V7 作为默认参考，独立 reference 仍 `BLOCKED_NO_LOCAL_CELLXGENE_CENSUS_CLIENT`，不能计作
  P4 生物学标签升级或 P7 正式 release 完成。
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

## 2026-08-16 07:00 依赖纠正快照

- 上述 03:37 记录的约 2.50 GiB 是历史快照；最新精确 `USER002` group-quota 余量为约 101.588 GiB。该值仍低于剩余 CM Stage01 的 129–157 GiB 外推和 200 GiB 安全门槛，因此 P5/CM 状态仍为 `blocked`，恢复入口仍为零基索引 580 的新 append-only attempt。
- P6/免疫治疗仍由 11 个 reference-insufficient inferCNV skip 的验收语义和缺失 final audit 阻断；P1/下载 retry v3 已终止，但 GSE201347 独立 recovery 与下载后 guardrail 复核仍待执行。
