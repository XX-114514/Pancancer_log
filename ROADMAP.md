# 项目路线图

状态值：`completed`、`in_progress`、`blocked`、`planned`、`待确认`。完成状态必须有可检查证据。

| 阶段 | 目标 | 输入 | 输出 | 依赖 | 完成标准 | 当前状态 |
| --- | --- | --- | --- | --- | --- | --- |
| P0 数据来源与库存 | 明确数据集、对象、矩阵语义和用途分支 | Project_v1/2/3 既有对象与官方元数据 | 数据集库存、对象库存、分支清单 | 可读取的对象元数据 | 数据集来源、物种、样本键和用途有可审计记录 | completed |
| P1 数据恢复与标准化 | 将可用来源转换为可追踪的 raw/counts 对象 | GEO/既有对象 | 每数据集 rawdata/checkpoint | 输入语义可确认 | 输入标识、矩阵语义和转换验证齐全 | in_progress |
| P2 逐数据集 QC | 对主发现队列进行一致且保留数据集差异的 QC | 冻结的 discovery manifest | QC 对象与验证表 | P0、P1 | 当前主线 33 个数据集全部通过验证 | completed |
| P3 统一整合 | 构建全量共同基因对象并完成多方法整合 | QC-pass human datasets | concat、Harmony、BBKNN、scVI、CellTypist、scIB | P2、基因空间审计 | 最终对象和三种整合均通过检查 | completed |
| P4 细胞注释 | 结合图聚类、marker、参考模型和污染审计生成主/亚型标签 | P3 BBKNN 对象 | 注释对象、证据表、复核清单 | P3 | 17 个 lineage 完整且交付审计通过 | completed |
| P5 肿瘤细胞与分支分析 | 进行恶性识别、CoVarNet 及必要敏感性分析 | P4 注释对象 | 恶性标签、CoVarNet 与分支结果 | P4 | 各分支输入、运行和结果验证完成 | in_progress |
| P6 跨队列验证 | 处理免疫治疗、snRNA、spatial、processed-only 分支 | 冻结的验证/敏感性清单 | 独立验证结果 | P0、P4/P5 | 分支边界明确，不与 discovery 泄漏 | in_progress |
| P7 报告与发布 | 汇总方法、结果、局限和可复现证据 | 已验证结果 | 报告、图表、共享记录 | P4-P6 | 结论可追溯到 run、输入、参数和验证 | in_progress |

## 当前里程碑

- P3：31 个 GSE、1,963,745 个细胞的全量对象完成 Harmony、BBKNN、scVI、CellTypist 和 scIB 评价。
- P4：17 major lineages、136 lineage clusters、106 final annotations，最终审计 `PASS`。
- P5：CoVarNet 发现分支完成固定 `K=9` 的 9 个模块；恶性/CNV 分支的完成状态仍需单独核验。
- P6：三个免疫治疗队列的固定模块投影完成，但缺少可支持确定性疗效/生存结论的独立临床 endpoint。
- P7：本地正式运行记录已补齐；远程同步受 GitHub 认证阻塞。

## 关键依赖

- `GSE278694` 的移除原因必须形成可审计决定，否则最终发现队列范围不能算完全冻结。
- `GSE162498` 的 extreme-scale 处理决定影响主队列是否扩展。
- low-confidence、unresolved 和 ambiguous annotations 需要人工复核后才能作为发布版标签。
- 外部验证必须继续使用冻结 mapper 和固定 CoVarNet 权重，避免在验证队列重新拟合造成信息泄漏。
- 每个阶段的权威证据通过 `inventories/runs.tsv` 和正式 run record 链接，不在路线图复制完整日志。
