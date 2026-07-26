# 项目路线图

状态值：`completed`、`in_progress`、`blocked`、`planned`、`待确认`。完成状态必须有可检查证据。

| 阶段 | 目标 | 输入 | 输出 | 依赖 | 完成标准 | 当前状态 |
| --- | --- | --- | --- | --- | --- | --- |
| P0 数据来源与库存 | 明确数据集、对象、矩阵语义和用途分支 | Project_v1/2/3 既有对象与官方元数据 | 数据集库存、对象库存、分支清单 | 可读取的对象元数据 | 数据集来源、物种、样本键和用途有可审计记录 | completed |
| P1 数据恢复与标准化 | 将可用来源转换为可追踪的 raw/counts 对象 | GEO/既有对象 | 每数据集 rawdata/checkpoint | 输入语义可确认 | 输入标识、矩阵语义和转换验证齐全 | in_progress |
| P2 逐数据集 QC | 对主发现队列进行一致且保留数据集差异的 QC | 冻结的 discovery manifest | QC 对象与验证表 | P0、P1 | 当前 human 主线数据集全部通过验证 | completed |
| P3 统一整合 | 建立可扩展 archival shards 与 capped modeling branch | QC-pass human datasets | archival shards、HVG/PCA/Harmony、UMAP/Leiden | P2、基因空间审计 | 所有计划 shard 可读，模型分支通过验证 | completed |
| P4 细胞注释 | 结合 marker、参考模型和数据集证据生成主/亚型标签 | P3 模型对象 | 注释对象、证据表、复核清单 | P3 | 标签规则、证据、置信度和人工复核均可追踪 | in_progress |
| P5 肿瘤细胞与分支分析 | 进行恶性识别、CoVarNet 及必要敏感性分析 | P4 注释对象 | 恶性标签、CoVarNet 与分支结果 | P4 | 输入样本数足够，运行和结果验证完成 | planned |
| P6 跨队列验证 | 处理免疫治疗、snRNA、spatial、processed-only 分支 | 冻结的验证/敏感性清单 | 独立验证结果 | P0、P4/P5 | 分支边界明确，不与 discovery 泄漏 | planned |
| P7 报告与发布 | 汇总方法、结果、局限和可复现证据 | 已验证结果 | 报告、图表、共享记录 | P4-P6 | 结论可追溯到 run、commit、输入和验证 | planned |

## 关键依赖

- `GSE162498` 的 extreme-scale 处理决定影响主队列是否扩展。
- P4 注释完成前，不应把 P5/P6 的结果称为最终泛癌结论。
- 每个阶段的权威运行证据应通过 `inventories/runs.tsv` 链接，而不是在本文件复制完整历史。
