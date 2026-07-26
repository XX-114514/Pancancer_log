# 术语表

| 术语 | 类型 | 当前定义 | 证据/备注 |
| --- | --- | --- | --- |
| PanCancer | 项目 | 当前泛癌单细胞转录组分析总项目 | 项目根 README 与用户确认 |
| Project_v1 | 目录/历史阶段 | 较初始的分析，保留重要数据与脚本 | 用户确认 |
| Project_v2 | 目录/历史阶段 | 后续分析，含免疫治疗数据、空间转录组相关内容、其他分析和 CoVarNet 脚本 | 用户确认 |
| Project_v3 | 目录/当前阶段 | 当前主要分析数据集和泛癌单细胞主流程记录 | 用户确认 |
| CoVarNet | 方法/代码 | 基于样本与细胞亚群组成进行模块/网络分析的工作流 | 精确定义和版本仍应由对应方法与代码 commit 固定 |
| discovery branch | 分支 | 用于主模型发现的队列，必须与验证分支隔离 | 当前主 run |
| validation branch | 分支 | 用于独立验证、不得泄漏进 discovery 的数据 | 具体成员待 inventory 完整记录 |
| sensitivity branch | 分支 | 用于 snRNA、空间或其他语义差异的敏感性分析 | 具体规则待确认 |
| processed-expression-only | 分支 | 没有可靠 raw counts、仅有处理后表达矩阵的数据 | 不应默认进入依赖 counts 的流程 |
| archival shards | 产物类型 | Phase 07 对完整整合矩阵按行分片保存的可审计对象 | 最新主 run 记录 39 shards |
| modeling branch | 产物类型 | 从 archival plan 按数据集/全局上限抽取的建模对象 | 最新主 run 记录 249,990 cells |
| sample unit | 统计单位 | runner 识别的样本键对应单位 | 不同数据集列名不同，必须逐数据集审计 |
| matrix semantics | 数据性质 | 矩阵是否为 raw counts、log-like、空间或其他表达语义 | 进入统一流程前必须记录 |
| GSE162498 | 数据集 | 当前 extreme-scale pending 数据集 | 尚未观察到进入 32-dataset human integration |
| GSE211602 | 数据集 | 在 Phase 07 被识别为非人并排除出 human-only 整合 | 后续用途待确认 |
| run ID | 标识符 | 一次可独立复现和审计的正式运行标识 | 见命名规范 |
| ADR | 文档类型 | Architecture Decision Record，记录重要技术决策及理由 | `decisions/` |

术语定义变化时，应在新的日志或 ADR 中记录依据。待确认项不能擅自补全。
