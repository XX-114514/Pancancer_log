# 项目概览

## 项目目标

PanCancer 项目旨在整理多癌种单细胞转录组数据，建立来源可追溯、矩阵语义明确、可分支验证的统一分析流程，并支持细胞注释、肿瘤细胞识别、CoVarNet 及跨队列验证。

本记录仓库当前以 `Project_v3` 为主。`Project_v1` 是较早期分析，保留重要数据和脚本；`Project_v2` 包含后续分析、免疫治疗数据、空间转录组相关内容以及 CoVarNet 代码；`Project_v3` 是当前主要分析数据集、正式运行目录和泛癌单细胞主流程记录。

## 主要分析问题

- 不同来源数据集的物种、样本单位、矩阵语义和基因空间是否适合进入同一发现队列？
- 如何在保留每数据集 QC 证据的同时，构建可扩展的统一整合对象？
- 如何为整合对象提供证据分层、可人工复核的细胞类型注释？
- 如何区分主发现队列、免疫治疗验证、snRNA 敏感性、空间分支和 processed-expression-only 分支？
- 如何使下游恶性识别与 CoVarNet 结果追溯到输入、代码、环境和参数？

## 输入

输入的权威实体保留在分析工作区，包括 GEO 数据、既有 AnnData、元数据、脚本和历史运行。本仓库只记录其稳定 ID 与逻辑路径：

- `${PROJECT_ROOT}/Project_v1`：早期数据与脚本。
- `${PROJECT_ROOT}/Project_v2`：后续分析、免疫治疗/空间分支与 CoVarNet 代码。
- `${PROJECT_ROOT}/Project_v3`：当前主要数据集、脚本、运行和报告。

具体数据集状态见 `inventories/datasets.tsv`。任何数字必须链接到可检查的 manifest 或 run record。

## 当前主要阶段

1. 数据来源与对象库存。
2. 数据恢复、矩阵语义和样本键审计。
3. 逐数据集 QC。
4. 统一 human discovery 整合。
5. 细胞类型注释。
6. 肿瘤细胞识别、CoVarNet 及分支分析。
7. 独立验证、报告和发布。

阶段状态和完成标准以 [ROADMAP.md](../ROADMAP.md) 为准。

## 预期结果

- 一个来源、输入、代码版本和环境均可追溯的泛癌单细胞分析主线；
- 可复核的主/亚型细胞注释；
- 明确区分 discovery、validation 和 sensitivity branches 的结果；
- 可由 run ID 重建的 CoVarNet 和其他下游分析；
- 记录失败、纠正和限制的项目历史。

## 已知边界

当前最新证据支持 Phase 07 整合完成，但不支持“最终细胞注释完成”或“最终生物学结论已获得”。`GSE162498` 仍是 extreme-scale pending 分支，`GSE211602` 已从 human-only 主分支排除。
