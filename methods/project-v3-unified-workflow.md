# Project_v3 统一泛癌单细胞工作流

- 状态：当前方法摘要
- 最近证据日期：2026-07-26
- 权威运行：`20260725_120105_unified_pancancer_covarnet`

## 目的

从来源和矩阵语义不完全一致的多数据集对象出发，建立分支明确、逐数据集可审计、能够扩展到大型数据的 human discovery 统一分析对象。

## 当前认可流程

1. **库存与分支冻结**
   审计对象来源、物种、矩阵语义、样本键和用途。明确 discovery、validation、sensitivity、spatial 和 processed-only 分支，避免数据泄漏。
2. **逐数据集 QC**
   在数据集内部使用已确认的样本键进行 QC 和 doublet 处理，保留输入/输出细胞数、样本数、counts 来源与失败记录。
3. **基因空间审计**
   标准化 gene symbol，检查重复基因和跨数据集交集。非人数据不进入 human-only 主整合。
4. **归档分片**
   使用显式 row-shard plan 写入 archival shards，避免意外构建稠密全矩阵。
5. **建模分支**
   按冻结的 dataset/global caps 从归档分片构建 capped modeling branch。
6. **建模预处理与整合**
   执行 HVG、PCA、Harmony、neighbors、UMAP 和 Leiden；参数必须在正式 run 中固定。
7. **证据化注释**
   当前待执行。应结合 marker evidence、参考模型、数据集来源和人工复核，分别记录主类群、亚型、置信度和争议项。
8. **下游分支**
   注释通过后再进入恶性识别、CoVarNet 和独立验证；不得把中间聚类直接写成最终细胞类型。

## 已观察验证

- Phase 06：33/33 datasets 验证成功，2,062,230 QC-pass cells。
- Human-only Phase 07：32 datasets，2,040,048 cells，11,063 common standardized genes。
- Archival：39/39 shards 可读，共 2,040,048 cells。
- Modeling：249,990 cells；HVG/PCA/Harmony 输出 249,990 × 3,000；UMAP/Leiden 记录 46 clusters。

这些数字只描述对应 run 的验证检查点，不代表最终注释或最终生物学结论。

## 已知限制

- `GSE162498` 尚未进入当前 human integration，需要 extreme-scale 方案。
- `GSE211602` 已被识别为非人数据集并排除。
- gene symbol 标准化和重复基因聚合可能影响跨数据集可比性；必须保留映射和验证。
- 当前尚未在本记录仓库固定 annotation 参数、模型版本和 marker panel。

## 外部证据

- `${RUN_ROOT}/RUN_RECORD.md`
- `${RUN_ROOT}/07_unified_discovery_build/`
- `${PROJECT_ROOT}/Project_v3/run_catalog/run_index.tsv`
