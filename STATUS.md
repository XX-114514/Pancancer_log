# 当前状态

- 最近更新：2026-07-26
- 主要范围：`Project_v3` 泛癌单细胞转录组统一分析主线
- 当前阶段：Phase 07 统一发现队列整合已完成，细胞注释待进行
- 权威外部运行：`20260725_120105_unified_pancancer_covarnet`

## 已完成

以下为从最新运行记录直接观察到的事实：

- Phase 06 已验证 33/33 个数据集通过逐数据集 QC，共 2,062,230 个 QC-pass cells、797 个 sample units。
- Phase 07 将非人数据集 `GSE211602` 排除出 human-only 整合分支，形成 32 个数据集、2,040,048 个细胞的整合输入。
- 32 个 human datasets 经基因名标准化后共有 11,063 个共同基因。
- 39 个 archival shards 已全部验证，覆盖 2,040,048 个细胞。
- capped modeling branch 已生成，包含 249,990 个细胞和 11,063 个基因。
- HVG/PCA/Harmony 已成功生成 249,990 × 3,000 的模型对象。
- UMAP/Leiden 已成功完成，记录为 46 个 Leiden clusters。

## 正在进行

- Phase 07 整合后细胞类型注释方案与执行。
- 将新的正式运行按本仓库模板建立轻量级记录和结构化索引。

## 阻塞或边界

- `GSE162498` 仍属于 extreme-scale pending 分支；已有分块计划，但最新检查中未观察到其进入上述 32-dataset human integration。
- `GSE211602` 被识别为非人数据集，已从 human-only 主整合分支排除；后续用途待确认。
- 当前记录仓库已配置用户指定的 GitHub 远程；首次 push 状态见初始化 Codex 日志。
- 当前模型结果尚不能据此声称完成最终细胞注释或得到最终生物学结论。

## 最新验证检查点

外部证据：

- `${PROJECT_ROOT}/Project_v3/runs/20260725_120105_unified_pancancer_covarnet/RUN_RECORD.md`
- `${PROJECT_ROOT}/Project_v3/run_catalog/run_index.tsv`

最近可见检查点为 Phase 07 model branch UMAP/Leiden：状态 `success`，对象形状 249,990 × 3,000，46 个 Leiden clusters。

## 下一步行动

1. 为 Phase 07 annotation 建立正式 run record，并冻结注释输入、marker evidence、模型版本和完成标准。
2. 确认 `GSE162498` 的 extreme-scale 分支是否进入后续主分析。
3. 核验 GitHub 远程为 private 并完成首次普通 push。
4. 每次有意义任务同步更新 Codex log；只有当前状态改变时才更新本文件。
