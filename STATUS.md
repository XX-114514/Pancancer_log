# 当前状态

- 最近更新：2026-07-31
- 主要范围：`Project_v3` 泛癌单细胞转录组统一分析主线
- 权威运行：`20260725_120105_unified_pancancer_covarnet`
- 当前阶段：全量多方法整合、证据化注释、CoVarNet 发现分析及 3 个固定模块外部验证均已完成；当前无运行中任务
- 运行结论：`completed_with_audit_warnings`

详细处理、参数、失败恢复和证据路径见 [正式运行记录](runs/20260725_120105_unified_pancancer_covarnet.md)。

## 已完成

以下为直接由运行产物、脚本和日志核验的事实：

- Phase 06 共 33/33 个数据集通过逐数据集 QC：2,252,118 个 QC 前细胞，2,062,230 个 QC-pass cells，移除 189,888 个细胞，覆盖 797 个 sample units；695/695 个满足条件的样本成功运行 Scrublet。
- human-only 整合计划排除非人数据集 `GSE211602`，形成 32 个数据集、2,040,048 个细胞、39 个可读 archival shards 和 11,063 个共同基因。
- 最终全量拼接对象实际使用 38 个 shards、31 个 GSE、1,963,745 个细胞和 773 个样本；`GSE278694` 在这一过渡中被移除，但现有代码和日志未记录科学理由。
- Harmony、BBKNN、scVI 三种全量整合和 CellTypist 参考标签均已完成；最终对象为 1,963,745 个细胞。
- 证据化注释使用 BBKNN 图，完成 17 个 major lineages、136 个 lineage clusters 和 106 个 final annotations；最终交付审计为 `PASS`。
- CoVarNet 在过滤后 755 个样本上以 `K=9` 完成 nsNMF，得到 106 × 9 的固定权重矩阵和 9 个细胞模块；网络含 101 个节点、374 条边。
- `GSE195832`、`GSE123813` 和 `GSE169246` 均使用冻结的 106-state mapper 与 106 × 9 权重矩阵完成外部投影，三个 `FINAL_AUDIT.json` 均通过。
- 最近一次已核验成功写入为 2026-07-30 11:18；进程检查未发现相关 Python、R、SLURM worker 或 CoVarNet 任务仍在运行。

## 结果解释边界

- scIB 100k 子集综合分数以 scVI 最高；全图 batch mixing 指标 iLISI 以 BBKNN 最高。不存在所有评价维度上的单一最佳方法，证据注释选择 BBKNN 是下游分支选择，不等于全指标冠军。
- CellTypist 使用 `Immune_All_Low.pkl`，它是免疫参考模型；在混合泛癌对象中的标签应视为辅助证据，不能替代 marker 和人工复核。
- 三个外部验证均为固定模块投影而非重新拟合。BH 校正后没有显著模块，当前只能报告探索性方向，不支持确定性的疗效或生存结论。
- `GSE123813` 的主分析仅使用 7 个 sort-matched 配对；全部 15 个配对仅作为受分选混杂影响的敏感性分析。
- `GSE169246` 的 treatment arm/response 来自 legacy metadata，仍需回到原始论文或官方补充材料核验；`Prog` 表示 progression timepoint，不应写成 recurrence 或 confirmed PD。

## 阻塞或审计缺口

- `manifests/stage_status.tsv` 仍把后续阶段标记为 `pending`，与最终产物和审计文件不一致；本记录采用产物审计后的真实状态，源 manifest 尚未改写。
- `GSE278694` 已从 32-dataset human plan 中移除，但现有 notebook/日志只有过滤操作，没有科学理由或 ADR；在补齐依据前应视为未决范围变更。
- `GSE162498` 仍停留在 extreme-scale planning-only 分支，未进入本次最终 31-dataset 主整合。
- 最终注释仍有 56 个 low-confidence、10 个 unresolved 和 5 个 ambiguous clusters，需要人工复核。
- GitHub 远程已配置，但服务器认证仍阻塞普通 push；本地记录仓库可继续维护。

## 下一步行动

1. 对齐源 `stage_status.tsv` 与最终审计状态，并补充 `GSE278694` 排除决定及依据。
2. 人工复核 low-confidence、unresolved 和 ambiguous clusters，冻结可发布的注释版本。
3. 决定 `GSE162498` 是否通过独立 extreme-scale 分支纳入后续分析。
4. 寻找包含可靠 response/survival endpoint 的真正独立队列，验证冻结的 9 个 CoVarNet 模块。
5. 在 GitHub 认证恢复后执行普通 push；push 前再次运行敏感信息和仓库一致性检查。
