
# 当前权威状态（2026-09-12）

- 当前细胞注释与 CNV 恶性参考：V7，状态 `FROZEN_PROJECT_REFERENCE`。
- 权威 run：`20260905_v7_cnv_rerun_freeze_v1`。
- V8 full-metadata software candidate 的 rerun2 已通过完整性审计，状态
  `PASS_FULL_METADATA_SOFTWARE_CANDIDATE_V7` / `NOT_FROZEN`；它保持 4,676,787
  cells、1,322 samples、54 GSE、43 cancers，以及 128 个冻结 V7 字段合同和 V7
  cell-index hash。该候选仍以 V7 为默认项目参考，不自动晋升，也不允许正式 downstream
  mutation。
- 上述 V8 结论仅为字段、schema、hash、round-trip 和内部一致性的软件候选验证；它不是
  独立准确率、真值、classifier execution 或生物学完成版结论。V7 及既有 V5/V7-derived
  downstream 均未被本候选覆盖。
- scATOMIC retry2 supporting salvage 已记录为 `PASS` 的只读后处理：没有重跑算法、联网、
  安装、覆盖或患者合并，且不能作为扩展资源上界或准确率证据。其最终可复现 supportive audit v2
  为 `PASS_SUPPORTIVE_ONLY`、scope-unverified，且不反序列化 candidate pickle。
- 细胞宇宙：4,676,787 cells、1,322 samples、54 GSE、43 cancers。
- 严格大类可用：4,057,641（86.76%）；严格亚型可用：2,455,424（52.50%）。
- exact Unknown/Ambiguous：129,843（2.78%）；非严格大类总计 619,146（13.24%）。
- primary malignant：515,031；uncertain：484,919；candidate not evaluable：45,284。
- inferCNVpy：938 success、384 explicit skips、0 failed。
- Copykat_python：668 success、654 explicit skips、0 failed。
- CNS 策略 44 samples 中 inferCNVpy/Copykat_python 成功 38/33；黑色素策略
  62 samples 中成功 27/18；其余主要因单样本候选或参考不足。

权威入口：

- [V7 release](releases/annotation_v7_20260905/README.md)
- [V7 run record](runs/20260905_v7_cnv_rerun_freeze.md)
- [V7 method](methods/pancancer-5m-v7-annotation-cnv.md)
- [V8 development candidate index](releases/annotation_v8_development_20260912/README.md)
- [V8 candidate run record](runs/20260912_v8_full_metadata_candidate_v7_rerun2.md)

论文证据链：

- 2026-09-06 通讯—空间—临床证据链已形成报告、汇总表和 10 组 PDF/PNG 成图。
- 该证据链复用 V5 CoVarNet/LIANA，不是 V7 重跑结果；新分析的身份/CNV 默认用 V7。
- [Communication release](releases/communication_evidence_chain_20260906/README.md)

Git 状态：

- 继续维护已有 `Pancancer_log` 历史仓库，不在 TB 级分析父目录初始化第二个竞争仓库。
- Git 仅保存轻量、已审查、可复现/可投稿内容；重型对象由逻辑路径和 SHA-256 管理。
- 当前主机没有 Git LFS，单文件上限按 5 MiB 执行。

---

以下 2026-08 状态作为历史快照保留；涉及当前注释和 CNV 时以上述 V7 段落为准。
# 当前状态

- 最近更新：2026-08-08 06:16 +08:00
- 主要范围：`Project_v3` 泛癌单细胞转录组主线及 `Project_v2` 四队列免疫治疗扩展
- 权威运行：`20260725_120105_unified_pancancer_covarnet`
- 当前阶段：主发现运行已完成且带审计警告；免疫治疗四队列进入最终验收阻断，恶性/通讯 v2 与 24-GSE 下载重试仍在运行
- 运行结论：主运行 `completed_with_audit_warnings`；扩展任务为 `blocked_final_acceptance`、`in_progress`、`in_progress_with_failure`

主运行的详细处理、参数、失败恢复和证据路径见 [正式运行记录](runs/20260725_120105_unified_pancancer_covarnet.md)。

## 三个活跃扩展任务

### 四队列免疫治疗整合

- 记录：[免疫治疗四队列运行记录](runs/20260803_immunotherapy_4cohort_analysis.md)。
- sc/snRNA 元数据、Scrublet、1,440,403-cell Harmony 图谱、43-label 分层注释、GSE316195 严格恶性合并、GSE273952 19-slice 空间整合和 40-factor cell2location reference 均已完成。
- cell2location 已完成 19/19 张切片；71,398 spots、40 factors 精确合并，missing/extra spot 均为 0，finite fraction=1、negative=0、zero-sum=0。
- Phase 6 Stage 13 于 2026-08-06 21:09 失败：22 个 GSE316195 文库中 11 个因样本内参考不足而 `skipped_insufficient_cells`，但 CopyKAT 均成功。上游 malignancy audit 同时记录 `samples_with_method_failure=0`，说明这是最终验收策略与允许跳过语义不一致，不是 CNV 程序崩溃。
- `FINAL_DELIVERABLE_AUDIT.json` 未生成且当前无该分支进程，因此整体状态为 `blocked_final_acceptance_policy_conflict`。

### 恶性识别与 CM 细胞通讯

- 记录：[恶性与细胞通讯运行记录](runs/20260731_malignancy_communication.md)。
- 当前 active v2 manifest 为 721 个生物学样本、1,849,413 个细胞；GSE155698/GSE116256 下游排除和 GSE166555 从 2 个技术组重组为 25 个生物学样本均有机器审计。
- Stage 01 已顺序完成 367/721 个样本（50.9%），第 368 个 GSM5573499/sample34 正在运行 CopyKAT；inferCNV 已成功，R 进程持续高 CPU 并按分钟写 heartbeat。
- 当前 v2 的 Stage 02–08 canonical artifacts 均不存在，会在 Stage 01 后重算；旧 attempt 的阶段历史不作为本轮完成证据。

### 新数据下载

- 记录：[24-GSE 下载与重试记录](runs/20260805_external_geo_download.md)。
- 初始 122-task 下载因 NCBI FTP connection refused 仅 4 verified、118 failed；失败历史保留。
- active Slurm retry v3 为 141 tasks，当前 90 个 terminal rows（63.8%）：76 verified、7 skipped_verified、6 unavailable_upstream、1 failed。
- 唯一失败为 GSE201347 大型 RDS，14,413,317,040/19,071,156,087 bytes 后达到 transfer attempt limit；它不是 `unavailable_upstream`，需单独恢复。
- 第 91 个 GSE274229 PCADT2 matrix 在 06:16 heartbeat 为 101,860,201/194,354,151 bytes（52.4%）；完整性验证前不能计入 verified。Slurm 26090 仍为 RUNNING，数据盘 98% used、约 1002 GiB 可用。

## 主运行已完成

以下为直接由运行产物、脚本和日志核验的事实：

- Phase 06 共 33/33 个数据集通过逐数据集 QC：2,252,118 个 QC 前细胞，2,062,230 个 QC-pass cells，移除 189,888 个细胞，覆盖 797 个 sample units；695/695 个满足条件的样本成功运行 Scrublet。
- human-only 整合计划排除非人数据集 `GSE211602`，形成 32 个数据集、2,040,048 个细胞、39 个可读 archival shards 和 11,063 个共同基因。
- 最终全量拼接对象实际使用 38 个 shards、31 个 GSE、1,963,745 个细胞和 773 个样本；`GSE278694` 在这一过渡中被移除，但现有代码和日志未记录科学理由。
- Harmony、BBKNN、scVI 三种全量整合和 CellTypist 参考标签均已完成；最终对象为 1,963,745 个细胞。
- 证据化注释使用 BBKNN 图，完成 17 个 major lineages、136 个 lineage clusters 和 106 个 final annotations；最终交付审计为 `PASS`。
- CoVarNet 在过滤后 755 个样本上以 `K=9` 完成 nsNMF，得到 106 × 9 的固定权重矩阵和 9 个细胞模块；网络含 101 个节点、374 条边。
- `GSE195832`、`GSE123813` 和 `GSE169246` 均使用冻结的 106-state mapper 与 106 × 9 权重矩阵完成外部投影，三个 `FINAL_AUDIT.json` 均通过。
- 主发现运行最近一次已核验成功写入仍为 2026-07-30 11:18；2026-08-06 的进程检查确认三个后续扩展任务均活跃，不能再概括为“当前无运行中任务”。

## 结果解释边界

- scIB 100k 子集综合分数以 scVI 最高；全图 batch mixing 指标 iLISI 以 BBKNN 最高。不存在所有评价维度上的单一最佳方法，证据注释选择 BBKNN 是下游分支选择，不等于全指标冠军。
- CellTypist 使用 `Immune_All_Low.pkl`，它是免疫参考模型；在混合泛癌对象中的标签应视为辅助证据，不能替代 marker 和人工复核。
- 三个旧外部验证均为固定模块投影而非重新拟合。BH 校正后没有显著模块，当前只能报告探索性方向，不支持确定性的疗效或生存结论。
- `GSE123813` 的主分析仅使用 7 个 sort-matched 配对；全部 15 个配对仅作为受分选混杂影响的敏感性分析。
- `GSE169246` 的 treatment arm/response 来自 legacy metadata，仍需回到原始论文或官方补充材料核验；`Prog` 表示 progression timepoint，不应写成 recurrence 或 confirmed PD。
- 免疫治疗空间分支的 cell2location 结果是 spot-level 多状态丰度；跨癌种参考只支持保守组成解释。
- CM 分支当前只完成逐样本 CNV 的 50.9%，不得把旧 attempt 的 CoVarNet/LIANA 结果写成本轮 v2 完成结果。
- 下载中文件增长只证明传输活跃；只有完整性验证通过才可计为 verified，更不等于已进入核心表达分析。

## 阻塞或审计缺口

- `manifests/stage_status.tsv` 仍把后续阶段标记为 `pending`，与主运行最终产物和审计文件不一致；本记录采用产物审计后的真实状态，源 manifest 尚未改写。
- `GSE278694` 已从 32-dataset human plan 中移除，但现有 notebook/日志只有过滤操作，没有科学理由或 ADR；在补齐依据前应视为未决范围变更。
- `GSE162498` 仍停留在 extreme-scale planning-only 分支，未进入本次最终 31-dataset 主整合。
- 最终注释仍有 56 个 low-confidence、10 个 unresolved 和 5 个 ambiguous clusters，需要人工复核。
- CM active run 的 CopyKAT 实际参数记录为 `genome=hg20`，命名异常，正式发布前需要确认是否为有意约定。
- 免疫治疗最终审计把 11 个可解释的 `skipped_insufficient_cells` 当作 required-method failures，而 malignancy audit 记为 0 method failures；必须先统一验收策略并以新 attempt 重跑 Stage 13。
- 下载分支仍面临数据盘接近容量上限和部分 SRA RunInfo 上游空响应。
- GitHub 远程已配置，但服务器认证仍阻塞普通 push；本地记录仓库可继续维护。

## 下一步行动

1. 继续监控三项扩展任务；只在各自最终 audit/summary 通过后标记完成。
2. 免疫治疗分支明确参考不足的 inferCNV skip 是否为允许的 terminal 状态，记录决策后以新 attempt 重跑 Stage 13；不得改写已通过的 19-slice 丰度结果。
3. CM 通讯分支完成剩余 354 个样本，重算合并、CoVarNet、LIANA、NicheNet，并确认实际 CopyKAT `genome=hg20` 参数是否符合预期。
4. 下载分支继续 141-task retry，为 GSE201347 设计单独可恢复 attempt，完成完整性校验和最终 summary，并在预处理前人工复核 core/sidecar 与各 GSE guardrail。
5. 对齐源 `stage_status.tsv` 与最终审计状态，并补充 `GSE278694` 排除决定及依据。
6. 人工复核 low-confidence、unresolved 和 ambiguous clusters，决定 `GSE162498` extreme-scale 分支，并在认证恢复后普通 push。
