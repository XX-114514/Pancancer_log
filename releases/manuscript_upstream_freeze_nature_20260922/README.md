# Nature manuscript upstream freeze (2026-09-22)

- Release ID: `manuscript_upstream_freeze_nature_20260922`
- Status: `FROZEN_MANUSCRIPT_UPSTREAM_SNAPSHOT`
- Target venue: *Nature*, Article
- Freeze date: 2026-09-22 (Asia/Shanghai)

本 release 固定“现在可用于写稿的注释上游”，并建立 Nature 初投稿文件骨架。它是
投稿证据快照，不是新的全局 V8 生物学发布版。

## 权威边界

1. 全局细胞身份和 CNV 恶性证据只使用冻结的
   [`annotation_v7_20260905`](../annotation_v7_20260905/README.md)。
2. V8 仅有五个谱系的规则/分类学冻结：Myeloid/DC、T/NK、Fibroblast、
   Endothelial 和 B/Plasma。它们没有形成统一的细胞级全局 V8 写回。
3. 2026-09-19 V9 的状态是 `FROZEN_CANDIDATE`；15 个 L3 cluster 和 632 个
   建议重跑 CNV 的样本尚未关闭，因此隔离为候选，见
   [V9_CANDIDATE_QUARANTINE.md](V9_CANDIDATE_QUARANTINE.md)。
4. 当前通讯、空间和临床证据仍来自
   [`communication_evidence_chain_20260906`](../communication_evidence_chain_20260906/README.md)
   的 V5 CoVarNet/LIANA 分支，不能描述为 V7 或 V8 重算结果。
5. V8 full-metadata candidate 仍为 `NOT_FROZEN`，不进入正式结果主分母。

详细冻结规则见 [UPSTREAM_FREEZE.md](UPSTREAM_FREEZE.md)，外部对象、规则和哨兵的
大小与 SHA-256 见 [EXTERNAL_ARTIFACTS.tsv](EXTERNAL_ARTIFACTS.tsv)。V7 的 raw-count
Zarr、身份父对象、坐标、实际 runner、动态 CNV 引擎和最终对象锁见
[V7_PUBLICATION_LOCK.tsv](V7_PUBLICATION_LOCK.tsv)。

## 已准备的投稿文件

- [Nature Article 主稿骨架](manuscript/NATURE_ARTICLE_DRAFT.md)
- [Methods—证据来源映射](manuscript/METHODS_EVIDENCE_MAP.md)
- [初投稿清单](submission/NATURE_SUBMISSION_CHECKLIST.md)
- [Cover letter 模板](submission/COVER_LETTER_TEMPLATE.md)
- [Reporting Summary 工作表](submission/REPORTING_SUMMARY_WORKSHEET.md)
- [Data/Code availability 草案](submission/DATA_CODE_AVAILABILITY_DRAFT.md)
- [主图与 Extended Data 规划](figures/FIGURE_INVENTORY.tsv)
- [逐 panel 来源模板](figures/FIGURE_PANEL_SOURCE.tsv)
- [Supplementary 目录](supplementary/SUPPLEMENTARY_CONTENTS.md)
- [Supplementary table 清单](supplementary/SUPPLEMENTARY_TABLE_INVENTORY.tsv)

## 当前是否可直接投稿

否。文件结构已就绪，但提交前至少要完成或由作者明确降级范围：

- 决定继续完成 V8 12-GSE sealed benchmark，还是以 V9 路线取代；不得混合两套门控；
- 完成 V9 的 15 个 L3 cluster 人工复核、632 样本 CNV 决策和全局晋升审计；
- 全局 V8/V9 是否晋升的正式决定；如不晋升，全文必须按 V7 主参考写作；
- 依赖身份的 CoVarNet/LIANA 在 V7/V8 上重跑，或明确保留 V5-derived 历史分支；
- Copykat_python `genome='hg20'` 参数核验；
- 每个主图 panel 的源表、筛选、统计量、分母和脚本闭环；
- 作者、单位、伦理批准、基金、数据/代码仓库 accession 和利益冲突信息；
- Nature Life Sciences Reporting Summary 官方 PDF。


已在本轮完成 14.2-GB raw-count Zarr 的确定性目录根锁，但 V7 历史精确环境版本、
54-GSE publication source/licence manifest 和 Project_v3 Git commit 仍无法从已有证据
补写，必须在稿件中诚实标为缺失或另行补档。
Nature 当前允许初投格式相对灵活，但建议使用带行号、正文和图合并的单一 Word/PDF。
本 release 按更严格的终稿方向提前组织，官网入口见
[Initial submission](https://www.nature.com/nature/for-authors/initial-submission) 和
[Formatting guide](https://www.nature.com/nature/for-authors/formatting-guide)。
