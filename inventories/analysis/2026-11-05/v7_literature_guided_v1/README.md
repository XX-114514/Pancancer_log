# V7逐类注释审计证据

冻结输入SHA256：`e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e`。4676787细胞、26个major、106个最终标签。

后续当前入口：[Mast扩展复核与带目标名称的15群图](Mast_source_review_v1/README.md)；[Tcell_6全子集来源敏感性诊断](T6_diagnostic_v1/README.md)。下列27群是保留的首轮初诊背景，新旧群编号不可互换。

- [106标签原因与处理](V7_106标签原因与处理表.tsv)：含原cluster计数、旧门槛、父类修订、CellTypist大类聚合及当前处理；完整106行计数守恒。
- [父类计数与待重建数](v7_major_reason_audit.tsv)：旧strict字段仅为历史判定，不能代替本轮验收。
- [机器可读逐标签证据](v7_label_reason_audit.tsv)与[旧cluster crosswalk](v7_label_legacy_cluster_crosswalk.tsv)。
- [A阶段终端](stage_a_complete.json)、[B原生输入终端](stage_b_export_complete.json)、[D诊断终端](stage_d_complete.json)。矩阵及细胞级metadata均不在仓库。

## 已完成的Mast首轮复查

全部45036个候选在95038原生基因上重新计算HVG/PCA/邻接图、8个分辨率及RNA marker，Seurat真实exit0。resolution0.4得到27个探索群，来源效应明显，不能称为27种亚型。

| 原V7情况 | 实测原因与证据 | 本轮处理 |
|---|---|---|
| Mast_unresolved，20559 | 两个旧cluster分别98.8%/99.2%细胞至少两个Mast marker阳性，未形成旧规则认可的亚型相对差异 | 支持先确认Mast父类；无需为消灭unresolved强造亚型，仍逐群核QC |
| Mast_coarse_provisional，1408 | raw≥2 UMI时279细胞同时至少两个Mast和两个T marker；存在竞争程序 | 重点双细胞/环境RNA/混群复查，不能直接只改后缀 |
| Mast_coarse_subtype_pending，825 | 父类修订后未重建；污染指标中位较高，来源与表达需复查 | 已纳入新的谱系模型，仍保留QC review |
| 旧Activated两标签，22244 | HDC/LTC4S等本身有广泛Mast表达 | 激活状态须另核完整程序，旧标签不自动当已确认亚型 |

工作性分流是38350父身份支持候选、6686竞争程序/QC待复查候选；不是已通过最终QC的计数，也不是最终排除。cluster16的748个细胞有Mast与增殖程序，作为状态候选。

- [27群review表](mast27_cluster_review.tsv)、[来源可测性修正后的marker表达](mast27_measured_marker_expression.tsv)、[每群top10 RNA DEG及测量范围](mast27_top10_RNA_DEG.tsv)。低来源覆盖marker不能作阴性判据。
- [Mast marker图](../../../../assets/annotation/20261105_v7_review/Mast_marker_review.png)：点大小为可测细胞中检出率，颜色为按基因标准化的log1p平均CP10K；×表示来源覆盖不足80%，*表示身份/QC复查群。
- [C2 Seurat终端](Mast_Seurat_terminal.json)、[C2参数](Mast_Seurat_complete.json)、[F review终端](Mast_review_complete.json)、[完整算术复核](Mast_independent_arithmetic_audit_v2.json)。

## Tcell_6专项核查

Treg的6个panel基因和3个经典锚点已通过，但旧程序先选高均分的TIGIT单基因CD8_exhausted候选，随后因无锚点拒绝；该排序遗漏完整Treg候选。本轮补查841/842样本支持，CellTypist Regulatory T cells=142908/156631；支持Treg富集，不能整群156631细胞直接赋Treg。见[E专项终端](stage_e_T6_candidate_complete.json)及[模型聚合旁证](T6_reference_model_aggregate.json)。

当前未发布新最终标签。481271细胞历史字段要求父类内重新分析；Mast首轮复查和T6候选验证已完成。方法与限制见[执行合同](../../../../methods/V7-literature-guided-reannotation-20261105.md)，具体阶段见[run](../../../../runs/20261105_v7_literature_guided_reannotation_v1.md)。

公共TSV由经过审查的聚合JSON重建，浮点数可能仅显示有限小数位；服务器原始artifact SHA在阶段manifest中按逻辑run路径登记，公共表文件字节以当前Git blob为准。图为服务器原始PNG，SHA精确比对通过。细胞级overlay和所有矩阵均留服务器，没有公共导出。
