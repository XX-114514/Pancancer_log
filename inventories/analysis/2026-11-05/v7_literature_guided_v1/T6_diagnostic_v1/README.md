# Tcell_6全子集诊断：聚类完成，marker复核继续

输入为冻结V7中cluster_uid=Tcell_6且major=Tcell的全部156631细胞，95038原生基因；K已逐行精确核验。52个实际来源中至少80%可测的17252基因作为模型候选，2000HVG、50PCA、30图PC；完整RNA文库归一化，不按CellTypist Treg预测预筛选。

已观察到raw res0.4的16群；GSE-Harmony敏感性res0.2/0.4/0.6/0.8分别7/12/15/18群，全部计数156631守恒。这些是诊断分区，不是成熟亚型数量。[各分辨率计数](cluster_resolution_counts.tsv)、[来源集中度](cluster_source_concentration.tsv)、[输入精确范围](input_scope_verified.json)。泛化exporter的uns中scope措辞错误已在K侧车纠正，实际矩阵mask始终正确，原矩阵不改。

K已于2026-11-05 05:11:00 UTC / 13:11:00 UTC+8真实exit0，RNA差异marker、私有metadata和完整对象均已生成；上述群计数快照早于该终端。群11只有2个细胞，未执行差异检验，不能当作成熟稀有亚型。拟合过程中有Quick-TRANSfer上限警告，后续须审阅稳定性及其影响。原CellTypist Regulatory T cells=142908/156631是辅助证据，未作为模型真值。旧winner排序问题与[panel语义审计](../../../../../methods/V7-panel-semantics-audit-20261105.md)和[run](../../../../../runs/20261105_v7_literature_guided_reannotation_v1.md)分开记录。

后续按T身份、CD4/CD8背景、Treg核心、活化/增殖/IFN等独立panel核查，并显示marker对应目标。此子集不覆盖整个T/NK/ILC父群，不能将其结果当作全图谱Treg总量。细胞行和矩阵不公开。

## RNA marker审阅已完成：工作分区与下一步

L2于05:17:09 UTC / 13:17:09 UTC+8完成101基因/1212条marker记录；独立算术审计05:18:22 UTC / 13:18:22 UTC+8 PASS。M于05:20:19 UTC / 13:20:19 UTC+8写出分层工作overlay并精确读回，全部156631细胞保留，旧字段不变。以下是群层面的工作注释，不能将群内每个细胞直接赋为对应亚型。

| 工作群 | 细胞数 | marker与证据驱动的当前处理 |
|---|---:|---|
| K-00 | 49866 | Treg富集T，停在粗粒度；不强分resting/effector |
| K-01 | 28692 | Treg富集T；在当前严格DE门槛下主要差异不足以命名独立亚型 |
| K-02 | 23856 | Treg富集，FOXP3/IL2RA及CCR8/TNFRSF4/TNFRSF9相关激活RNA程序 |
| K-03 | 18818 | CD4背景，CXCL13/TCF7/CXCR5/BCL6相对富集；保留marker描述及Tfh/Tph候选核查 |
| K-04 | 17666 | CD4背景，CXCL13/GZMA/KLRB1及PDCD1相关程序；不直接硬赋Tph或Tex |
| K-05 | 9330 | Treg富集，IFIT1/IFIT3/MX1等IFN状态 |
| K-06 | 3167 | T身份与Ig恒定区RNA；B/浆细胞、环境RNA或双细胞复核 |
| K-07 | 2232 | GZMK/GZMA/CCL5程序，保留T+marker/状态；不凭细胞毒程序判CD8或NK |
| K-08 | 1837 | Ig可变区、来源集中（最大来源82.5%）且浆细胞相关RNA待复核 |
| K-09 | 1027 | T与CSF1R/CD14/MS4A7等髓系程序；污染指标中位0.243，245个联合doublet候选，未确诊 |
| K-10 | 138 | marker支持偏弱，深度与可测性复核；不能据某个相对高表达基因强命名 |
| K-11 | 2 | 差异检验未执行；原记录total_counts中位290981.5，即使跨分辨率不变也不构成新稀有亚型证据 |

- [身份/候选亚型图](../../../../../assets/annotation/20261105_v7_review/T6_K_v5/T6_K_identity_subtype_v5.png)
- [非互斥状态程序图](../../../../../assets/annotation/20261105_v7_review/T6_K_v5/T6_K_state_programs_v5.png)
- [竞争谱系与QC图](../../../../../assets/annotation/20261105_v7_review/T6_K_v5/T6_K_competing_lineage_QC_v5.png)
- [Chu 2023完整14基因Treg功能签名对照](../../../../../assets/annotation/20261105_v7_review/T6_K_v5/T6_K_reference_states_v5.png)
- [marker对应目标和来源](marker_to_type_state_map.tsv)、[可测分母表达表](T6_measured_marker_expression.tsv)、[工作注释](T6_working_cluster_annotations.tsv)、[RNA top10](T6_top10_RNA_DEG.tsv)、[逐群QC](T6_cluster_marker_QC_summary.tsv)、[模型旁证](T6_cluster_reference_top5.tsv)、[深度/QC](T6_cluster_QC_depth.tsv)、[多分辨率重叠](T6_resolution_overlap_summary.tsv)、[独立审计](independent_audit.json)。同名SVG可编辑。

### 为什么不能整群转Treg

K-03/04合计36484细胞，FOXP3检出率约8.1%/10.3%、IL2RA约8.2%/13.8%，而CXCL13约48.5%/50.1%；与K-02的FOXP3约80.6%、IL2RA约71.2%明显不同。CellTypist在K-03/04仍分别将86.0%/72.4%预测为Treg；scATOMIC则以CD4/helper候选为主。模型粒度、训练与共享程序会影响预测，二者都不是本数据的独立真值。这些结果支持拆开检查CD4 helper与Treg，不支持按模型多数票整群转标签。

同细胞阳性规则下，35287细胞有FOXP3及IL2RA/IKZF2/IKZF4至少一个≥1 raw UMI，并同时满足保守T组件支持；提高Treg核心到≥2 UMI为15650。它们是阳性RNA证据数量，受深度/测量范围影响，不能作为Treg总数或以未通过者为Treg阴性。全部T6有845个联合竞争程序+高doublet分位候选，也未据此自动删除。

N于05:22:13 UTC / 13:22:13 UTC+8补充来源测量范围并读回：25753个细胞因至少一组必要marker组件不可测，不能评估上述跨组件阳性规则；另25942个可评估但未达到保守T阳性支持，仍不等于T身份阴性。104936个通过T组件支持，其中35287有Treg核心阳性、69649未通过该Treg正证据规则。未测组件与低正证据已独立编码，不能合并成“缺少T/Treg marker”。[测量可评估性计数](T6_marker_evaluability_counts.tsv)、[N终端](N_complete.json)。

K-05在res0.2/0.6/0.8的最大重叠Jaccard约0.899/0.898/0.832，支持IFN分区的分辨率稳定性；K-00/01/02部分边界较易变化，继续保留粗粒度。该比较使用同一图，不是留出来源复现，也不能消除Quick-TRANSfer警告或保证生物学正确。K-08来源集中、K-11高计数等说明“聚类稳定”本身不足以通过注释。

### 文献适配和验收边界

对照使用[Chu作者代码](https://github.com/Coolgenome/TCM)及其Table S6完整Treg签名；其余panel明确是本项目按文献整理的审阅组合，不冒称作者逐字分类器。Lodi的谱系内重查/QC/重新建模用于分析顺序，不能复制未知Technology或未定义回归参数。本轮未回归IFN/癌种、未按参考Treg预筛选，未硬分Tfh/Tph、resting/effector Treg或CD8耗竭亚型。

下一步先复查K-03/04的CD4 helper边界、K-06/08/09混合RNA以及小群测量/QC，再回到整个T/NK/ILC父群共同重建。当前细胞级overlay仍是可审阅修订层，最终图谱注释未发布；原始V7完整保留。
