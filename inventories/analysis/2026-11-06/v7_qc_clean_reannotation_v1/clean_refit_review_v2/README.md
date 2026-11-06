# V7实际清理后：Mast与T6大类工作注释及marker复查

本轮已经实际生成排除/暂存名单、从原始计数重新建模，并写出逐细胞工作注释侧表且读回一致。最新QC为v2；最初25文件提交84adfc6记录的是v1启动节点，原证据仍保留。**这里完成的是旧V7 Mast/T6候选池的大类工作判断与QC/状态复查，不是全V7细分注释或发表终验。**

| 范围 | 原n | 主分析QC剔除 | QC暂存 | 本轮模型n |
|---|---:|---:|---:|---:|
| 冻结V7 | 4,676,787 | 2,003 | 236 | 全局保留视图4,674,548，未重拟合全局 |
| T6旧候选池 | 156,631 | 783 | 62 | 155,786 |
| Mast旧候选池 | 45,036 | 1,186 | 173 | 38,923；另4,754个QC保留细胞待父类／竞争RNA复核 |
| 旧K-09（包含于T6） | 1,027 | 227 | 18 | 782 |

2,003＝783 T6＋1,150 Mast初轮＋36 Mast追加NK组件联合高分＋34原未评分单元补测。236＝62 T6＋173 Mast＋1补测数值警告候选。它们是操作性QC排除／暂存，不是有生物学真值验证的doublet。没有整簇按配额删除，也没有删除原文件或改写冻结V7。

## unresolved与coarse的原因已经分开

| 原因 | 本轮直接证据 | 已采取处理 |
|---|---|---|
| 候选仅打标签，未进入分析掩码 | 原T6/Mast联合候选仍在旧对象 | 已建立新掩码并用于实际计数子集与重新建图 |
| TCR在输入轴上不可用 | 7来源／25,753个T6细胞没有TRAC/TRBC1/TRBC2，保留CD3D/E/G | 保留unknown测量状态，补查CD3；其中16,635个有两个CD3基因各≥2 raw UMI |
| 把逐细胞保守阳性规则当作大类必需条件 | 补充CD3后单细胞T阳性证据由104,091升至132,698（clean T6） | 该数值仅作证据；大类按簇多marker和来源范围人工复核，不要求每个细胞满足同一阈值 |
| 异源RNA／剩余混杂 | T07恒定区Ig、T08可变区Ig、T09髓系RNA；Mast部分簇带淋巴或组织RNA | 保留逐细胞竞争程序／来源／污染率，单独标复核，未把上下文自动叫作新亚型 |
| 分辨率产生分区，但没有充分细分依据 | 当前Mast M01/M05/M07无符合本设置的正向DE；T6 T00/T07/T10/T11同样如此，T11仅2细胞 | 保留coarse或复查，不为每个计算簇强造亚型 |
| 状态与身份混用 | IFN、增殖、热休克等可跨亚型 | 分开记录大类、候选细分、状态关联与QC |

CD3补充规则为原CD3＋TCR阳性，或≥2个CD3D/E/G各≥2 raw UMI；这是项目正向诊断条件，不是从论文复制的统一分类器。原文[Lodi代码](https://github.com/lambrechtslab/Pan-cancer_Lodietal/blob/main/Cell%20subtype%20annotations/Tcell_analysis.R)查看CD3D/CD3E、CD4/CD8A及竞争谱系marker，没有统一TRAC必需门槛。来源缺测既不代表生物学不表达，也不直接决定Treg阴性。

## 已落实的大类工作注释

- T6：154,974个获得簇级多marker支持的T大类工作标签；662个保留T／髓系竞争RNA复查，150个为低信息T候选。154,974中仍包含Ig复查群和2细胞的小分区，不能解读成全部通过精细注释终验。单细胞132,698个直接阳性证据另列保留，未把其余细胞一概叫unresolved。
- Mast：38,923个获得Mast大类工作标签；其中2,150个处于额外RNA／来源上下文复查簇。原候选池另4,754个尚需父类／竞争RNA复核，不从其他父类招募，所以这也不是全V7最终Mast普查。
- 所有工作侧表保留原V7、CellTypist和scATOMIC字段。分类器仅作参考，没有作为本轮建图输入或最终答案。细分最终标签本轮仍未接受，cluster状态文字不是逐细胞功能实验结论。

## T6全部12个当前诊断簇

| 新簇 | n | 工作注释／表达背景 | 复核marker（含非DE） | 最大来源占比 |
|---|---:|---|---|---:|
| T00 | 52,004 | Treg-enriched T (coarse) | CD3D, CD3E, FOXP3, IL2RA | 11.2% |
| T01 | 24,979 | Treg-enriched ／ IKZF2 / FOXP3 | IKZF2, FOXP3, IL2RA | 17.2% |
| T02 | 21,477 | Treg-enriched ／ activation context | FOXP3, IL2RA, TNFRSF4, TNFRSF9 | 19.7% |
| T03 | 20,531 | CD4 T ／ CXCL13 / TCF7 context | CD4, CXCL13, TCF7, GNG4 | 17.2% |
| T04 | 19,120 | CD4 T ／ CXCL13 / GZMA context | CD4, CXCL13, GZMA | 14.0% |
| T05 | 9,569 | Treg-enriched T ／ IFN context | FOXP3, IL2RA, IFIT1, MX1 | 20.8% |
| T06 | 2,700 | T ／ GZMK / CCL5 context | GZMK, GZMA, CCL5 | 15.9% |
| T07 | 2,371 | T + Ig-constant RNA [review] | CD3D, CD3E, IGHG1, IGLC3, JCHAIN | 23.2% |
| T08 | 2,221 | T + Ig-variable RNA [review] | CD3D, CD3E, JCHAIN, IGLV7-43, IGLV2-11 | 80.4% |
| T09 | 662 | T + myeloid RNA [review] | CD3D, CD3E, FCGR2A, CD14, SPI1 | 11.2% |
| T10 | 150 | T candidate ／ low support [review] | CD3D, CD3E | 30.0% |
| T11 | 2 | T ／ tiny partition [review] | CD3E, TRAC | 50.0% |

T07的IGHG1、IGKC检出率分别约73%、96%；T08的IGLV2-11约66%，且80.4%细胞来自单一来源。Ig RNA不等于B细胞身份；须结合MS4A1/CD79A或MZB1/DERL3等独立组件。旧K-09保留782个：533个进入当前T09，249个分散到其他簇。当前T09总662个，其中186个仍触发髓系程序、172个有髓系锚点，不能作为已确认的T亚型用于生物学解释。

当前编号没有继承旧编号的生物学含义。当前T11的2细胞来自旧K00，原旧K11高深度细胞去向另见[完整转移审计](receipts/stage_k_model_transition_audit_v1/T6_old_K09_K11_to_clean_clusters.tsv)。

## Mast全部18个当前诊断簇

| 新簇 | n | 工作注释／表达背景 | 复核marker（含非DE） | 最大来源占比 |
|---|---:|---|---|---:|
| M00 | 9,945 | Mast ／ heat-shock / IEG-enriched | HSPA1A, HSPA1B, DNAJB1 | 27.3% |
| M01 | 9,219 | Mast coarse | TPSAB1, TPSB2, CPA3, KIT | 23.2% |
| M02 | 6,266 | Mast ／ HLA-II-enriched | HLA-DPA1, HLA-DRA, CD74 | 13.5% |
| M03 | 3,330 | Mast ／ NF-kB-associated | BCL2A1, NFKB1, RELB | 27.9% |
| M04 | 2,562 | Mast ／ TXNIP context | TXNIP | 14.9% |
| M05 | 1,588 | Mast coarse | TPSAB1, TPSB2, CPA3, KIT | 24.7% |
| M06 | 1,166 | Mast ／ IFN-enriched | MX1, IFIT1, ISG15 | 10.7% |
| M07 | 1,012 | Mast coarse | TPSAB1, TPSB2, CPA3, KIT | 14.8% |
| M08 | 948 | Mast + lymphoid RNA [review] | CD2, IL7R, KLRB1, NKG7, CCL5 | 32.9% |
| M09 | 762 | Mast coarse ／ PKM context | PKM | 27.6% |
| M10 | 652 | Mast + inflammatory RNA [review] | S100A8, S100A9, LYZ, IL1B | 12.9% |
| M11 | 614 | Mast ／ cycling-enriched | CCNB2, TOP2A, UBE2C, MKI67 | 14.7% |
| M12 | 309 | Mast ／ ITGAV / NFKB1 context | ITGAV, NFKB1 | 21.0% |
| M13 | 234 | Mast + acinar RNA [review] | PRSS1, PLA2G1B, CEL, CTRC | 48.3% |
| M14 | 111 | Mast ／ source-only partition [review] | MATR3, CHST2, EGR1 | 100.0% |
| M15 | 94 | Mast + endocrine RNA [review] | CHGB, PCSK1N, CPE | 100.0% |
| M16 | 69 | Mast ／ heat-shock / source [review] | HSPA6, BAG3, HSPA1A | 87.0% |
| M17 | 42 | Mast + neural RNA [review] | MAPT, TACC2, CACNA1A | 100.0% |

Mast追加36个排除使用不同NK受体（KLRD1/NCR1）、效应组件（GNLY/PRF1）和NKG7各≥2 raw UMI，同时原样本／方法内分位≥.95且校准有效。其余宽泛淋巴上下文不自动删除；[pan-MC原代码](https://github.com/xiaoyanger01/pan-MC/blob/main/Fig.%201.R)的参考C9本身含CCL5/IL32/GZMB，因而不能只见细胞毒基因就判doublet。当前M08仍需复核，未称新的细胞毒Mast亚型。

## 图怎么展示

图的上方均标对应身份／候选亚型／状态，下方是marker，左侧为工作注释＋n。点大小使用来源可测细胞分母；覆盖<80%标×；颜色是基因内相对表达，不能跨基因直接比较。10张PNG全部经过模型范围、SHA、SVG文本读回及人工看图；[图核验](figure_readback_audit.json)。

| 角色 | T6 | Mast |
|---|---|---|
| 大类与候选细分，报告正文必要 | [身份图](figures/T6_clean_identity_subtype_v7.png) | [身份图](figures/Mast_clean_v2_identity_subtype_v7.png) |
| 状态与身份分开，报告正文必要 | [状态图](figures/T6_clean_state_programs_v7.png) | [状态图](figures/Mast_clean_v2_state_programs_v7.png) |
| 标准竞争谱系检查，保留为QC证据 | [QC图](figures/T6_clean_competing_lineage_QC_v7.png) | [QC图](figures/Mast_clean_v2_competing_lineage_QC_v7.png) |
| 本次混杂原因，复查报告必要 | [Ig上下文](figures/T6_clean_extended_RNA_context_v7.png) | [扩展RNA上下文](figures/Mast_clean_v2_extended_RNA_context_v7.png) |
| 文献状态对应，仅作补充 | [参考状态](figures/T6_clean_reference_states_v7.png) | [参考状态](figures/Mast_clean_v2_reference_states_v7.png) |

## 验收与边界

新掩码互斥守恒、全原始计数子集／10x精确一致；T6与Mast-v2终端均0；HVG与DE测量范围、全部3,396条marker记录和全部分辨率／来源计数已独立重算。原冻结V7的SHA再次核到不变。逐细胞侧表201,667行覆盖两个完整原候选池，已写出并读回一致。[侧表收据](receipts/stage_h_working_review_overlays_v1/complete.json)／[当前分析视图合同](receipts/stage_g_Mast_NK_cleanup_v2/analysis_view_contract.json)。完整细胞行、样本标识、矩阵及R对象保留私有。

T6保留10条Harmony初始化Quick-TRANSfer警告；稀有2细胞分区无法DE；Mast-v2有未过min.diff.pct的分区。计算完成、统计读回通过、数值敏感性与生物学身份终验分开，不能称已达发表级别精细图谱。全基因输入仍是raw counts，没有假定存在decontX校正层；污染率不等于校正后的表达。其他谱系及其余T候选池尚未完成本轮逐类复审。

下一步先复核T09、T07/T08与Mast剩余父类／上下文候选；在可测范围内确认身份后再决定Treg/CD4等细分，并补来源／参数敏感性。现有CoVarNet输入未自动替换；新TME数据工作线与本旧V7展示／排查分开。
