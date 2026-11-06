# 完整轴、未纳入与选择保留

全部原轴4,676,787行保留。既有QC-v4为4,668,312保留／8,237排除／238暂存；本轮没有追加doublet剔除。正式身份冻结3,918,715，占QC保留细胞83.94%；749,597个QC保留细胞未纳入该身份版本，占16.06%。这不是已完成全V7发表级注释的宣告。

| 审计状态 | QC动作 | 细胞数 |
|---|---|---:|
| broad_identity_marker_frozen | primary_QC_keep | 3,918,715 |
| identity_evidence_incomplete | primary_QC_keep | 578,604 |
| outside_primary_QC | QC_holdout | 238 |
| outside_primary_QC | exclude_primary_QC | 8,237 |
| positive_competing_lineage_programmes | primary_QC_keep | 170,637 |
| related_parent_identity_not_separated | primary_QC_keep | 356 |

审计状态只用于说明未纳入原因，不出现在正式身份列中。保留足够大类证据而缺少细分的细胞可以进入正式大类＋marker层；父身份仍混合/低信息的记录不通过改名纳入。正向竞争不等于工具doublet真值；相关谱系共享程序有专门例外。

[54来源保留率](stage_e4_release_retention/source_retention.tsv) · [历史父类保留率](stage_e4_release_retention/historical_parent_retention.tsv) · [来源×历史父类](stage_e4_release_retention/source_by_historical_parent_retention.tsv) · [历史癌种语境](stage_e4_release_retention/historical_cancer_context_retention.tsv)。均为聚合计数，没有患者或样本行。

冻结集合经身份选择，各来源、父类损失不同。它可以用于展示其支持范围内的身份表达，不能未经重新规划分母直接用于无偏细胞丰度、CM发现或替换新TME/CoVarNet输入；原2026年10月CM（v20261010）及更早CM继续标记其真实输入版本。

[全部342条未纳入分组×原因及正向组件汇总](stage_g_outside_identity_diagnostics/outside_identity_group_diagnostics.tsv)。组件可重叠，不是细胞类型百分比。最大未纳入记录为旧Myeloid_boundary/Unknown分区58,020个，其次Unknown分区38,941、Ambiguous分区35,445；这是继续来源内重建的优先级，不是正式身份标签。
