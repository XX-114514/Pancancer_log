# V7 旧 panel、命名与判定逻辑审查

本次结论：用户指出的问题成立。旧规则把父类身份、细胞状态、邻近谱系线索和污染结论混为同级候选，且存在先选最高分再验证的算法缺陷。它能解释一部分 unresolved，也会使已经获得 high/medium 的名称表达过度。不能把旧置信度直接解释为生物学正确率，更不能把所有 unresolved 当作父类未知。

本轮完成静态定义审查与聚合证据复算；没有改写 V7 细胞标签，没有重新运行全体细胞注释。修订意见是下一轮注释的约束，不是新参考真值。

## 范围与可追溯产物

- 读取服务器原始 `pancancer_evidence_marker_specs.py`：17 个谱系、162 个候选 panel。
- 读取初始 V5 亚型判定、后续 expanded-panel rescue、大类正负 panel、别名修复和聚类脚本。17 个大类各有正向和排除定义，均逐项审查。
- 162 个候选中有 43 个 contamination panel、28 个 IFN/Cycling panel。并非 162 个互斥细胞类型。
- 基于历史 2,052 条 cluster×panel×gene 聚合证据重建 40 个群的 373 个候选；基因分数、通过条件和 40 个第一候选均与历史结果一致。
- [全部 panel 逐项审查](../inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1/all_panel_reviews.md)：保留旧 marker、旧 canonical、问题及建议角色。另提供 [162 行机器可读表](../inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1/all_162_panel_reviews.tsv)和[17 行父类审查](../inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1/all_17_major_panel_reviews.tsv)。
- [候选复算](../inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1/candidate_reconstruction.tsv)、[逐群规则审查](../inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1/cluster_rule_audit.tsv)、[源码校验和](../inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1/source_manifest.json)。

162 个 panel 是完整配置库，并非全部实际用于 V7 亚型赋值。已核实有效策略中 Epithelial 使用逐细胞 CNV 状态、不做该 marker 亚型判定；Platelet 因 rare-coarse 策略禁用细分。未触发的 panel 缺陷属于配置风险，不计为已发生误标。已量化的候选遮蔽结论仅来自明确复算的 40 群。

## 已确认的主要问题

| 位置 | 旧规则或名称的问题 | 本次处理原则 |
|---|---|---|
| Tcell / NK_contamination | GNLY、NKG7、KLRD1 是共享细胞毒相关基因；不能与 CD8_GZMK 作互斥身份竞赛 | 转为状态/边界证据；TCR/CD3 身份与细胞毒状态分开 |
| Tcell / CD4_naive、CD8_GZMK、CD8_cytotoxic、CD8_exhausted | panel 本身未要求相应 CD4/CD8 身份；IL7R、GZMA、CCL5、TIGIT 等不提供分支特异性 | 先确认 T，再比较 CD4/CD8/非常规分支；最后命名状态 |
| NK / Adaptive_KLRC2 | 把 FCER1G 作为正向 canonical，且混入 TRDC | 暂停直接 adaptive 命名，审查有方向的 adaptive 程序及 γδ T 边界 |
| NK / NKT_like | CD3/TRAC 加细胞毒基因不等于 NKT，更不等于 iNKT | 回到 T/NK 身份边界；TCR/抗原限制证据不足时用描述性候选 |
| B / Naive；Plasma / IgA、Plasmablast | MS4A1 不特异 naïve；JCHAIN 不特异 IgA；MZB1 不区分 plasmablast 与成熟 plasma | 身份、分化阶段、同种型分轴；名称需对应直接证据 |
| B/Plasma、Monocyte/Macrophage 的 contamination | 同谱系分化或共享表达被当作污染 | 区分连续分化、错入父类、真实混合、环境 RNA 与 doublet |
| Granulocyte/DC / Myeloid_contamination | 粒细胞/DC 的共享髓系表达不是通用污染证据 | 改为具体邻近谱系对照，禁止笼统髓系否定 |
| Mast / Activated_IL1RL1 | Mast_2 由 HDC、LTC4S 通过；IL1RL1 未通过仍可套用此名称；activation 主张缺乏独立验证 | 保留 Mast 身份；按实际 marker 描述，激活程序另查 |
| Endothelial / Blood、Venous_ACKR1 | Blood 是更粗层级；PLVAP 可充当 ACKR1 命名的 canonical 替代 | 修正层级；ACKR1 名称须有 ACKR1 特征支持 |
| Fibroblast / CAF、Pericyte、SMC、Mesothelial | CAF 状态/肿瘤背景与邻近谱系混在同一层级 | mural/SMC/mesothelial 单独审查；CAF 需来源和身份支持 |
| Granulocyte / PMN_MDSC | OLR1/ARG1 RNA 不能单独证明免疫抑制功能 | 无功能验证时使用 marker-associated 或 MDSC-like 候选 |
| Neural / Melanocytic | 子 panel 以 PMEL/MLANA/TYR 为支持，但父类规则又排除同一程序；大类偏 Schwann/glia | 独立 melanocytic 谱系；神经/胶质分支及恶性状态另判 |
| 所有 IFN/Cycling | 状态与身份竞争唯一标签 | 状态允许与身份、分化阶段共存 |

这些判断与文献的关键区别有据可查：Qian 的 CD8 群可表达 GNLY 或 KLRD1/FGFBP2，支持“共享程序不能直接视为污染”；不能反过来照搬该文的 NKT-like 解释为本项目已确认 NKT。[Qian 2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf)

经典 adaptive NK 研究描述部分群体 FCER1G/SYK 等降低，故旧 panel 把 FCER1G 正向升高列为 canonical 不合适；这也不意味着所有 adaptive NK 必须 FCER1G 阴性。[Schlums 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4612277/)；[Yang 2019 单细胞 NK 图谱](https://pmc.ncbi.nlm.nih.gov/articles/PMC6718415/)

Ma 的泛癌 DC 工作先用不同分支组合识别 cDC1/cDC2/pDC/LAMP3 DC，再合并细分，而不是将 CD1C 一个共享基因当成 DC3 终点。[Ma 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/)。PMN-MDSC 的功能性名称需按相应鉴定标准解释，不能由转录 panel 自动证明抑制作用。[MDSC 命名与鉴定建议](https://pmc.ncbi.nlm.nih.gov/articles/PMC4935811/)

## 算法缺陷与实际影响

### 1. 先排序，后验收，只验第一名

初始 V5 只用 `canonical` 列表计分；后续 rescue 改为全部 panel 基因计分并要求至少一个 canonical。两个阶段都先选第一名，再检查基因数、样本数与分差，没有对所有候选分别完成验收。不能把两个阶段混称为同一算法。

在本次 40 个 rescue 目标群中，下列 9 个群的第一名不满足旧规则最低基因条件，但后续存在满足“≥2 通过基因且≥1 canonical”的非 contamination 候选：

| 群 | V7 细胞数 | 旧第一候选 | 后续最低基因条件通过候选 |
|---|---:|---|---|
| Dendritic_0 | 5,001 | DC3，仅 CD1C | cDC2 |
| Dendritic_3 | 4,881 | DC3，仅 CD1C | cDC2 |
| Endothelial_4 | 6,158 | Venous_ACKR1，仅 PLVAP | Blood、Tip_angiogenic |
| Endothelial_5 | 15,117 | Venous_ACKR1，仅 PLVAP | Blood、Capillary、Tip_angiogenic |
| ILC_3 | 5,228 | ILC3，无 canonical 通过 | ILC2 |
| Plasma_1 | 12,534 | Bcell_contamination，仅 CD37 | Plasmablast、Mature_SDC1、IgG |
| Plasma_2 | 69,779 | Bcell_contamination，仅 CD79A | Plasmablast、Mature_SDC1、IgA |
| Plasma_3 | 48,340 | Plasmablast，仅 MZB1 | Mature_SDC1 |
| Tcell_6 | 156,631 | CD8_exhausted，仅 TIGIT | CD4_Treg、Tissue_resident、CD4_Tfh、IFN_T |

合计 **323,669 个细胞**。这是“需要复核的规则暴露范围”，不是纠正了 323,669 个细胞。最低基因条件仍使用有缺陷的旧 canonical，本身不保证生物学合理；本轮未为所有备选重算样本复现，也不按第二名自动改标签。Tcell_6 的 Treg 样本与模型复核已有上一轮独立证据，不能外推到其他八群。

### 2. 通过基因均值会偏爱窄信号

旧分数为通过门槛基因的 `pct_in × max(log_difference,0) × max(delta_pct,0)` 的均值。未通过基因完全不进入均值，所以一个高分共享基因可能压过多个一致但中等强度的身份基因。增加一个中等支持基因甚至能降低均值。canonical 也只是手工名单，不代表特异性得到验证。

Tcell_6 正是实际反例：TIGIT 单基因的 CD8_exhausted 为 0.476662；FOXP3、IKZF2、IL2RA、TNFRSF4、CTLA4、TIGIT 全部通过的 CD4_Treg 为 0.328864。程序先选前者，随后因其 canonical 不足拒绝命名。此前若把原因概括成“这个群缺少 canonical”，应更正为“被选中候选缺少 canonical，其他候选未被完整验收”。

### 3. 伪互斥分差与未校准门槛

Tcell_4 的 CD8_GZMK 0.697529 与 NK_contamination 0.690298 差 0.007230，小于旧 0.015。这只是两个共享表达 panel 的分差，不是 T/NK 身份概率差。已核验群内 91.0% 细胞检测到 CD3D/CD3E/TRAC 至少两个，且 59.1% 同时有该 T 模块及至少两个 NKG7/GNLY/KLRD1；因此不能把共享细胞毒表达直接解释为两类细胞混合。

0.010/0.015 随谱系丰度分档，是项目启发式阈值，并未证明是论文采用或误注率校准后的门槛。修订方向不是把 0.015 降低，而是先移除不该互斥的比较。

### 4. “未富集”不等于“未表达”，也不等于“无身份”

旧门槛比较一个 cluster 与同父类其余细胞。广泛的谱系身份基因可能在两边都高表达，因差异不足不通过。40 群中有 15 群所有候选均无通过基因；排序仍保留配置顺序的第一项。这个 top_subtype 不是可信的生物学第一候选，应显示 `no_differential_candidate`。

旧字段 `log2fc` 实际是 `(mean(log1p CP10K)_in − mean(log1p CP10K)_out)/ln(2)`，不是线性均值表达比的 log2 fold change。后续必须正确命名统计量，不能直接套用论文采用另一种 logFC 的数值阈值。

### 5. 名称、high/medium 与父类门槛均需独立审查

初始 V5 的后缀 marker 取通过列表第一个，受 panel 排列顺序影响，不是差异排序第一。high/medium 根据样本支持数与分差生成，不是准确率评估。后续 rescue 以名称是否含 contamination 决定是否拒绝，说明行为绑在字符串上；初始 V5 该段代码反而没有同样的显式 contamination 拒绝，两个阶段语义不一致。本次初始结果表未发现成功发布的 contamination 标签，不把此代码风险写成已发生误标。

大类代码注释声称单 marker 不足，但实际 `len(present)≥2` 检查的是该来源存在多少 marker 列，并不是每个细胞检测到多少 marker；随后仍取表达均值。全局大类得分、排除得分、原注释与 CellTypist 组合均需核验，不能以代码注释替代实际判断条件。反过来，父身份也不能只由一次全局 cluster 的标签继承到所有细胞。

## 数据覆盖和别名：本轮不把它们误报为已证实主因

162 个 panel 的请求基因经旧 KIM1→HAVCR1 映射后，全部存在于历史 359-gene 对象的列中；没有本轮 panel 基因被旧 blacklist 实际移除。**这只证明对象有列，不证明每个来源都测得。**旧对象使用 per-GSE 可测性 mask，未测得不能按阴性解释。

还核实 `04b_expand_marker_aliases.py` 已处理 SELENOP↔SEPP1、ACKR1↔DARC、JCHAIN↔IGJ、CXCL8↔IL8 的来源特异回退。因此不能因为 rescue 脚本只写 KIM1 一个映射就声称其他别名都没处理。后续原生全基因重注释仍需 stable-ID/版本审计，禁止随意剥离后缀合并不同 feature。

## 修订后的工作约束

下列是本项目拟采用的结构化审查约束，并非宣称某一篇论文逐字使用这些字段。

1. **父身份**：`parent_identity`、`parent_support`、`parent_conflicts`。使用多基因组合、同细胞表达、来源可测性及模型辅助；边界不明才进入混合/父类待定。
2. **分支/分化身份**：`identity_subtype` 与 `subtype_status`。同一父类内比较真正可区分的身份候选；先验收每个候选，再判断冲突。多个候选都不够时保留父类。
3. **状态**：`state_programs` 可多值，保留 cytotoxic、IFN、cycling、inflammatory、resident-like 等，与身份共存。
4. **质量与混合**：`qc_review`、`admixture_hypothesis`。异源 marker 仅触发复核；结合表达是否来自同一细胞、独立亚群、环境 RNA、doublet 与来源偏倚，不自动删除或改父类。
5. **命名**：marker 后缀要求该 marker 的真实表达/相对特征支持，并报告判据。描述性 marker-enriched 群可用，但不能伪装成已确认功能或成熟亚型；不要求每个细胞都有该 marker。
6. **终点**：区分 `resolved_at_major`（大类是合理终点）、`subtype_pending`（亚型证据待核）、`parent_ambiguous`（父身份冲突）、`qc_pending`。红系、HSPC、ILC 不是必须细分，也不是天然不能细分；证据支持到哪一层，结论停在哪一层。HSPC/ILC 尤其要先确认大类本身。

Lodi 的公开脚本包含同谱系提取、QC/doublet 检查、重新建图、找 marker 和人工 RenameIdents，但亚型 marker 示例存在占位内容；不能称旧固定 panel/阈值是严格复现 Lodi。[Lodi 作者代码](https://github.com/lambrechtslab/Pan-cancer_Lodietal/tree/main/Cell%20subtype%20annotations)。Chu 的 T 细胞研究可用于 CD4/CD8/非常规与状态分层参考，仍需在本项目的表达与来源中验证。[Chu 作者代码](https://github.com/Coolgenome/TCM)

下一步先处理 T/NK/ILC 的共同边界：冻结历史标签作为对照，审查 Tcell_4 和 Tcell_6，并同时抽查已有 `CD4_naive_IL7R` 与 `CD8_cytotoxic_CCL5` 群。这些名称也没有被旧 panel 直接证明 CD4/CD8。必要时在确认身份的子集内重建图/聚类；不因降低门槛而“清空 unresolved”。最终对照原标签、人工证据、CellTypist/其他参考、跨来源稳定性与混合细胞比例，才发布新版本。
