# TME 从头注释 v2：按 Lodi et al. 2025 设计（2026-11-05）

本文件取代 [v1 方案](tme-de-novo-annotation-v1-20261105.md)。v1 保留为历史版本，簇13的 v1 探索性运行只作诊断，不作为最终注释。

用户要求（2026-11-05）：
- 借鉴 V7 panel 并与本项目 v1 panel 合并，给出文献依据；
- 上传完整的注释基因面板和逻辑；
- 严格学习 Lodi et al. 的代码和论文设计注释逻辑；
- 先从头注释，再与2026年10月 CM（v20261010）的标签对照。

对象：TME v1 run v2（`${PC_ROOT}/tme_v1/run_v2/tme_v1.zarr`）。共3,347,659个细胞，CoVarNet v1 纳入2,980,723个。

## 依据

| 来源 | 本文采用的内容 | 核实方式 |
|---|---|---|
| Lodi F et al., *Cell Rep Med* 2025, doi:10.1016/j.xcrm.2025.102416（PMID 41086810） | Methods：QC、各癌种分别确定大类、跨癌种合并后分谱系聚类、回归变量、低质量/双细胞簇判据、丰度阈值 | 2026-11-05 通读 PMC12629804 Methods |
| [Lodi 作者代码](https://github.com/lambrechtslab/Pan-cancer_Lodietal) @35f054b | `Cancer type analysis.R`（大类）、`Cell subtype annotations/*.R`（T、B、DC、EC、Mac/Mono）、`Master_files/marker_genes.txt`、`Minor_2.0_Intermediatecelltype_annotation.csv` | 2026-11-05 克隆后逐文件阅读；注册表中每项都记录了文件和行号 |
| [V7 panel repair v1](../resources/annotation/v7_panel_repair_v1/README.md)（项目内，2026-11-05） | 112个证据模块：parent/branch/state/isotype/boundary 分层、鉴别对照、命名 marker、器官条件、逐基因依据 | 原样并入，未修改 |
| [合并面板 v2](../resources/annotation/tme_annotation_panel_v2/ALL_PANELS.md) | 179个模块 = V7 112 + Lodi 作者代码20 + Lodi 发表的 DEG 谱21 + 本项目 v1 补充26；共651个基因、49条文献 | 用 V7 `evidence_engine.validate_registry` 校验通过；V7 测试15/15 通过 |

文献核实（2026-11-05）：
- v2 新增文献均在 PubMed esummary 或 Crossref DOI 上核对了标题和作者。
- Lodi 代码中写的"Chen et al, Cell 2021"（肺泡巨噬）没能对应到具体论文，记为 unresolved，不作为独立依据。
- v1 补充模块的文献只核实到模块层面，逐基因出现证据待核，状态记为 `module_reference_verified_gene_level_pending`。
- Lodi 面板与 DEG 谱逐字取自代码和表格，状态为 `author_code_exact_gene_list` / `author_DEG_table_exact`。

## Lodi 的设计（照原文和代码整理）

1. **样本级 QC（在确定大类之前）**
   - 过滤：基因数 >200 且 <6000，UMI >400，MT <25%。
   - 癌种例外：GBM UMI >600 且 MT <20%；HCC MT <50%；HNSCC 和早期 BC MT <15%。
   - 每个样本跑 DoubletFinder（nExp = 3.9%，pN 0.25，pK 0.01，PC 1:20）。只记录分数，不按单细胞删除，留到谱系内判断双细胞簇时使用。
2. **各癌种分别确定大类**
   - 流程：LogNormalize 1e4；vst 2000 HVG；ScaleData 回归样本 ID、UMI、MT%、S/G2M；PCA，PC 数按肘部图选；多分辨率 Louvain。
   - 差异基因：`FindAllMarkers(only.pos=TRUE, min.pct=.25, logfc.threshold=.25)`。
   - 用29个标准 marker 人工 `RenameIdents`。
   - 肥大细胞和 DC 常与髓系共聚，需在更高分辨率下分出。
3. **癌种特有群体先在单个癌种内识别，并在跨癌种合并前移除**：HCC 的 Kupffer 细胞和肝窦 EC，GBM 的小胶质细胞，NSCLC 的肺泡巨噬和肺泡 EC。所用基因见注册表 `tissue_specific_population`。
4. **按谱系跨癌种合并，两轮分析**
   - **第一轮（QC 轮）**
     - 特征与回归：vst 2000 HVG；ScaleData 回归样本 ID、UMI、MT%、S/G2M、癌种、IFN（MSigDB BROWNE），T/B/EC 另加 stress，Mac/B 另加 hypoxia。
     - 嵌入与聚类：PCA 按肘部图取 PC；Harmony 只校正 Technology（3′/5′）；多分辨率 Louvain（0.2–2.8）。
     - 差异基因：`FindAllMarkers(min.pct=.10, min.diff.pct=.25, logfc.threshold=.25, return.thresh=.01, only.pos=FALSE)`。
     - 每簇检查29个大类 marker（他谱系污染）和 DoubletFinder 预测。簇被命名为谱系、`Low quality`、`Doublets`、`IFN` 或 `Proliferative`。
     - 低质量：UMI/基因数低、MT 高、几乎没有基因表达。
     - 双细胞：表达其他谱系的 marker，且 DoubletFinder 双细胞率高于 10x 预期（每回收1,000个细胞0.8%）。
     - 删除低质量和双细胞簇。
   - **第二轮（干净对象重建）**
     - 3000 HVG，去掉血红蛋白基因（B 细胞另去掉 IGHV/IGKV/IGLV）后取前2000。
     - 同样回归；PCA 50，取40个 PC；Harmony（Technology，`epsilon=-Inf`，不提前停止）。
     - 分辨率网格 0–6；在主分辨率按经典 marker 人工命名；稀有亚型在更高分辨率下补标（例：T 细胞在0.7命名，γδ 用1.8、CD8 TN 用4.5、MAIT 用3.6、TFH 用2.2）。
     - 再删除低质量簇。
   - **B 细胞例外**：Harmony 不足以去除批次，改用 Seurat CCA（按 Technology 拆分）。
   - **用 LISI 确认批次已去除。**
5. **不做跨癌种细分**：上皮/肿瘤细胞（"Cancer/epithelial"）和成纤维细胞只作为大类。器官特异群体按单个癌种注释。
6. **状态簇**：增殖、IFN、缺氧、MT 高的簇作为命名状态保留（例：Prolif T、IFN Mac、Hypoxic Mac、MT Mac、IFN EC、Prolif EC）。
7. **下游丰度阈值**：
   - 只纳入细胞数 >500 的样本（第10百分位）；
   - 亚型丰度分析要求每个样本 T 细胞和巨噬/单核 ≥20 个，B、DC、EC、CD8 T 或 Treg ≥10 个；
   - 亚型比例以所属大类为分母，例如 CD4 T / 全部 T/NK。

## 本项目的执行流程（v2）

| 步骤 | 做法 | 与 Lodi 的差异和理由 |
|---|---|---|
| A0 样本级 QC | 沿用 integration v2 的 QC 和 decontX（污染 <0.5）。新增：对全部样本补算统一的双细胞分数。R 环境没有 DoubletFinder，先用 Scrublet（已安装）或安装 DoubletFinder；参数和预期率按 10x 指南写入运行记录 | 现有来源的双细胞判定继承自各原始 store，口径不一，因此要统一重算一个只作旁证的分数。不按单细胞删除，与 Lodi 一致 |
| A1 恶性细胞 | 恶性细胞已用 scATOMIC + inferCNV 共识排除；uncertain 留在 TME 并带标记，不进入 CoVarNet v1 | Lodi 把肿瘤和上皮合并为一个大类；本项目研究对象是 TME |
| A2 各癌种大类 | 对28个 `analysis_main` 分别做 vst 2000 HVG、ScaleData 回归（样本、UMI、MT、S/G2M）、PCA（肘部图）、多分辨率聚类，FindAllMarkers 参数与 Lodi 大类一致 | 与 Lodi 相同。全局 run v2 的区室标签只作对照 |
| A3 大类命名 | 用 Lodi 的29个 marker（`Lodi_major_*`）加 V7 parent 模块人工命名；肥大细胞和 DC 在高分辨率下从髓系分出；上皮簇用器官条件模块确定组织来源（v1 补充模块和 V7 的 Epi_*） | 与 Lodi 相同，另加 V7 的边界对照 |
| A4 癌种特有群体 | 在单个癌种内识别：肝癌（HCC、iCCA）的 Kupffer 和肝窦 EC，胶质瘤的小胶质细胞，肺癌（NSCLC、LUAD）的肺泡巨噬和肺泡 EC。识别后从跨癌种谱系池中移出，作为独立类型保留 | 与 Lodi 相同；本项目的肝癌和肺癌各包含多个子类型 |
| B1 谱系池 | T/NK、B/浆细胞、巨噬/单核（含中性粒）、DC、EC 按 Lodi 跨癌种合并 | 与 Lodi 相同 |
| B2 第一轮 QC | 2000 HVG；回归变量按 Lodi 各谱系脚本；Harmony；多分辨率；FindAllMarkers 参数按 Lodi 谱系脚本；逐簇检查29个大类 marker 和双细胞旁证分数；人工判定 LQ/Doublets/IFN/Prolif | **批次变量用来源而不是 Technology**：本项目70个来源来自不同实验室，去批次比较 v1 显示 Harmony 按来源 θ2 最好（同来源近邻0.154 对比未整合0.628）。化学（3′/5′/BD/Singleron/Drop-seq）嵌套在来源内，另记录并用 LISI 检查 |
| B3 回归的实现 | 样本 ID 作为分类回归量，按样本对每个基因去均值（Frisch–Waugh：与连续协变量联合 OLS 等价）；再回归 UMI、MT、S/G2M、IFN（BROWNE，MSigDB 2025.1.Hs，67个基因）、hypoxia（Lodi 15基因）；癌种嵌套在样本内，样本去均值后已被吸收 | 1,200个样本的哑变量在百万级细胞上做 OLS 不可行。数学上与 Lodi 的回归等价；用小样本与 Seurat ScaleData 对比验证 |
| B4 stress 评分 | Lodi 的 `stress.genes.csv` 不公开。使用 V7 的 `Stress_response` 模块（HSPA1A、HSPA1B、DNAJB1、HSPH1、HSP90AA1），作为替代并明确标注 | 不能声称与 Lodi 的 stress 分数相同 |
| B5 第二轮 | 3000 HVG，去掉血红蛋白基因（B 细胞另去 IG V 基因）后取前2000；PCA 50 取40；Harmony 不提前停止；分辨率 0–6；在主分辨率命名并在高分辨率补稀有亚型；用 LISI 检查 | 与 Lodi 相同。B 细胞若 Harmony 后 LISI 仍不足，按 Lodi 的做法改用 Seurat CCA（scrna_r_full 中有 Seurat 5.3.1） |
| B6 命名 | 每簇按[合并面板 v2](../resources/annotation/tme_annotation_panel_v2/ALL_PANELS.md)并列给出 V7 的证据报告（`evidence_engine.py`，不选赢家）。名称用 Lodi Minor_2.0 的词汇（[lodi_subtype_vocabulary.tsv](../resources/annotation/tme_annotation_panel_v2/lodi_subtype_vocabulary.tsv)）；Lodi 没有的类型用 V7 的命名。Lodi 发表的21个 DEG 谱作为参照，计算相似度 | Lodi 是人工 RenameIdents，本项目同样人工判定，但证据表全部留档 |
| C1 上皮与成纤维 | 按 Lodi 不做跨癌种细分：非恶性上皮按单个癌种标为器官上皮类型；成纤维细胞、周细胞、平滑肌只保留大类 | 若 CoVarNet 需要成纤维状态，单独提出，属于对 Lodi 设计的扩展，需另行决定 |
| D1 下游阈值 | CoVarNet v1 的样本纳入采用 Lodi 的丰度规则（>500个细胞；谱系内亚型以谱系为分母，T/Mac ≥20，B/DC/EC ≥10） | 与 Lodi 相同 |
| E1 与十月注释对照 | 新注释完成后，在有十月标签的46.3%细胞上做列联表和 ARI/NMI，列出拆分、合并和新增的类型 | 十月标签不参与新注释 |

## 对簇13（约6.8万个非恶性细胞）的处理更正

按 Lodi 的设计，上皮/肿瘤不做跨癌种细分。v1 对全局簇13做的跨癌种子聚类（job 30667416）只保留为诊断。正式流程：
1. 簇13的细胞回到各自癌种，参加 A2/A3 的大类分析。
2. 非恶性上皮用器官条件模块确定组织类型。
3. 簇13中混入的其他谱系、低质量和双细胞，在各癌种的大类分析里识别。

## 尚未完成

- 统一的双细胞旁证分数（A0）。
- 28个癌种的大类分析（A2–A4）。
- 各谱系的两轮分析（B1–B6）。
- 与十月注释的对照（E1）。
- v1 补充模块的逐基因文献出现证据。

## 参数的理由与本项目取值（2026-11-06 补充）

用户要求（2026-11-06）：不要照搬 Lodi 的参数，要理解每个参数的用途。下表左侧是 Lodi 的取值和它要解决的问题（来自代码注释、Methods 原文或方法学常识；属于推断的已注明）；右侧是本项目针对自身数据做的选择。

| 环节 | Lodi 取值 | 目的 | 本项目取值与理由 |
|---|---|---|---|
| 细胞 QC | 基因数 200–6000，UMI >400，MT <25%，按癌种调整 MT | 去掉空液滴和破损细胞；6000上限去掉多重液滴；肝细胞本身 MT 高，所以 HCC 放宽 | 已在 integration v2 统一做过（按样本 MAD 的 MT 上限、Scrublet 或继承的双细胞判定、decontX <0.5）。72个来源由不同的实验室和芯片产生，基因数分布差异大，统一的6000上限在高深度来源会误删单细胞、在低深度来源无效，因此不设。多重液滴交给双细胞旁证和簇级判断 |
| DoubletFinder 预期率 | 固定 3.9% | 对应每样本约5,000个细胞时，按 10x 指南每1,000个回收细胞约0.8%（Methods 中写明） | 按每个样本的实际细胞数计算：预期率 = 0.008 × n/1000，上限0.15。很多来源已经作者过滤过，预期率会偏高，所以分数只用作簇级旁证，不按单细胞删除。DoubletFinder 未安装，先用 Scrublet，参数记录在运行日志中 |
| 回归样本 ID | ScaleData 回归 orig.ident | 去掉同一癌种内样本之间的整体偏移 | **不回归**。在细胞水平回归样本，相当于按样本对每个基因去均值；如果样本之间细胞组成不同（例如某样本 Treg 多），会把组成差异也去掉，而这正是 CoVarNet 要研究的信号。批次改由 Harmony 按来源校正；来源内样本之间的混合程度用单位 iLISI 检查 |
| 回归 UMI、MT | 回归 | 让深度和质量不成为主要成分 | 大类阶段不回归，改为识别并标记由 QC 驱动的低质量簇（Lodi 也做簇级识别）。谱系第二轮如仍出现深度驱动的簇，再加回归，并记录前后对比 |
| 回归 S/G2M、IFN、stress、hypoxia | 回归（Mac/B 另加 hypoxia） | Methods：不回归会出现增殖、IFN 高、解离应激的"亚型"。目标是让亚型对应细胞身份而不是状态 | **采用同样的目的**：在谱系分析中先回归这些状态分数再聚类，同时为每个细胞单独记录状态分数（V7 的 state 轴）。状态作为附加标记，不作为亚型。簇13 的诊断结果说明了必要性：最大的几个子簇由应激、糖酵解、增殖程序驱动 |
| 回归癌种 | 回归 | 去掉以癌种为主轴的分群 | 不单独回归。癌种嵌套在来源内，Harmony 按来源校正已覆盖；组织特异群体按 Lodi 的做法在各癌种内先识别并移出 |
| 批次 | Harmony 按 Technology（3′/5′），B 细胞用 CCA | 他们的数据主要来自同一实验室，最大的技术差异是化学 | Harmony 按来源 θ2（去批次比较 v1：同来源近邻0.154，次优 scVI 为0.212，BBKNN 过度混合）。B/浆细胞若来源 LISI 仍低，再比较 CCA |
| HVG | 2000（第二轮取3000，去掉 HB 后保留2000）；B 细胞去掉 IGHV/IGKV/IGLV | 避免环境 RNA（血红蛋白）和体细胞高突变的 V 基因主导分群 | seurat_v3 HVG 按来源分批计算后合并排名（避免单一来源的特异基因主导）。排除技术基因（MT-、RPS/RPL、MALAT1）、HB*，所有谱系都排除 IG V 基因，T/NK 另排除 TRAV/TRBV/TRGV/TRDV（克隆扩增会按克隆型分群，属于推断，是常见做法）。数量：大类2000，谱系第二轮3000 |
| PC 数 | 按肘部图人工选（10–40） | 保留信号成分，丢弃噪声 | 自动肘点：50个 PC 方差曲线的 Kneedle 拐点，下限15、上限40，取值写入运行记录 |
| 分辨率 | 多个分辨率，人工选主分辨率；稀有群在更高分辨率补 | 不同类型的最适粒度不同 | 计算多个分辨率，用相邻分辨率间的 ARI 和 marker 可解释性选主分辨率；稀有类型用更高分辨率补，并记录每个补标对应的分辨率和簇号 |
| 差异基因 | 大类：only.pos、min.pct .25、logFC .25；谱系：min.pct .10、min.diff.pct .25、logFC .25、return.thresh .01 | min.diff.pct 要求簇内外检出比例差，用来找特异 marker，而不是表达量高的共享基因 | Wilcoxon（每簇最多取2,000个细胞），报告 logFC、pct.1、pct.2、pct 差和 FDR；命名用的 marker 须满足 pct 差 ≥0.25，并且在贡献细胞 ≥50 的来源中分别重复上调（跨来源可重复性，是 Lodi 没有的检查） |
| 低质量簇 | UMI、基因数低，MT 高，表达很少 | 去掉破损细胞 | 同样的判据，相对于同一癌种或谱系中其他簇的分布判断（例如中位基因数低于该对象的第10百分位且没有特异 marker），不用固定阈值 |
| 双细胞簇 | 表达其他谱系 marker，且 DoubletFinder 双细胞率高于预期 | 两个条件同时成立，避免把过渡态当作双细胞 | 同样要求两个条件：其他谱系的 V7 parent 证据，加上该簇双细胞旁证分数高于同一对象其他簇。只满足其中一个时标记为待复核，不删除 |
| 丰度阈值 | 样本 >500个细胞（第10百分位）；T/Mac ≥20，B/DC/EC ≥10 | 比例估计需要足够分母，避免少数细胞造成的离群值 | 按本项目自己的样本细胞数分布取第10百分位，并按谱系设最少细胞数，用于 CoVarNet v1 的样本纳入；先报告不同阈值下保留的样本数，再定值 |

## 簇13（68,374个非恶性细胞）的诊断结果与处理（2026-11-06）

诊断运行（job 30667416，只用于诊断）：簇13共290,546个细胞，在分辨率1下分为74个子簇。

1. **身份清楚的子簇（约1.5–2万个非恶性细胞）**：uncertain 占比低，marker 明确，例如：
   - 乳腺基底上皮（KRT15、SFRP1、GABRP；GSE161529，4,445个）；
   - 肺 AT2（SFTPB、NAPSA，1,480个）；
   - 杯状细胞（MUC2、SPINK4，1,266个）；
   - 前列腺腔面（KLK3、MSMB）、乳腺腔面（ANKRD30A、ESR1）、club（SCGB3A1）；
   - 间皮（ITLN1、C1S；GSE253973）。

   处理：回到各自癌种的大类分析，用器官条件模块命名为非恶性器官上皮或间皮。
2. **程序驱动的跨癌种混合子簇（约4.8万个非恶性细胞）**：
   - 每个子簇约70个来源，uncertain 占63–84%，CNV 相关0.43–0.60；
   - top 差异基因是应激（JUN、FOS、HSPA1A）、糖酵解/管家（GAPDH、RPLP0）、增殖（TOP2A）或缺氧（NDRG1、BNIP3）程序，不是身份 marker；
   - v1 的自动标签（如 plasma_IgG、Kupffer）由簇间相对 z 分产生，是错误的，作废。

   处理：回到各自癌种，在回归状态程序后的大类分析中判定：
   - 若与该癌种的恶性/uncertain 细胞同簇、又没有器官上皮身份证据，标为 `epithelial_tumour_adjacent_unresolved`，建议不进入 CoVarNet v1（待用户决定）；
   - 若显示为其他谱系、低质量或双细胞，按相应判据处理。
3. **双细胞和低质量（约3–4千个）**：例如中性粒和鳞状上皮共表达的子簇、MT 高的子簇。在各癌种大类分析中按判据确认。

## CoVarNet v1 的细胞纳入范围（用户决定，2026-11-06）

本节是 CoVarNet v1 细胞范围的权威说明。决定记录见 [source_decisions_v1.tsv](../inventories/tme_v1/source_decisions_v1.tsv)，计数见 [inventories/annotation_v2/](../inventories/annotation_v2/)。

| 细胞 | CoVarNet v1 | 说明 |
|---|---|---|
| 免疫、基质（成纤维、周细胞/平滑肌）、内皮谱系的非恶性细胞 | **纳入** | 须通过后续的谱系 QC 和节点丰度阈值 |
| 上述谱系中的 uncertain 细胞（67,690个） | **纳入（回收）** | 见下文；保留 `rescued_uncertain_v1` 标记 |
| 非恶性上皮细胞 | **暂不纳入** | 保留注释，不进入 CoVarNet v1 |
| 肝细胞、黑色素细胞、胶质/神经、间皮 | **暂不纳入** | 保留注释，不进入 CoVarNet v1 |
| 上皮及上述谱系中的 uncertain 细胞 | **不纳入** | 不回收 |
| 恶性细胞 | 不纳入 | 不进入 TME |
| 低质量、双细胞、红细胞、未分配的增殖细胞 | 不纳入 | 按判据剔除 |
| GSE179994、GSE182434、GSE189357 的3个 AIS 样本 | 不纳入 | 见 scope_exclusions_v1.tsv |

### uncertain 细胞的回收

大类复核后发现，带有明确免疫、基质或内皮标签的簇中有67,690个 uncertain 细胞，按恶性判定证据分为：

| scATOMIC \| inferCNV | 细胞数 |
|---|---|
| 正常 \| CNV 阳性 | 63,460 |
| 肿瘤 \| CNV 阴性 | 4,143（其中 T/NK 2,864） |
| 未运行 \| CNV 阳性（scATOMIC 失败的样本） | 87 |

最多的是 B/浆细胞（30,411个）：在检查过的簇中，scATOMIC 判为 Plasmablast 或 B Cell，只有 inferCNV 判阳性。解释（未验证）：高表达的 Ig 基因可能造成 CNV 假阳性。

判断依据是簇级的谱系身份，加上两种方法没有同时支持恶性。这些细胞在恶性判定中没有被改写：`malignancy_call` 仍为 uncertain，只在 CoVarNet 纳入层面回收，并按证据类别记录，可以单独做敏感性分析。上皮及其他实质谱系中的 uncertain 细胞（约29万个）不在回收范围内。

### 成纤维细胞的跨癌种亚型

用户决定（2026-11-06）对成纤维细胞做跨癌种亚型。这是对 Lodi 设计的扩展：Lodi 只把成纤维细胞作为大类。基质池（Fibroblast、Mural）因此也进入第二轮，命名使用面板中的 Fib_* 和 Pericyte/Smooth_muscle 模块。骨肉瘤、软组织肉瘤等 CNV 主导癌种中，恶性细胞已被判定移除，但间质样肿瘤细胞可能有残留，需逐簇核查癌种组成。

## 第二轮参数（2026-11-06）

Lodi 在干净对象上把 PC 从18提高到40，是为了分辨更细的亚型。本项目的 Kneedle 肘点在 QC 轮的各个池中都落在下限15，对亚型来说太粗。因此第二轮把下限设为30（上限40），HVG 候选取3,000个（候选集已排除技术、HB、IG V、TCR V 基因），分辨率 0.3、0.5、0.7、1、1.5、2、3，再按相邻 ARI 和 marker 可解释性选主分辨率。

QC 轮的逐簇判定见 [lineage_r1_reviewed.tsv](../inventories/annotation_v2/lineage_r1_reviewed.tsv)。只有其他谱系 marker 而没有双细胞旁证（或相反）的簇记为 `set_aside_review`：不算双细胞，但不进入第二轮。低深度但 MT 和污染都正常、只是没有 marker 的簇保留并标记。
