# V7 整套 panel 修复与文献引导的注释流程

本轮交付的是完整规则修订、可执行证据流程及验证，尚未完成全V7逐细胞重新注释。权威资源为[registry](../resources/annotation/v7_panel_repair_v1/registry.json)，[全部112个panel显式表](../resources/annotation/v7_panel_repair_v1/ALL_PANELS.md)与[162项迁移表](../resources/annotation/v7_panel_repair_v1/old162_to_new.tsv)覆盖旧配置。旧字段与历史结果保留。

## 为什么必须改整套规则

旧panel同时混合了父身份、细胞亚型、功能状态和所谓污染。在部分规则中，只对通过阈值的基因取均值、先选择最高分再判断是否合法；共享细胞毒程序、TIGIT等可能遮挡真正需要检查的分支。此次取消跨角色的赢家排名和统一margin，不是简单调低0.015或改一个NK名称。

| 旧问题 | 新处理 | 证据边界 |
|---|---|---|
| T的NK_contamination含GNLY/NKG7/KLRD1 | Cytotoxic与T_NK_boundary，T/CD8/NK身份各自检查 | 共享效应不证明NK污染，T+NK程序不自动叫NKT |
| CD8_exhausted可仅由TIGIT获高分 | CD8A/CD8B分支与PDCD1/HAVCR2/TOX/LAG3状态拆开 | TIGIT/CTLA4也见Treg；耗竭功能不能由RNA直接证明 |
| Treg被其他高分panel挡住 | FOXP3与IL2RA/IKZF2/IKZF4组件，另保留Chu完整14基因功能签名 | 全部候选并列；不能把旧T6整体改成Treg |
| IFN、cycling占据亚型位置 | 合并为跨谱系独立状态 | 不抹去其生物学变化，也不回归后再声称不存在 |
| 43项contamination触发冲突 | 全部退出亚型排名，转邻近谱系/边界审阅 | 不能据此删除细胞或确诊doublet |
| B/ASC中CD79A、CD37或JCHAIN被当独占marker | B、ASC组合及B_ASC_transition；Ig同种型另轴 | 单个XBP1/JCHAIN不决定浆细胞，分泌/增殖不决定成熟度 |
| adaptive NK正负方向混合 | KLRC2候选；FCER1G/SYK/SH2D1B/ZBTB16降低证据另列 | 必须同NK背景对照，不要求所有adaptive NK都降低 |
| Mast Activated_IL1RL1由HDC/LTC4S替代 | Mast身份与IL1RL1/CMA1/CPA3表达描述分开 | 未测IL1RL1不得以该基因命名；HDC不证明激活 |
| Pericyte/SMC/Mesothelial混在Fibroblast亚型 | 独立入口与比较模块 | ACTA2不独立证明myCAF；CAF需肿瘤来源背景 |
| PMN-MDSC、M2、HSC/MPP、GMP/MEP等功能/谱系断言 | marker相关方向或祖细胞priming | RNA不能独自证明抑制功能、多潜能或正常造血身份 |
| 泛癌组织上皮共享名称过宽 | 器官条件、专属组合及恶性独立轴 | KRT19不等于所有器官luminal；melanocytic不归上皮 |

## 文献中借鉴什么，哪些属于本项目实现

| 来源 | 已核实逻辑 | 本轮采用及限制 |
|---|---|---|
| [Lodi作者代码](https://github.com/lambrechtslab/Pan-cancer_Lodietal) | 谱系内检查、人工剔除QC/doublet群后重建；HVG/PCA/图/多分辨率；RNA FindAllMarkers；人工RenameIdents | 采用流程顺序与DE参数轮廓；不复制未知Technology或回归变量。其MarkersToAdd占位不等于公开了完整自动panel |
| [Qian 2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) | 泛癌共享/组织特异类型；按组合marker检查并召回邻近父类中DC | 父群不是不可推翻真值；DC与T/NK/ILC边界并行检查。DEG摘录不冒称作者完整分类器 |
| [Chu 2023](https://www.nature.com/articles/s41591-023-02371-y)及[作者代码](https://github.com/Coolgenome/TCM) | T细胞分支、marker与独立功能签名 | 原样保留Table S6的14基因Treg功能签名；Treg身份组件为项目审阅组合，不是论文的通用二元阈值 |
| [Ma 2025 DC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) | cDC1/cDC2/pDC组合、稀有DC与成熟迁移程序 | cDC1用CLEC9A/XCR1/CADM1方向，cDC2用CD1C/FCER1A/CLEC10A；LAMP3/CCR7/FSCN1单列成熟状态 |
| [Rao 2017](https://www.nature.com/articles/nature20810) | RA中PD-1高、CXCR5低的CD4 helper群及功能验证 | 在TME只作为Tph-like候选逻辑；不由PDCD1判耗竭或由CXCL13判Tfh/Tph |
| [Aizarani 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6687507/) | 肝细胞与胆管方向在EPCAM群内亦不同 | ALB/ASGR1/HNF4A方向与KRT19/CFTR方向；仍需器官背景 |
| [Tirosh 2016](https://shaleklab.com/wp-content/uploads/2024/03/Tirosh-et-al-Science-2016.pdf) | melanoma的MITF与AXL相关状态 | TYR/PMEL/MLANA支持分化程序，低表达不能排除去分化来源；不判良恶性 |
| [Henry 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6411034/) | 正常人前列腺细分后用DEG和蛋白/组织验证 | KLK3/KLK2、ACP3/MSMB方向；不能覆盖全部肿瘤可塑性 |
| [UCell](https://doi.org/10.1016/j.csbj.2021.06.043)及[官方实现](https://www.bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_vignette_basic.html) | 单细胞内秩评分 | 可选辅助；全来源可测基因宇宙内计算、记录版本和maxRank。没有统一接受阈值，不称细胞类型概率 |

其余NK、ILC、Mast、成纤维、HSPC、肺及MDSC命名依据见[来源表及范围](../resources/annotation/v7_panel_repair_v1/references.tsv)。全部634条模块—基因记录均给出[逐基因依据状态](../resources/annotation/v7_panel_repair_v1/gene_reference_trace.tsv)：266条有本轮已核验作者材料出现证据，368条仍待进一步逐基因核实。该比例不是正确率；也不能从出现证据推导作者验证过项目组合。

## 新计算流程

1. **范围与测量先确认。** 使用全RNA原始计数、细胞元数据、基因×来源可测矩阵。旧父标签用于寻找入口，同时召回邻近谱系；不可测不补成阴性。仅显式版本化alias映射，别名碰撞必须先按原始feature ID解决。
2. **谱系内QC。** 初轮保留全部候选，检查深度、来源、双谱系同细胞RNA、doublet工具旁证与环境RNA。只有人工审阅的逐细胞keep/reason清单才能进入下一轮过滤；本程序不冒充自动doublet识别器。剔除后从counts重建，旧全局父标签不直接继承为终判。
3. **重新构建分群。** 全文库CP10K/log1p归一化；以来源覆盖率筛模型候选，重新选HVG、PCA、图和多个分辨率。群数量由图聚类产生，没有目标数量。若技术批次确实已知，可启用Harmony；默认不自动回归癌种、IFN或应激。分辨率一致只是诊断，还需要来源/患者稳定性与marker可解释性。
4. **RNA差异证据。** 使用Seurat RNA data/Wilcoxon，`min.pct=.10,min.diff.pct=.25,logfc.threshold=.25,return.thresh=.01,only.pos=FALSE`。这些为Lodi式DE输出参数，不是身份阈值。汇总DE只测试所有来源均可测基因，其余在来源内测试；归一化分母始终保留全文库。来源内DE不等于患者独立复现。
5. **身份、状态、同种型和边界分别解释。** 所有panel都输出，不选单赢家。每个核心组件显示可测基因和检出比例；`min_component_max_detection`只是边际检出描述，不是概率、同细胞比例或最终接受分数。命名marker必须本身可测；检出与有差异富集证据分别记录。
6. **同细胞证据独立计算。** raw≥1和raw≥2 UMI两档检查各组件至少一个基因阳性；为保证完整组件可评估，要求该来源测到该组件全部备选基因。保守分母可能小，但不可把未通过者算身份阴性或把阳性证据数量等同类型总数。这是项目诊断，不声称文献规定了这两个阈值。
7. **可选UCell。** 仅调用官方R包；每个来源完整测到某个签名时才评分，缺基因另列；保留完整来源测量宇宙而不是marker子矩阵。不同来源基因宇宙仍可能不同，不能直接对其分数设共同阈值。降低方向不靠零表达自动判定，须人工选择同谱系对照。本轮R环境未安装UCell，因此该可选分支未运行验证。
8. **终判与原因拆开。** 记录父身份、候选分支、状态、命名marker、QC及人工决定；大类已支持但亚型未解决应明确写出。建议原因字段区分measurement_gap、insufficient_branch_evidence、shared_program_only、mixed_lineage_review、QC_review、reference_scope_gap。证据不足时允许停在大类或“大类+marker相关群”，不强制HSPC/ILC/红系亚型。

补充参数边界：证据引擎在已提供的DE表上额外使用`p_val_adj<.01`作为项目保守的正向报告条件；它不等于FindAllMarkers的`return.thresh`，不冒称Lodi作者相同过滤，也不决定细胞身份。Seurat提供的校正值保留原名，不称为BH/FDR。模型候选覆盖率小于1时仍有结构性缺测影响，真实运行需要与共同可测基因模型作敏感性比较；本轮合成测试使用覆盖率1。

CellTypist/scATOMIC作为独立参考输出，与上述表达证据并列。新版脚本不重新运行它们，不能把既有模型多数票当训练真值或最终注释。

## 已执行的验证与实际限制

- Python 15项反例测试通过：包含T/NK共享marker、TIGIT遮挡、命名marker缺失、signed adaptive NK、alias碰撞、DEG别名、错误分母和边际/同细胞区别。
- 服务器Seurat5.3.1合成测试：360细胞、2100feature、2来源；按测试QC清单剔除4细胞后，356细胞完整重建出3群。fixture/workflow/independent verify均真实exit0。1227条表达记录、6条Treg同细胞记录和来源缺测分母独立核对通过；CP10K最大绝对误差约4.69e-13。该测试只证明程序行为，不是生物学准确率。
- [真实聚合回放](../inventories/analysis/2026-11-05/v7_panel_repair_v1/replay_v2/manifest.json)：旧40群1073743细胞范围→4480条panel证据；T6最新12群156631细胞范围→1344条panel证据。两范围重叠，不能加总。新panel缺的基因不伪造；T6 top10 DEG明确是截断输入，未列出不等于不富集。
- T6 K-02的FOXP3/IL2RA检出约80.6%/71.2%；K-03约8.1%/8.2%，K-04约10.3%/13.8%，后两群CXCL13约48.5%/50.1%。这些已有观测要求分别审阅，不能将旧T6全部归Treg。其既有全矩阵分析与本轮聚合回放分开记账。

全V7的409基因重新提取、各谱系真实全矩阵重建、来源留出复现及人工最终命名尚未完成。Oligodendrocyte、Neuron和近端小管方向仍标记参考适用性缺口。当前适合开始逐谱系修订，不应把资源和测试通过宣传成已经获得发表级全图谱标签。

可执行代码与输入合同见[脚本说明](../scripts/v7_panel_repair_v1/README.md)；实际运行、哈希与退出码见[run](../runs/20261105_v7_panel_repair_v1.md)。
