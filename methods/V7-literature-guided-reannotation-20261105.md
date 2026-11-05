# V7逐类注释复查：文献做法与执行合同

本轮输入是冻结V7的4676787细胞，原标签和CNV结论只读保留；新增证据与修订另存。完整注释指每类有可解释的身份层级、状态层、证据和处理结论，不预设必须消灭所有coarse/unresolved。

## 文献与可复现边界

- [Lodi 2025作者代码](https://github.com/lambrechtslab/Pan-cancer_Lodietal)：初始样本QC/DoubletFinder之后得到大类，在同谱系集合中再查看QC、竞争谱系表达和已有双细胞证据；清理后重新做谱系HVG/PCA/批次处理/图/多分辨率聚类，以RNA marker和手工核对命名。T脚本将T与NK一起复查。
- [Qian 2020正文](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf)：逐癌种和谱系复查、marker支撑与适当分辨率选择；共享与组织特异群分别评价，不能以批次混合最好作为唯一目标。
- [Chu 2023作者代码](https://github.com/Coolgenome/TCM)：T细胞CD4/CD8/非常规T/增殖分支，类型与功能程序分开；参考投射作为旁证。

以上决定分析顺序和证据结构。Lodi公开代码有未定义regression_conditions、MarkersToAdd占位和版本依赖；不能宣称逐参数复现未知部分，也不将原论文某个cluster编号、分辨率或预期双细胞比例直接套到V7。实际参数与适配要写入每次run。

## 第一步：冻结与原因归因

同时读取V7 major/final/confidence、strict eligibility、requires_recluster、原cluster、父类修订原因和原始计分/扩展panel/描述性marker表。建立每个最终标签及其原cluster的crosswalk，计数守恒。历史原因与当前新证据分开存储。

| 情况 | 处理 |
|---|---|
| rare_lineage_no_forced_subtype | 主动保留粗类；复查身份后可作为终点 |
| major_refined_subtype_requires_recluster | 新归入父类的细胞未在对应谱系重新建模；纳入新谱系分析 |
| coarse_provisional或strict_major为false | 父身份仍需复查，不能只去掉标签后缀 |
| expanded_signal_without_canonical_anchor | 旧相对表达/亚型panel条件未满足，不等同父类marker不表达 |
| top_expanded_panel_is_contamination_conflict | 查看同细胞共表达、父类关系、参考与QC；不自动定为doublet |
| marker-only或no_robust_de_novo_marker | 描述性命名/无稳定谱系内差异，复查分辨率、状态和来源 |

## 第二步：父身份、QC和测量范围

1. 先查父类marker表达比例和完整panel；共享程序或单一异源marker不直接判双细胞。
2. 按来源基因可测性计算分母。未测基因不是阴性；不把共同10432基因作为后续HVG/marker默认上限。全基因对象须精确验证与V7细胞轴一致，只读取counts及证据，不继承后加标签为答案。
3. 按来源/样本和谱系看nCount、nFeature、mt、已有doublet/ambient指标，保留原QC记录。新的排除记录原因及敏感性；稀有谱系不能按全局阈值机械剔除。
4. CellTypist保留模型、范围、置信度及来源重叠审计，作为辅助。免疫模型不能否决非免疫身份；稀有类不因参考缺失而被否定。

## 第三步：谱系内重新建模与命名

有父类调整、混群或亚型分辨不足的集合，使用全基因counts重选HVG、PCA、邻接图，检查多分辨率marker/QC/来源结构。批次变量由真实技术/来源元数据决定，记录未校正和校正敏感性；不直接回归癌种、IFN或研究目标状态。

命名同时依赖父身份、正marker、竞争marker、谱系内DEGs、跨样本/来源复现及分辨率稳定性。增殖、IFN、应激、细胞毒单列状态。小群单独检查，不靠增加分辨率强造亚型；旧阈值命中不等于发表级验收。

红系、血小板、Mast和HSPC可在身份明确时保留大类或大类+稳定marker群；HSPC若涉及肿瘤/髓系邻近程序继续审阅。ILC先核对TCR/CD3阴性证据的可测性、ILC阳性程序和NK/肥大/髓系竞争信号，不能默认其父类已成立。

## 验收与交付

每一谱系交付旧标签→原因→证据→修订标签表、marker表达/DEG、QC和来源分布、必要嵌入图及分辨率对照。细胞级overlay留服务器；公共仓库只存审阅后汇总、通用代码和脱敏日志。原标签不覆盖，终端manifest及SHA后才能称运行完成；生物学验收另列。

## Marker图的固定展示要求（用户2026-11-05补充）

图中必须同时展示marker和对应的身份/亚型/状态：列上方按panel分组写明目标名称，底部保留基因名；行标签给出cluster编号、工作性身份或状态、review标志及细胞数。不能只给一排marker或只给数字cluster。亚型尚未确认时写candidate，状态与父身份分别标示；共享/竞争谱系marker明确写对应谱系。

Mast的MC_T/MC_TC相关RNA panel参考[JCI2025原研究](https://www.jci.org/articles/view/174981)，仅为候选表型关联，不能凭CMA1未检出定义MC_T。[跨器官Mast原研究](https://rupress.org/jem/article/220/10/e20230570/214245)提示人类转录状态不能简单按传统二分法硬切。另以[泛癌pan-MC作者Fig1代码](https://github.com/xiaoyanger01/pan-MC/blob/main/Fig.%201.R)的10个三基因panel作为参考图；Ref C0–C9不等于本项目cluster编号，也不直接迁移作者标签。panel到目标的映射和来源单独输出表。

图分为身份/候选亚型/状态、竞争谱系QC、文献参考panel三张，保留逐来源可测性分母。颜色为各基因内部标准化的log1p平均CP10K，点大小为可测细胞中检出比例；来源覆盖不足80%用×标识，不把未测量当0。

Mast本轮继续使用原始样本与doublet_method分层的完整背景校准doublet分位；高分与多谱系共表达的联合条件只产生doublet候选。H敏感性采用pan-MC代码中的counts/feature/mt数值QC、跨功能Mast阳性marker及无强竞争程序，并以49个实际Mast来源中至少80%可测的基因作为模型候选。完整RNA assay和marker轴保留；按GSE做Harmony是明确的来源敏感性，不能冒称已复现未知Technology变量，也不能以混合更好作为生物学准确率。没有任何整群自动删除。

## Mast扩展复核与层级输出

H的初始竞争panel只有T、髓系、间质和上皮。来源敏感性RNA DEGs仍发现B/Ig程序，所以I补充B、浆细胞、内皮、血小板并扩展T/髓系panel。Ig单独不作为浆细胞身份；ITGA2B单独不作为血小板身份；共享TYROBP不作为异源髓系冲突。8类panel出现至少3个基因各≥1 raw UMI，或至少2个基因各≥2 raw UMI，触发复核。阈值是保守、可追溯的审阅规则，未训练成身份概率或双细胞准确率模型。

Mast父类支持要求TPSAB1/TPSB2至少一个，加CPA3/KIT/MS4A2至少一个，各≥2 raw UMI；低于规则记录low signal，不据此否定真实Mast。竞争程序与原记录样本/方法内doublet分位≥95%联合时才记joint doublet candidate，仍不确诊双细胞或自动排除。记录样本单位不默认独立患者或经验证的独立建库。

细胞级新overlay将`parent_identity`、`parent_support`、`parent_conflicts`、`identity_subtype`、`subtype_status`、`state_programs`、`qc_review`与`admixture_hypothesis`分开，保留全部旧字段。对应[panel语义审计](V7-panel-semantics-audit-20261105.md)的身份/状态/QC区分，`resolved_at_major`仅表示此轮规则下可停在Mast大类，不等于最终主分析纳入、谱系纯度或发表级验收。MC_T/MC_TC阳性RNA关联非互斥，仍是候选；未给任何细胞硬赋成熟亚型。

[当前Mast证据与带目标名称的新版图](../inventories/analysis/2026-11-05/v7_literature_guided_v1/Mast_source_review_v1/README.md)包括全部45036细胞分流和36007细胞敏感性模型的15群。模型子集和全体复核的分母不同，群富集状态与单细胞双marker阳性数量也不同，不能互换。

## Tcell_6的工作注释与阳性证据

原T6全156631细胞单独重建只是诊断范围，不替代全部T/NK/ILC共同父类分析。Treg核心与共享checkpoint、CD4/CD8背景、Tfh/Tph候选和状态分别画图。完整Chu2023 Table S6 Treg功能签名作参考，其余组合注明本项目适配，不逐字冒称作者判定器。FOXP3加IL2RA/IKZF2/IKZF4的同细胞阳性、T身份支持、RNA DEGs和来源结构共同用于review，不凭TIGIT或自动模型多数票转Treg。

阳性raw规则只给证据下界：严格T组件要求至少一个CD3D/E/G和至少一个TRAC/TRBC1/2各≥2 UMI，再查看FOXP3及至少一个第二核心marker。未测组件与低于正证据阈值分开记录，未通过不作Treg阴性，也不抹除其他T身份证据。功能/状态可以与多个身份共存；CD8耗竭必须先有CD8身份背景。跨分辨率最大重叠Jaccard仅为同图敏感性，不等同独立来源留出复现。

K-00/01不为增加亚型数而强分resting/effector；K-03/04先用CD4+marker描述，保留Tfh/Tph候选；K-06/08/09保留独立混合程序/QC复核。K-11虽跨分辨率不变，但只有2细胞、无有效DE检验且原记录计数异常高，不能因稳定而命名新亚型。[完整表和带目标标题的四张图](../inventories/analysis/2026-11-05/v7_literature_guided_v1/T6_diagnostic_v1/README.md)。
