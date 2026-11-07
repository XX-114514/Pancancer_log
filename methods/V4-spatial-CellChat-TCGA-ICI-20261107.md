# 当前V4重点CM：成像空间CellChat、TCGA免疫及ICI关联

分析对象固定为2026-11-06 V4主v2的14因子、45非上皮节点、164样本/13来源，优先CM03、备选CM06；原注释SHA和W/H/网络不改。完整[配置及原始机会](../inventories/analysis/2026-11-07/v4_spatialCellChat_bulk_v1/README.md)先于本轮空间/临床效应冻结。现有53张Visium、CellCharter和25组CPDB继续，不抢占GPU。

## 成像数据与身份资格

当前研究服务器实际输入为5个Xenium单元（乳腺技术重复2、结直肠/肺/卵巢公开IO各1）及8个CosMx NSCLC归档。CosMx8份包含技术重复，不能称8位独立患者；其他目录空或仅metadata的来源不虚构为可分析表达对象。公开仅来源/标签汇总，个体轴、坐标和所有模型私有。

每个原生面板与既有同版50标签参考相交；参考/留出split保持，分别训练fine50和major的L2 multinomial logistic，C1、balanced、seed2026110703。留出precision和F1均≥0.7、support≥20才允许该标签参与当前面板测量；这是原标签重放诊断，不是独立空间真值。查询细胞至少50原生基因计数、20共同基因检出、预测概率≥0.7且前两类差≥0.2，重点谱系再要求多marker阳性。fine还要求命名marker实测且表达、父类一致；不能靠模型把未测IGHD/XCL1补成阳性。

当前5个Xenium面板实际均缺IGHD及XCL1，因此不能声称精确复现CM03的三个原亚型。可继续支持的Bcell/Endothelial/NK等大类空间背景分析；结果用broad-context命名，不能冒充原CM成员。全部14原CM的精确身份机会单列，不可测与阴性分开。CM06也须同时核对其原B回退标签和IL1B巨噬身份。

## 空间通讯与独立邻近检验

使用现有CellChat2.2.0的`datatype=spatial`，原生全基因面板CP10K/log1p，RNA不做PPI补值。配体、受体复合体及数据库非空调节亚基必须全部实测，主分支仅蛋白信号。小面板使用`identifyOverExpressedGenes(do.DE=FALSE,min.cells=10)`，`truncatedMean/trim0.1`、raw.use TRUE、population.size FALSE、k.min10、100表达置换，seed固定。身份小于10细胞的几何单元不推断，不填成阴性。

Xenium原坐标为微米；指定CosMx legacy global pixel坐标乘0.12028转为微米，禁止将局部FOV坐标直接拼接。不跨原FOV计算邻接，各组织/FOV分别推断；技术重复不扩充患者数。这些输入要求、单位和小面板接口依据[CellChat官方教程](https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.Rmd)及[官方空间FAQ](https://github.com/jinworks/CellChat/blob/main/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.Rmd)，实际函数签名已在安装版核对。

距离全部2×2组合固定为diffusible100/250μm、contact10/20μm，tol5μm。为满足安装版区域距离缩放要求，scale.distance只由几何计算：1.5/同一输入各合格FOV最小正有限非对角区域距离，全FOV/四变体共享；不跨输入直接比较绝对权重。nominal p<0.05候选与每个完整模型假设轴的BH候选并列，所有变体保留。

CellChat的表达置换P不自动检验空间共定位。另固定6近邻且≤100μm的视野内图、999次FOV内标签置换，保持谱系数和几何，报告原14CM对应父类背景的富集及完整14槽BH；单父类/空网络或细胞不足为NA。它只检验支持大类细胞中的局部组成背景，不证明精确亚型、细胞边界接触、空间环境已全面控制或因果。整数标签编码在乘法前显式int64，原int8源码保留、新SHA另记。

## bulk测量先于结局

借鉴[CoVarNet原研究](https://doi.org/10.1038/s41586-025-09053-4)的CMT pseudo-bulk基因程序思路：同560,871细胞的原生CP10K按164样本求均值，再log1p；用既有CMT对其他样本做Seurat5.3.1样本级Wilcoxon，正log2FC>0.25、校正P<0.05、取最高FC10基因，High至少3样本。Seurat5实际FC为组内反变换值之和加1再除组样本数，独立审核遵循已读回源码，不错用“组均值加1”。原网络为空或不足10基因不做CM临床推断，全部原槽仍保留。

本轮临床连续得分为癌种/队列内逐基因患者百分位后，完整10基因等权均值；这是明确的项目扩展，**不是作者AddModuleScore分类的精确复现，也不是原H或45节点组成的精确投影**。不将作者GTEx的0.2阈值移植到TCGA，不优化生存/应答切点。CMT程序可反映整个样本背景，不能保证三个CM成员同时存在；尤其CM07/09/11仍受原先组织/分解诊断限制。

实际充分程序仅CM03/04/07/09/11五个非空因子；CM06/08/10/12无足够基因、CM13仅3、CM14 High只有2样本，原空CM01/02/05单列。CM06本轮bulk不可评估，不能改为一两个marker或放宽门槛凑满。

## TCGA免疫及临床

只读现有Toil/Xena表达、TCGA-CDR和C1–C6免疫亚型。原发实体瘤01材料，按患者合并重复、临床冲突排除；不把不同癌种混成一个免疫/生存P。33癌种×14因子×PFI/OS×基础/分期敏感性共1,848槽；基础模型为程序+pan-T+年龄、可变性别，分期模型另加明确I–IV整体分期（按有序项），覆盖须≥70%。n≥60、事件≥max(20,10p)、满秩且条件数≤1e6；缺失事件/分期不补0。终点可用性/随访按[TCGA-CDR原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC6066282/)保留限制，不因事件达门槛就称所有癌种终点均充分。

免疫部分每癌种/CM固定10槽：C1–C6 omnibus Kruskal，加9个RNA panel/基因（CYT、PTPRC、内皮和PDCD1/CD274/CTLA4/LAG3/TIGIT/HAVCR2）Spearman，共4,620。模块签名与某RNA终点重叠则该槽不作独立免疫证据；至少60患者、足够变异，亚型各组至少10。C1–C6为同一TCGA表达衍生关联背景，依据[Thorsson等原研究](https://pubmed.ncbi.nlm.nih.gov/29628290/)，不当独立测量验证。尚无ABSOLUTE纯度校正，不将RNA相关称纯度独立效应。

## ICI及完整统计家族

预设IMvigor210、Riaz2017、Hugo2016、Jung2019、Rose2021、SU2C-MARK六来源；只用真实临床字段。Riaz/Hugo排on；CR/PR=应答、SD/PD=非应答，NE/未知缺失。Rose档案取样时点未逐人核实，不称治疗前预测。SU2C按作者Keep和精确RNA/患者key连接；公开全TPM实际152 RNA列，不用393全研究人数充当RNA数。

SU2C作者完整公开TPM及临床成员已校CRC/SHA，完整源本机保留。按已冻结签名+免疫控制精确投影88基因/全部152 RNA列，逐值不变，派生SHA与签名SHA绑定后传至服务器；未复用旧34基因子集。慢速完整传输和失败/分块历史保留不用，公开数据来自[作者Zenodo v3](https://zenodo.org/records/11179623)，不申请受控原始数据。

ICI应答84槽（14×6）、OS/PFS168槽，主要程序+pan-T，n≥30、应答少数类≥10、生存事件≥max(10,5p)。各完整家族分别BH/BY，不能将不同家族q当一个总q；不可评估以P1进入校正且展示NA。TCGA/ICI临床总2,100槽；直接Cox Efron/logistic似然和全部拟合设计另独立复算，患者得分/设计绑定、免疫统计及FDR另审核。未经验证的患者重叠、临床/纯度/治疗混杂保留；ICI内关联不等于治疗特异获益或可用预测器。

## 执行与验收

全新`${PROJECT_ROOT}/runs/20261107_v4_spatialCellChat_bulk_v1`，两个有限CPU分支：13成像原生机会逐个准备→邻近背景→四参数空间CellChat→独立数组/计数/BH；并行TCGA/六ICI完整家族→测量/设计绑定及数值复算。原成功准备复用，失败分版本保留；不改共享环境、原注释、GPU队列，不自动公开推送。来源/细胞/患者/坐标行、矩阵、临床设计和模型私有，公开仅方法、通用代码及审阅后的汇总。

实施核对补充：安装版computeRegionDistance的预检分组字段为meta$group，不能只传labels，否则会错误得到0×0矩阵。CellChat_v2保留原基因、距离、置换和细胞，仅修接口；旧不可评估收据不当生物阴性。模型最大迭代警告须单列收敛资格，fine/大类任一未收敛时精确CM审阅不晋级，不将留出重放叫空间真值。

独立表格审核按字面字符串读患者索引，并使用float_precision=round_trip，避免CSV自动类型/浮点解析破坏原键或并列秩。原免疫亚型实际是“Wound Healing (Immune C1)”等六个完整标签；按原字符串精确映射、独立正则核对，原RNA P不变，当前完整4620槽BH/BY按修正后类别P重算。该字段修正不放宽任何组别/样本/基因或P阈值。

在已观察到1000次上限未收敛后，后续同一L2/C1目标增至5000次；同面板缓存只有已收敛才直接复用，否则相同参考/标签/基因顺序下warm-start续优化。该数值求解预算不根据临床或LR效果调参，不覆盖已产出预测；原版本/预算及新SHA保留，精确身份仍须通过收敛、留出重放和原生marker资格。

08:53:57 UTC取得首乳腺大类背景聚合表：133517模型细胞、11类、15完整蛋白互作，主参数273名义候选/原软件BH265。08:57:18 UTC／16:57:18 UTC+8全四变体非零Monte Carlo下限诊断通过：100次重采样原P0经(100P+1)/101和全1815假设BH后各0候选，最小q0.0593–0.0678。原数组/P未改、无细胞替换/模型重跑，不将候选或软件q0称稳健通讯。通用诊断R、聚合表和SHA v3保留；此诊断不是空间零模型或完整tie校准，精确CM03仍未验证。
