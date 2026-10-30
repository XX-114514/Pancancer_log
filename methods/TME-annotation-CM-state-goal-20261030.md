# 新TME：自动注释、CoVarNet比较与状态模块实施计划

2026-10-30。用户要求先规划、随后作为持续goal实施。范围为当前来源资格通过的新TME；历史泛癌/免疫CM保留其日期与版本。本计划不预设新CM数目，不以恢复旧CM02或得到更小临床P值选择注释、rank或状态。

## G0｜固定对象和分析单位

复用当前接受的原生行、采用RNA层、来源资格v7和标准图。模型拟合细胞与资格通过但未拟合的细胞分别登记；未拟合者不自动获得新cluster。来源资格、细胞QC、身份及CM输入资格分别验收。

先建立私有细胞—来源资格工作单位—生物材料—物理文库—患者桥。沿原文/已审元数据判定合并技术重复；不同组织、区域或治疗阶段保留各自材料。患者缺失时按明确的材料/文库替代单位报告，不把工作单位数称独立患者数。不能从名称模式或细胞组成猜分选、患者或材料。

主发现使用合格、未经谱系富集的肿瘤组织scRNA材料；原发与转移分别组织，正常/癌旁、血液、产品、富集材料和snRNA保持各自用途。来源未知/分选未知不进入主发现资格；来源与癌种重合时作为可观察情境说明。单癌种或单来源模块仍可有价值，另列适用范围和证据等级。

产物：完整行轴/SHA合同、分析单位桥、QC/材料/测量mask、分来源覆盖与排除原因。仅发现集计算选项的完成不代表独立临床验证。

QC继承[既有历史验收](TME-standard-reintegration-20261029.md)：原5M的逐样本QC/Scrublet、新来源及补充分支的doublet/环境RNA验收与采用层分别恢复。carried字段缺失不等于没有做过QC；本轮不为字段覆盖率重复过滤或重估doublet/环境RNA。准备标志、历史处理完成、实际采用及当前用途资格分开。

## G1｜自动身份注释与多方案

1. **CellTypist逐细胞方案**：完成已启动的官方预测及独立审核；分别保留原始标签、粗身份、细身份、状态提示、模型适用范围及不可测原因。高/低粒度模型是同一工具的粒度比较，不当作两种独立算法。
2. **SingleR独立参考方案**：核实可用现代版本与缓存公共参考。先在本群全部成员的采用RNA汇总上运行官方聚合注释，保留pruned/未定标签；同层级比较CellTypist聚合与SingleR聚合，避免把工具差异完全混同于逐细胞/群级差异。群均值不能揭示混合群的内部组成，争议群继续用逐细胞结果和本群marker审阅。逐细胞扩展按资源与必要性另登记。
3. **cluster-marker与自动共识方案**：以本群多基因模式和适用参考形成可追溯大类；保留数值分群及精确成员SHA。细身份不足时回退大类，不借单个展示marker断言功能、恶性或细状态。共识不是独立真值。

建立共同身份词表及成员crosswalk。原始自动标签中的cycling、activation等状态词另存，先生成身份归并视图，再在G4正式测量状态。没有合适参考的非免疫群保留有来源标记的marker大类/未定，不把免疫模型输出冒充全TME参考注释。

**两个比较轴分别固定：**

- 注释引擎比较：在相同细胞/聚合单位、共同可测基因和可比身份层级下，比较实际改变归属/合并的方案。
- 粒度比较：同一引擎下比较粗细身份及保留cluster的方案；标准图不同分辨率或Harmony/scVI属于另一项分群敏感性。

只换展示名而成员、父类、分母均不变的方案，数值输入及CM应完全一致，作为负对照。相同名称后缀不同不直接合并；相同规范身份可以另设合并视图并记录改变的成员集合。

G1在查看CM结果前冻结映射和不确定性处理。主要比较使用共同可评估细胞/材料集，同时报告各方法完整可覆盖域及未定比例，避免只保留共识细胞掩盖覆盖差异。

执行补充：匹配群均值CellTypist仅作同输入的工具诊断。其标签若与同来源同群的逐细胞投票不同，不将群均值细标签广播给全部成员；以完整逐细胞结果及相应群投票建立主视图。SingleR保留完整来源查询批次，独立复算反转查询列顺序以检查标签、pruned结果和delta的轴不变性；不能抽几个查询群后假定其自适应pruning阈值应与完整批次相同。

身份方案分为原始参考表型、身份归并、身份×数字群、固定群CellTypist投票、固定群SingleR及marker/自动共识。原始参考表型可能包含状态词，按其真实含义报告，不改称纯身份。参考主视图保留适用父类内的官方best-match；不将只在prob-match模式中有意义的0.5阈值强加到best-match。唯一正提议另作置信处理敏感性，不解释为校准后的生物学准确率。父类冲突、范围外或细身份未定保留明确回退节点；不从完整父类分母删除。数字群定义的粗身份子群与无细身份的笼统回退节点分开，前者仍可提供可检验的组成区分。

参考词表在CM输入前以完整精确字典解释，禁止运行时用子串猜大类。例如Mast不因字母结尾与T-cell文本相交而改为T；cytotoxic、naive、cycling等原参考提示与规范身份/正式状态评分分开。修正版语义层不覆盖原始模型标签或概率。固定群CellTypist投票与SingleR群级参考的实用流程比较仍包含参考和聚合路径差异，不能把其差异全归因算法，也不据此声称身份准确率优越。

性能实施：函数级计时定位到object基因名字数组的重复成员匹配。无损转换名字数组表示后，仍调用同一官方Classifier构造和未修改的celltype；不修改表达、名称/顺序、模型数值、缩放、scope或投票规则，也不猴补丁全局NumPy/安装包。先核完整名称membership及普通/新路径的标签、decision和probability逐位相同，接受后记录确切owned PID和完整来源检查点再受控切换。只复用完整来源，旧产物及切换原因保留，后继全轴审核不省略。[NumPy成员匹配](https://numpy.org/doc/1.26/reference/generated/numpy.isin.html)与[官方CellTypist推理](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html)为接口依据；性能增益由本项目实测，不当生物准确性结果。

覆盖补充：原项目要求模型特征覆盖至少80%，该门槛不是官方分类器的推理前提。对事前登记且RNA可测、已有大量模型匹配基因的受影响来源，保留原严格视图，另按官方匹配特征接口补全主参考预测。原有预测、模型数值与表达层不改；基因覆盖率和严格范围敏感性继续报告。普通路径与加速路径须在相同输入及相同批大小比较；不同批大小间的浮点末位差异不能当模型改变。生物学准确性不由数值等价推出。

首轮实际冻结五个共同区室/相同数字图视图：共同区室×cluster对照、CellTypist原参考表型×cluster、规范身份×cluster、SingleR HPCA群提议×cluster和SingleR Monaco免疫群提议×cluster。非免疫部分共享数字cluster分区；器官专用参考另设适用测量面板，不在未应用来源填零。免疫参考父类冲突/未定留显式回退并保留完整父分母，该回退不作为具名细亚群发现节点。大类内可解释但未细分的规范身份与参考冲突分开。SingleR仍是来源×群的实际参考提议，不把广播标签说成逐细胞推断或纯工具准确率实验。此首轮是完整工作单位的探索发现，后续来源/患者块留出使用独立冻结的重拟合和投影。

## G2｜组成矩阵与CoVarNet

采用作者CoVarNet的子群/所属大类内比例、minmax、nsNMF和相关网络基线。输入先显式处理测量mask，再调用作者接口；不调用会将缺失大类NA置0的默认频率转换。大类充分捕获但子群未见可记0；大类未捕获或身份/材料未审记NA。

初轮沿用已有项目可解释的门槛作为覆盖预案：材料至少100细胞，父类至少20；候选节点总计至少100细胞且至少4材料观察到，比例有变异；拟合面板至少12材料/8分析块。它们是本项目门槛，非原论文统一标准。完成覆盖审计后、在CM拟合前冻结实际面板和分母；主比较使用各方案相同材料掩码。

先按覆盖构建明确的完整面板及谱系内面板，不要求稀有大类在所有材料出现，也不将缺失补零凑成大矩阵。仅有一个子群的大类，其大类内比例恒为1，不能用于该范围的共变发现；保留目录并报告不可识别。全TME总比例是单独敏感性口径，不能混称作者的大类内比例。

参数基线沿既有项目实现：rank从2至min(12,节点数−1,分析块数/4)，初轮每rank10次、最终30次nsNMF，保留训练侧cophenetic/RSS与是否触及搜索边界。选择规则属于工程启发式，另报共同rank和跨rank匹配。相同输入/rank/种子比较与各方案自选rank比较分开。缩放、节点筛选和rank选择只在发现/训练部分进行。

主发现同时交付总体组成主题和情境内边；来源/癌种/阶段不是必须全部回归掉的“坏信号”。报告来源等权敏感性及来源内方向、同癌种多来源/多癌种支持。不同癌种均值差不能单独证明癌种内个体协同。重复患者按块重采样；不同阶段不先平均成一个无时间标签患者。

组成输入器要求已接受的成员代码、原生轴和用途绑定，先取各方案共同可评估细胞域，再按预先冻结的父类/材料面板取相同观察单位。细身份未定使用明确回退节点；删去发现节点时不删其父类分母。完整面板在任一方案覆盖不足时共同标记不可比较，不自动为该方案更换样本。输入准备与独立审核分开，准备清单不能直接进入正式执行器。候选患者块只是重复取样敏感性依据，未据此新增独立患者认证。

作者代码核查发现，`pair_correlation`调用的当前`psych::corr.test`默认先作Holm调整，作者随后再对取出的上三角P值作BH。兼容基线保留该真实行为并明确命名；另输出由Pearson相关及观察数独立计算的原始P/BH敏感性网络，不把作者`pval_fdr`称为单次原始P的BH。本轮统计复算不以更换校正获得更小P值为目标。合成API验收已包含真实作者nsNMF、非均匀块权重的KL恒等式、有效WS重构、展示重命名及空网络序列化；它不代表人类TME发现。

节点×工作单位测量机会另存显式布尔表，先在实际可测单位选候选节点，再冻结所有方案共同的测量面板。未应用分类器的零输出不能当阴性频率。手工反例分别验证真实零保留、未测零触发所有方案共同剔除该单位。实际输入另以不同计数实现从细胞成员重算全部分子与父分母，并直接复算正值和真实零；验收前不启动人类数据NMF。

未解决单位、QC或必要测量门槛的范围只交付覆盖报告/明确的探索分支，不静默提高为正式CM输入。

执行顺序补充：首T/NK–髓系五方案先运行各自rank调查；另预定共同rank=4的五方案对照，再运行其余三个已验收父区室面板。共同rank属于事前可识别的计算对照，不称最佳生物学模块数。扩展队列最多并行两fit，依赖首轮真实数值验收；每次科学失败停止队列且不自动重跑，保持单项请求、进程、实际退出及原始产物。此有限科学队列不是无限轮询或按效应追加候选。

## G3｜比较CM而非比较名字

交付各方案节点/分母、W/S/H、完整相关/网络及覆盖。nsNMF按有效基WS对齐，同时保留raw W；不把表达cNMF和组成nsNMF混同。

比较内容包括：

- 共同父身份、共同cluster锚及精确细胞成员crosswalk上的载荷主题和成员重叠；成员拆分/合并单列。
- 一一匹配并允许未匹配，避免多个旧CM同时挑同一个新CM；软权重余弦、加权重叠、网络边重叠和相同材料活动相关分开报告。
- 多种子、患者/材料块bootstrap、来源等权、留一来源、跨rank及同细胞预算稀释稳定性。
- 来源限定、多癌种一致、同癌种多来源支持及无法比较的范围；不因来源依赖自动否定模块。

新版本使用日期/范围/方案标识，模块描述名在稳定成员审核后增加。旧CM02仅作为冻结历史主题的对应分析，不作为选择目标。发现患者/来源留出基于冻结图谱条件，不能称完全重新整合的端到端外部验证。

实际对比在同一训练材料、共同细胞域内建立五方案成员签名的共同细分原子，再按成员列归一化映射各拟合节点。保存有效WS与raw W两个口径、归一化特征载荷余弦/加权Jaccard和相同样本的H活动Spearman。匹配用一对一Hungarian及显式未匹配，首阈值0.6是项目描述门槛，不当生物验证。归一化特征载荷不是细胞数量；训练成员映射和种子一致性不替代来源内检验、来源留出及共同rank对照。

## G4｜身份之上的状态与模块

复用并核对现有[预定义状态程序接口](predefined-state-programmes-and-module-interface-20261030.md)。以固定文献字典和官方UCell为首轮坐标，T细胞可加入适用TCAT/starCAT参考；cNMF/Spectra用于经训练侧审核的发现补充。先谱系内解释，保留处理/解离/深度相关程序的诊断身份。

状态是多基因程序向量，同一细胞可多轴活跃。先审计原生基因覆盖、别名、父身份适用性、来源内拆半/稀释可靠性及检测深度影响。类型内均值/中位数优先；尾部分位和活跃比例须通过精度检查。字典和阈值在模块发现/临床评价前冻结，未测为NA。

资源准备复用已冻结字典与其Cytopus审校侧车，在当前全部来源的已接受原生基因轴上登记逐程序覆盖，同时保留共同基因空间口径。已有T/NK评分不复制为全TME状态。可测覆盖不等于适用谱系、原始文献审核、功能验证或模块准入；原生统计阶段的验收与后续Wilcoxon阶段终态分别绑定。

T/NK、B/Plasma等共享组成区室不能替代状态评分需要的生物身份适用域；状态分析先分辨实际T、NK、B及Plasma等可支持身份。高状态比例的阈值在训练侧固定，不能每样本各自强制取前10%而制造恒定比例。多轴活跃仍保留连续评分；若需CoVarNet输入，采用明确的身份内计算分区，不能把重叠状态阳性集合冒充互斥子群并重复计入分母。

两条状态模块路线分别定义：

测量候选进一步要求至少80%的20次拆半可估计，且至少4个合格单位观察到正程序分数；不能只从少数可估计拆半取高相关来宣称可靠，也不能由单个正单位驱动主特征。此处正分数表示实际测得程序表达排序信号，不等于功能被激活。多来源支持与单来源支持分列，后者仍保留相应情境价值；这一步不是最终状态特征冻结，不替代RNA深度、背景基因空间、程序重叠与源/单位稳定性审阅。

执行接口补充：采用冻结Spectra/Cytopus作者资源，官方pyUCell 0.8.0在每个来源的完整、唯一可测原生基因排序空间评分，不以HVG空间代替。采用RNA层、物理列RNA总量和原生行轴逐块核对；maxRank固定1500、不平滑，基因缺失先通过本项目至少5基因/70%门槛，再按官方skip策略评分。不可测程序保留NaN。每来源以独立正表达秩方程复算，并检查分块不变性；旧T/NK分数不复制为新全TME结果。[UCell作者接口与秩公式](https://github.com/carmonalab/UCell)支持评分方法；5/70%是本项目测量规则，不称论文验证标准。

原始程序分数、身份适用性和模块准入分开。主适用域按Cytopus逐程序的生物父类型与已支持的参考身份匹配，使用精确词表。Spectra默认文件的归属记录原应用分组，不能一概当通用生物边界：资源明确TNK的程序可用于适用T与NK上下文，资源B的程序可用于适用B上下文；原default与资源父类型交集保留为保守敏感性。leukocyte程序不能因被打包在global中而推广至上皮/基质，CD8-T专属资源也不推广为NK。未核实缩写父类、未匹配资源或缺少父类信息先hold，保留原始分数。初始v1交集视图保留，主层级视图v2在审阅状态特征/临床效应和训练状态模式前冻结，不改分数。[Cytopus官方说明](https://cytopus.readthedocs.io/en/latest/)明确按细胞类型提取程序，选择global集合应覆盖实际输入类型。[Spectra论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10958532/)提供方法、知识库和配套输入集合依据；资源匹配不冒充每项最初扰动论文或本项目功能验证。

新状态学习继承已有冻结来源角色，仅继承development/lockbox分组，不套用旧T/NK专用资格或旧患者数。固定程序量化及数值审核可覆盖两类；特征挑选、参数选择和状态模式学习只用开发来源，冻结后才审阅留出来源的生物学读出。图谱和身份框架此前在完整来源建立，因此这属于冻结图谱条件下的状态学习留出，不称端到端外部或独立患者验证。

适用性执行配置显式设置`scope_policy_mode=resource_hierarchy_primary`；旧交集另存敏感性矩阵，不能从某个历史终态字符串推断当前策略。后继读取验收清单中的`scope_policy_mode`、精确输入SHA和两份适用mask。

开发侧读出先应用来源基因测量mask及上下文适用mask，按资格工作单位×参考身份上下文保存均值。至少20细胞才释放该读出，少于20仍保留计数及NaN；单位必须通过既有材料面板和共同域至少100细胞才参与可靠性选择。每单位该上下文至少40细胞、同来源同上下文至少10单位时，进行20次确定种子的随机半分；20细胞左右的各半不跨材料，不把患者或不同身份混成一个读出。Spearman用标准平均秩，常量列记不可估计。中位相关0.6和正向比例0.8是本项目初轮描述性测量门槛，非功能验证或临床显著性阈值。RNA量、检测基因数关联和基因集合重叠另外交付，深度关联不自动剔除，后续调整/稀释及去冗余复核仍须完成。留出分数数组不参与这一步筛选。

1. **身份内状态模式分群→CoVarNet**：在可靠状态向量上形成可重现的状态模式分区，以身份×模式作为子群、身份为分母。保留完整连续向量和多轴高低解释；分区只是计算表示，不宣称生物状态互斥。与身份CM在相同材料和细胞预算下比较。
2. **连续程序→状态协调模块**：使用既有Pelka式程序读出和适用DIALOGUE/多细胞MOFA接口，保存观测mask。连续强度不是细胞比例；有符号残差不塞进NA置0/负值截断的CoVarNet入口。

检验状态模块相对身份组成/细亚型比例的新增信息，并检查身份与状态使用基因的重叠。细胞毒/调节等自动标签本身不能充当独立状态验证。空间、临床和LR作为后续冻结候选的检验，不回流挑状态字典。

## Goal交付与验收

| 里程碑 | 必须交付 | 验收 |
|---|---|---|
| G0 | 当前输入、单位桥、资格/QC/材料覆盖 | 行轴、分母、重复材料和排除原因可追溯 |
| G1 | CellTypist、SingleR及共识/粒度视图 | 官方执行、词表/成员crosswalk、不确定项和参考范围审核 |
| G2 | 冻结组成输入与CoVarNet各方案 | 相同分母、缺测不填阴性、真实退出与独立数值读回 |
| G3 | CM匹配、拆分合并、稳定性及情境支持 | 负对照一致、允许未匹配、块重采样与来源检查 |
| G4 | 状态字典、可靠读出及状态模块比较 | 身份/状态分开、相同可测域、相对组成的新增信息 |

持续维护同任务、daily/run日志和私有执行清单；公共仓库发布经审查的通用方法/技术节点，个体记录、矩阵与模型参数留在服务器。当前goal只有在上述交付和验收均完成后才完成。

依据：[CoVarNet原研究](https://www.nature.com/articles/s41586-025-09053-4)、
[作者代码](https://github.com/QiangShiPKU/CoVarNet/blob/main/R/CoVarNet.R)、
[SingleR官方聚合注释及其混合群边界](https://bioconductor.org/books/release/SingleRBook/advanced-options.html#cluster-level-annotation)、
[SingleR诊断](https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html)、
[UCell作者文档](https://github.com/carmonalab/UCell)。

### 状态学习执行补充

连续模型采用作者官方MOFA Python0.7.5，公开wheel按SHA冻结并放在独立任务vendor，不更新全局环境。主输入保留来源/癌种差异；来源中心化和同来源内RNA量、检测基因数及身份细胞数残差作为敏感性。三种估计口径使用同一观测mask，缺失为NaN。主rank10，预设5/15和两个种子对照，不以临床或既有CM效应挑rank。模型待原生基因空间/不重叠基因诊断及合成API验证成功后运行；拟合成功不等于状态功能或独立验证完成。参照[作者Python模板](https://github.com/bioFAM/MOFA2/blob/master/inst/scripts/template_script.py)和[方法FAQ](https://biofam.github.io/MOFA2/faq.html)。

身份内模式路线预定以来源/单位均衡的开发细胞参考，在固定程序向量上标准化、PCA、15邻居和Leiden（分辨率0.3/0.6/1.0、两个种子；0.6为主）。项目预算上限为每单位100、每来源2000、每上下文20000个训练细胞；不按程序方差排名选基因集。测量不全的特征先按开发侧原生基因覆盖处理，无可支持模型的身份保留原身份。固定图上的Leiden分区保存；其余细胞以固定PCA和15近邻距离加权投影，弱邻居一致性另标。模式先编号，不先命名为耗竭、调节或耐药。参照[Scanpy Leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)和[sklearn近邻分类器](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html)；项目采样与分辨率预算不冒称作者推荐的通用阈值。

CM共同K4稳健性采用100次候选患者/材料块重采样、来源等权和留一来源重新拟合。以冻结输入的minmax参数和nsNMF方法为条件，保存raw W/effective WS的一一匹配、未匹配及固定基KL投影活动。候选块不自动视为独立患者，冻结图谱与注释不重拟合，因此不称端到端外部验证。

### Common-cell budget and source-context recurrence

The common-rank first panel uses50 prespecified paired draws with20 cells per broad parent and qualification unit. All annotation arms use the exact same cells and retain fallback cells in full denominators. Author minmax is re-estimated in each draw; freezing a positive original minimum could create negative sampled coordinates, so this control does not silently clamp negative values. Identity, rank and node vocabulary stay fixed; sampled constant nodes remain explicit.

Frozen network edges are deduplicated within a scheme across shared modules, seeds and author/raw-BH variants. Source/cancer/material/phase Spearman effects require10 materials and8 candidate blocks, with2000 exact integer block-bootstrap draws and an80% finite-replicate gate. Weighted midranks reproduce explicitly expanded unequal-size blocks and tied zeros. Same recorded cancer/material/phase summaries and broader cross-context source support are separate. Bootstrap intervals are descriptive after selected-edge discovery; they are not relabelled as multiplicity-adjusted confirmation. Source dependence, unknown design and uncertified patient blocks remain recorded.

For state increment, follow the [frozen sender/receiver validation protocol](TME-state-module-increment-validation-20261030.md). Original training float32 matrix layout and thread budget are part of numerical PCA replay; model files and thresholds are retained when correcting an audit operation mismatch.
