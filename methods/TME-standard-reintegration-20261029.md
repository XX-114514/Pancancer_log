# 当前TME标准HVG/PCA/Harmony重整合

用户明确要求替换自定义方差/均值排名，核实历史doublet/环境RNA处理，并针对当前来源合格TME重新整合。新版本只写独立运行目录，旧表达、计数、身份、状态、来源资格、局部群及日期CM保持。

## 方法依据

主配方沿用项目既有5M整合脚本的标准`scanpy.pp.highly_variable_genes(flavor='seurat_v3', batch_key=来源)`，6500个候选中按标准HVG排名及批次数选4000个非技术基因；技术排除沿用MT、核糖体、血红蛋白、ERCC与MALAT1规则。中心化ARPACK PCA最多50维，Harmony以来源为batch，最多20轮，固定随机种子。Harmony显式使用theta=2、lambda=1：theta=2对应已安装2.0接口的默认多样性惩罚，lambda=1采用历史版本的固定岭惩罚。两者不能统称当前API默认：2.0中lambda=None会启用动态估计，历史0.0.10中theta=None解析为1、lambda=None解析为1。本次不是旧Harmony数值参数的原样复刻，差异在运行前冻结，不按新UMAP形状优化；实际接口与输出方向先核。线程限制是计算资源设置，不是新增生物学参数。保留原始PCA与校正PCA，两者不通过二维美观选赢家。

参考[Scanpy HVG](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pp.highly_variable_genes.html)、[标准教程](https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering.html)、[Harmony原研究](https://doi.org/10.1038/s41592-019-0619-0)及[公开Scanpy skill](https://github.com/K-Dense-AI/scientific-agent-skills/blob/main/skills/scanpy/SKILL.md)。skill按其主文档与workflow只作流程合同参考，不把示例QC阈值、细胞名或另版本API直接移植进本项目，不安装/升级现有共享环境。

## 历史QC与输入层

先前“未闭合”的措辞没有充分衔接历史全量技术验收，须纠正。原5M准备已有逐样本QC/Scrublet；新增9来源/逐文库doublet与环境RNA准备有完整最终验收。缺少新汇总字段不说明作者或旧流程未做处理；处理完成、模型可辨识、是否采用估计和是否进入当前用途分别记录。

本轮不重做Scrublet/DoubletDetection/DecontX，不为覆盖率重删细胞。当前已接受来源合格拟合掩码作为输入，未拟合保留记录继续留档。存在已采用`ambient_removed_counts`的存储按`counts - ambient_removed_counts`恢复分析计数；没有该层的legacy分支使用既有发布源整数counts并保留原处理历史，不伪称全部有统一环境RNA估计。每个来源记录实际采用层、差值非负和分析总RNA；原估计未采用时不会再次扣除。

标准v3 HVG使用真正整数counts，禁止将DecontX浮点估计四舍五入。PCA使用已采用分析计数归一化表达。按官方预处理顺序，从完整已采用RNA总量normalize_total(10000)、log1p后再取HVG；这与旧脚本先取HVG再按子集总量归一化的实现不同，明确登记，避免特征选择改变深度分母。HVG选择所用counts和PCA所用分析层单列。所有特征必须在当前每来源唯一实测，不将未测基因补为生物学零。全基因表达仍从原归档读取，不物化重复全量counts。

## 运行与验收

首轮阶段只完成标准HVG、全量当前TME的PCA与直接Harmony，不复用旧全档案参考或Symphony查询坐标。冻结输入合同、原有脚本SHA、参数、实际环境与编译LOESS来源。先在实际环境运行标准batch-aware HVG小型接口测试；其通过只证明API可执行。资源检查同时取物理空闲和组硬配额较小者，预留输出/审核/工作空间；内存门槛依据既有原生nnz与HVG矩阵规模，不把整个共同基因矩阵稠密化。

分来源逐块核原始RNA总量/整数count域、行与实测特征映射，记录完整拟合输入轴及精确范围。分析层逐块核扣除值非负有限、不增加计数及完整RNA归一化；不更改已有校正标志。标准HVG表、特征名、PCA载荷/方差和全体原始/校正坐标独立保存、SHA及数值读回。技术验收与身份/状态保留评价分开，不由Harmony收敛自动宣布生物学整合成功。

后继阶段在本次Harmony表示上建立明确近邻图，Leiden编号和UMAP共享该图，原生每群marker与正负标志审阅另立证据。不会沿用之前自动模块规则直接宣称最终类型。图谱中的小群、连续状态、来源特异肿瘤程序与批次混合共同评价；后续标签或CM变化另行验收。
