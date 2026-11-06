# V7新增QC与clean refit规则（2026-11-06）

范围和已观察结果以[登记](../inventories/analysis/2026-11-06/v7_qc_clean_reannotation_v1/README.md)为准。用户明确授权清理，新增视图保留原对象、原标签和历史。

处理顺序：旧联合候选→不同组件锚点复核→剔除/holdout/保留互斥掩码→独立输入审计→干净谱系HVG/PCA/图→RNA DE和marker复审→新版本身份侧表。竞争RNA可能来自doublet或污染，本操作性排除不声称真值。低表达、IFN、增殖、未测marker和缺少doublet分数均不是独立删除理由。

正向锚点：T用CD3D/E/G与TRAC/TRBC1/2；B用MS4A1/CD19/BANK1与CD79A/B；Plasma用MZB1/DERL3与JCHAIN/SDC1；Myeloid用CSF1R/CD14/CD163/C1QA与LST1/MS4A7；Stromal用COL1A1/2与DCN/LUM；Epithelial用EPCAM与KRT8/18/19/14/17；Endothelial用VWF/CLDN5与PECAM1/KDR；Platelet用PF4/PPBP与GP9/ITGA2B；Mast用TPSAB1/TPSB2与CPA3/KIT/MS4A2。每组≥一个marker各≥1rawUMI，加既有竞争程序门槛。共享单基因不当作完整谱系证据。测量缺失不当阴性。

原not_run：按完整记录QC单元（来源×样本×原检测方法）用Scrublet0.2.3补查，≥100细胞；seeds0/37/20261106、expected_rate=.04、min_counts2/min_cells3/min_variability85、PC=min(30,max(5,n/5))，精确近邻。3次成功、各可检测比例≥.1/调用率≤.15且极差≤.05后审阅直方图；3/3同细胞阳性才考虑剔除。有invalid-log警告的候选暂存。这是项目保守规则，没有把95分位或稳定调用解释成95%doublet概率。方法依据[Wolock 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/)、[原作者代码](https://github.com/swolock/scrublet)及[Scanpy](https://scanpy.scverse.org/en/latest/api/generated/scanpy.pp.scrublet.html)。同型doublet仍可能漏检。

干净拟合：完整RNA库LogNormalize1e4；来源≥80%可测基因中VST2000HVG（沿前诊断的globin候选过滤）、50PCA、30PC图；来源gse_id作为Harmony敏感性校正变量（不冒充实测技术批次）；raw res.4与source res.2/.4/.6/.8，seed20261106；不回归细胞周期/状态，不使用旧身份建图。RNA Wilcoxon DE仅用100%保留来源可测基因，min.pct=.1/min.diff.pct=.25/logFC=.25/return.p=.01。稀有或未共有的marker仍可用来源限定正向表达复查，不能因未进全局DE而认定不存在。

marker图必须上方标对应身份/候选亚型/状态，下方标marker，左侧标工作注释与n；点大小用来源可测细胞分母，颜色为基因内相对表达，覆盖<80%标×。MC_T/MC_TC蛋白定义只作为RNA候选；Treg需FOXP3等核心及T组件支持，分类器仅佐证；CXCL13或细胞毒程序不能单独决定Tfh或CD8身份。
