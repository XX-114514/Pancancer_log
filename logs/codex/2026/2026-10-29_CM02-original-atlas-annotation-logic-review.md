# 原免疫图谱注释流程与逻辑说明

用户要求先全面说明当前整合免疫图谱的注释流程、逻辑及结果。本轮明确范围为CM02审查包的原九来源免疫v1，不与另一工作流的广泛来源工作注释混用；没有恢复LR目标、改标签或启动新计算。

只读核对免疫v1方法、冻结图谱scope、完整群目录、实际应用面板、全基因marker候选、原生dotplot统计，以及归档的biology与marker-refinement实际代码。实际查看两页T群dotplot。原文件与并行工作保持。

关键流程核查：大类来自逐细胞原生marker面板打分与固定门槛，T/NK有专门覆盖／边界规则；全局粗聚类不等同大类标签。大类内重新建表示与聚类，子类面板群均值仅形成有门槛的候选字段；展示marker随后用本群全原生表达和大类内背景对比选择，未作正式DEG显著性。原生／整合坐标、身份／状态、群分割／细身份及候选／晋级分别说明。

本轮新识别的解释风险：细分面板是均值竞争，缺少部分身份必需锚点的硬门控；候选面板可能混有共享状态基因；低于可测基因门槛的面板分数在代码中填零，须区别dotplot的可测分母处理；major_only记录仍可保留最高分候选字符串，消费端不能忽略证据状态。展示的单marker、应激／性别等特征及群内对比不能直接升为新细胞类型。以上是代码逻辑审查与后续修订方向，不宣称已完成重新注释。

技术汇总与原目录逐项核对，计数合计、全部证据状态、细身份未晋级字段、展示marker来源角色和T群候选均以本机私有目录为据。结果数值与具体科学判读只在用户对话／既有私有审查包说明。未导出患者／细胞行或公开科学表。

下一先审大类未决／竞争与来源捕获，再把T细胞身份硬证据和状态程序分开，逐群审候选与major_only字段；若另行修订，保留原节点／标签并输出完整crosswalk，不能凭本次文字审阅直接宣称新身份已接受。建议提交：`Record original immune atlas annotation logic review`。

## 用户追问后的原始文献与流程核对

用户进一步要求研究判定逻辑的验证依据、Harmony按来源的具体含义、全量邻接图与抽样范围，以及cluster→marker→整群命名的参考流程。重读实际PCA、Harmony、neighbors、UMAP及谱系内建图调用，区分抽样拟合与全体细胞投影／建图；没有把近似kNN称为逐对精确距离计算。

实际阅读Seurat官方guided clustering流程的FindAllMarkers及RenameIdents步骤、Luecken/Theis2019最佳实践的cluster annotation章节、Harmony作者公开原论文PDF与官方接口说明，以及scIB原始benchmark。SCINA仅取得原论文索引／摘要信息，全文页浏览器验证未通过，不宣称全文已读。Nature部分HTML认证跳转不可用，Harmony使用正常公开的作者论文PDF；没有处理CAPTCHA或绕过受限数据权限。

依据：现有文献支持marker／参考辅助及cluster层人工审查，但不能验证本项目自定义均值分数、margin与阈值。现有记录未提供相应阈值的外部标注性能／校准验证；应保持初始候选身份。当前先按规则划父群、后在父群内分群，再补本群marker，无法独立纠正父群分错的问题。用户提出的cluster→本群marker→多证据整群命名有直接官方流程参考，应成为下一注释主线；群内冲突须先拆分／重聚类，不强制整群同名。marker检验用于解释群，聚类后同数据P值不作独立身份验证。

Harmony的dataset批次变量对应研究来源，各细胞连同来源标签一起进入同一次低维校正；不是各研究各跑后拼接。源码支持全归档细胞参与Harmony和全量kNN图；PCA和UMAP拟合各自有抽样，其余细胞投影。来源平衡抽样本身不是独立方法验证。下一比较未校正／校正表示，同时检查生物保留与批次混合；来源与癌种、组织、捕获混杂时不能一概视为技术噪声，也不盲目把每患者作为待消除批次。上述属于文献与源码审查和下一方案，未重跑或晋级注释。

文献入口：[Seurat官方流程](https://satijalab.org/seurat/articles/pbmc3k_tutorial)、[2019最佳实践](https://pmc.ncbi.nlm.nih.gov/articles/PMC6582955/)、[Harmony作者公开原文](https://jef.works/assets/papers/nmeth101038.pdf)、[Harmony官方接口](https://portals.broadinstitute.org/harmony/index.html)、[scIB原研究](https://www.nature.com/articles/s41592-021-01336-8)。仅补充本任务技术日志，不公开新的科学表或修改任何原注释。
