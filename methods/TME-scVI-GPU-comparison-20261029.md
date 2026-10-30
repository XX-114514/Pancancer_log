# 当前TME的GPU scVI平行整合

用户要求同时使用scVI和服务器GPU。本分支与已启动的Harmony并行，不串联两种校正。细胞轴、来源资格、既有QC与共享HVG从本次冻结输入继承；临床结局和CM效果不参与训练或调参。

依据[历史多方法运行](../runs/20260725_120105_unified_pancancer_covarnet.md)、[官方SCVI接口](https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html)和[官方流程](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html)。优先使用已有CUDA环境，先核实际版本、API与小GPU算例，不升级共享环境。

输入为当前相同拟合池和标准HVG的发布整数counts，不输入log值、PCA或Harmony坐标；不对浮点去污染估计取整。特征全部要求当前来源唯一实测。来源作为batch，患者、癌种、治疗和疗效不自动当作要消除的协变量。Harmony使用已有采用的分析计数，scVI似然使用整数counts，两条标准路线的输入层和库大小处理不同；因此差异不能全部归因算法。

沿用历史可核配方：latent30、2层、hidden128、batch2048、最多30epochs和固定seed；likelihood/dispersion、拆分和实际默认核接口后冻结。随机训练验证拆分是优化监控，不是跨来源独立验证。不按UMAP或临床P值追调。

实际旧脚本核到gene_likelihood=nb、默认gene dispersion和90%训练拆分，本轮据此冻结。训练计划复核发现当前默认KL warmup为400epochs，不能把30epochs称为已完成全部正则化训练。首轮预算仍为30epochs，显式10epochs线性KL预热，后20epochs完整权重；这是本轮预先登记的预算适配参数，不冒称软件默认或旧配方原样复刻。不启用提前停止，以保证完整目标训练阶段；保存每epoch的KL权重和损失，终态要求KL权重到1。原默认及修改依据见[官方TrainingPlan](https://docs.scvi-tools.org/en/stable/api/reference/scvi.train.TrainingPlan.html)。

使用单张可用CUDA卡，显式GPU训练并检查模型参数实际设备；不可用时不静默退回CPU。核调度方式及他人占用，保留所有进程。检查RAM、硬配额/物理空闲及既有预留；仅保存模型、历史和latent，不复制完整表达矩阵到磁盘。

完整latent轴和有限值、模型保存重载通过后建立scVI近邻图，用同图Leiden/UMAP并与Harmony采用一致查看口径。来源混合与生物保留共同审阅，大类边界、稀有群、肿瘤程序和T状态结构另审；二维美观不作胜出标准。原生表达用于marker，latent不作DEG表达矩阵。旧标签及日期CM保持。

## 2026-10-30实际GPU执行恢复

调度器宣告GPU占用与全局运行作业/物理显卡空闲记录不一致，原Slurm分支没有进入模型训练。用户允许取消较早任务，并在当前对话明确转述管理员允许直接GPU运行；仅取消被替代的本轮排队条目和用户指明的旧数组，保留文件与其他任务。新的直接执行分支使用单卡可见性、实际空卡检查和任务锁，不伪造Slurm作业ID，不改共享调度配置；数据、架构、NB似然、种子和训练计划不变。

首个直接执行版本的GPU小训练成功，保存重载时device=0被当前Lightning解释为零张卡而失败，历史保留。新版本按既有脚本用device=auto，单卡可见性不变；真实GPU训练/保存/重载算例通过后，全量正式训练已启动。实际CUDA设备、每轮损失、KL权重和显存记录用于验收，不用进程存在或显卡可用代替模型完成。

后继独立审计核全部latent、完整训练历史及每轮实际CUDA设备，从保存模型导出真实训练/验证暴露索引，并用原生计数在CPU重载复算有界探针。这里CPU仅用于独立推理读回，正式训练为GPU。随机细胞验证不能称独立患者或来源验证；生物整合效果和最终注释仍需后续审阅。
