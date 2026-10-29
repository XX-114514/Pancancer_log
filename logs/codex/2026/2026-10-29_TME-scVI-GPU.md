# 当前TME加入GPU scVI

用户要求后续校正使用scVI及服务器GPU。已读仓库状态、路线、待办、旧多方法运行和全范围方法，并核官方接口；保留已有修改。Harmony原分支继续。

计划先核现有CUDA环境、GPU占用和调度状态、实际API与小GPU算例；沿用历史latent30/layers2/hidden128/epochs30/batch2048。冻结同一细胞和HVG，整数counts独立训练，不串联Harmony。两者输入层与库大小处理不同，不称纯引擎因果比较。

只读资源/环境预检已发出，尚未声称真实scVI训练开始或完成。模型与latent留服务器，通用代码/方法/技术日志按节点维护；最终身份和新CM尚未发布。见[方法](../../../methods/TME-scVI-GPU-comparison-20261029.md)。

05:00UTC取证GPU空闲、Slurm分区idle；现有scvi环境1.4.1/Torch CUDA12.6可用且真实GPU算例通过。scvi_py311存在AnnData SparseDataset接口冲突，保留，不修共享环境。实际旧06_full_scvi.py已读，确认NB似然；实际TrainingPlan默认400epochs KL预热与30epochs预算不匹配。冻结10epochs预热/30epochs完整训练、逐epoch记录KL权重，不称默认或最优。

另核当前Harmony主阶段04:24UTC actual exit0；独立审计v1把浮点ambient差值用整数counts读取器校验而失败，同图阶段按门槛中止。新审计v2仅修正独立非负浮点读取，原counts整数约束保留；新同图v2等审计后继续，未重跑Harmony主计算。失败目录保留。

05:06UTC /13:06 UTC+8修复版独立审核实际通过；同图近邻v2已经开始。四个浮点delta/负值/非有限值/整数counts回归实际通过，新增脚本AST与diff检查通过。首个GPU提交观察SSH255无输出，独立新连接确认目录不存在后再提交，未重复作业。05:08UTC Slurm接受；05:11UTC仍为Priority排队，无实际GPU训练终态。队列/节点只读核对未改变其他任务；训练后同图Leiden/UMAP的有限等待链已启动。方法、代码、run及本日志待技术同步，来源行/矩阵/学习参数保持私有。
