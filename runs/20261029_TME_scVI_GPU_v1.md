# 当前TME GPU scVI平行整合

状态取证：2026-10-29T05:11 UTC /13:11 UTC+8。Slurm已接受GPU作业，当前PENDING Priority，尚未开始真实小模型或全量训练；现有CUDA环境及小算例已通过。排队并非CUDA失败，不静默退回CPU，不重复提交或取消其他任务。

输入为当前冻结拟合池和标准HVG的原生整数counts，与Harmony并行。NB似然、latent30/layers2/hidden128/batch2048；固定30epochs、10epochs KL预热，90%训练/10%验证，完整权重阶段不提前停止。使用既有scvi CUDA环境，不改变共享依赖。具体资源、轴、路径和作业ID在私有冻结合同；输出仅模型/历史/latent，全部表达留在原归档。

先真实GPU小训练及保存重载，再核完整原生输入后训练。逐epoch损失、KL权重、实际device和峰值显存保存；完整latent轴/有限值、文件SHA和模型重载后才进入后续同图Leiden/UMAP。该有限等待链已实际启动。最终整合评价和身份/状态命名尚未验收。

同期Harmony主计算actual exit0；修复版独立数值审核05:06UTC actual exit0，同图近邻阶段已继续。审计v1误用整数读取器验证浮点ambient差值导致失败，历史保留；新分支保留原counts整数检查。四个实际回归案例通过。详见[方法](../methods/TME-scVI-GPU-comparison-20261029.md)和[日志](../logs/codex/2026/2026-10-29_TME-scVI-GPU.md)。
