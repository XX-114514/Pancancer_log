# TME GPU scVI节点

05:00UTC /13:00 UTC+8：用户要求GPU scVI，核实际CUDA环境/API/旧脚本；旧环境依赖冲突保留，另一现有环境GPU算例可用。NB/架构来自实际旧代码，30epochs预算与默认KL预热不匹配，本轮显式10epoch预热。

05:06UTC /13:06 UTC+8：Harmony主阶段已真实完成，独立审核v2通过，同图近邻v2继续；原审计对浮点delta的整数约束错误保留并修复，四个count/delta回归通过。

05:08UTC /13:08 UTC+8：GPU scVI通过Slurm实际提交；05:11UTC取证PENDING Priority，未称真实GPU训练已经执行。单次SSH返回255后先确认无目录/作业，再新会话提交，避免重复。训练后同图链已启动等待。见[run](../../../runs/20261029_TME_scVI_GPU_v1.md)。
