# 20261106_V5_K9_CoVarNet_communication_v1

范围、版本、冻结参数、实际完成/待完成及复现偏差见[报告](../docs/reports/20261106_V5_CoVarNet_communication/README.md)。任务记录见[日志](../logs/codex/2026/2026-11-06_V5-CM-CoVarNet-communication.md)。

运行目录 `${PROJECT_ROOT}/runs/20261106_V5_K9_CoVarNet_communication_v1`。两套现有环境：scanpy用于输入/统计，cpdb的CellPhoneDB5.0.1用于主计算。数据库v4.1.0。命令：`SCANPY_PYTHON=... CPDB_PYTHON=... python scripts/supervisor.py`；监督器PID、逐阶段UTC起止/退出码留logs，分组complete.json包含输出SHA。

正式流程实际启动；不是仅计划。小试20置换只技术验收，正式1000置换不继承小试效果。已完成group可恢复；患者/细胞矩阵不离开研究环境。
