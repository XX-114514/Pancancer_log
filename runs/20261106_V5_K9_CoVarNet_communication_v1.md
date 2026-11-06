# 20261106_V5_K9_CoVarNet_communication_v1

范围、版本、冻结参数、实际完成/待完成及复现偏差见[报告](../docs/reports/20261106_V5_CoVarNet_communication/README.md)。任务记录见[日志](../logs/codex/2026/2026-11-06_V5-CM-CoVarNet-communication.md)。

运行目录 `${PROJECT_ROOT}/runs/20261106_V5_K9_CoVarNet_communication_v1`。两套现有环境：scanpy用于输入/统计，cpdb的CellPhoneDB5.0.1用于主计算。数据库v4.1.0。命令：`SCANPY_PYTHON=... CPDB_PYTHON=... python scripts/supervisor.py`；监督器PID、逐阶段UTC起止/退出码留logs，分组complete.json包含输出SHA。

正式流程实际启动；不是仅计划。小试20置换只技术验收，正式1000置换不继承小试效果。已完成group可恢复；患者/细胞矩阵不离开研究环境。

A全局与D细胞级响应已完成并独立读回，结果与复现偏差见报告。B完整细胞规模需跨日运行；旧supervisor计时进程由guardian接替，科学生产进程未重启。守护命令使用同SCANPY_PYTHON/CPDB_PYTHON环境执行`scripts/guard_existing_comparisons.py`。逐组SHA核验后自动读回、汇总和生成配对图；最终收据为`logs/guardian_terminal.json`及`audit/B_all18_acceptance.json`。当前没有完整B效果结论。
