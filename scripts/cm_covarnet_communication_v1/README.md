# 历史V5/K9通讯补分析执行

仅用于冻结20260825身份与CM，不自动重注释/重发现。先读[方法与状态](../../docs/reports/20261106_V5_CoVarNet_communication/README.md)。

1. 在研究环境设置PC_ROOT；运行qualify.py生成独立run与私有输入核验表。不要把private目录同步到代码仓库。
2. 将本目录脚本复制到RUN_ROOT/scripts，将报告配置复制至RUN_ROOT/config；公开原始参考资源按freeze_receipt.json中名称和SHA放入RUN_ROOT/sources。数据库官方v4.1.0，运行软件5.0.1。签名需原始S3和MGI表，通过prepare_signatures.py生成。
3. qualify_metadata_v2.py、prepare_global.py、audit_numerics.py；CPDB环境执行run_cpdb.py pilot，再用audit_pilot_v2.py验收。未通过不启动正式流程。
4. prepare_full.py在scanpy环境中运行；原生读取保留完整SHA，正式High/Low不抽样。设SCANPY_PYTHON及CPDB_PYTHON后执行supervisor.py，锁阻止重复提交，阶段/每组完成收据用于恢复。
5. 审阅原生基因背景、NA分母、全部CM和所有阴性项；summary_final和图件仅技术产物，需独立科学验收。不要将细胞对P值称患者P值。

该执行器单阶段超时24小时。若长正式比较超时，保留已完成各组与失败收据，核对进程确已结束、资源与方法冻结后恢复；不能把超时或缺少输出当阴性。只改执行顺序/资源也应记录，不能据新效应改组或阈值。
