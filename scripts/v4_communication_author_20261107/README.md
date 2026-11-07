# 当前V4的原文通讯分析脚本

[方法](../../methods/V4-CoVarNet-communication-plan-20261107.md) · [实际run](../../runs/20261107_V4_CoVarNet_communication_v1.md) · [冻结配置与代码SHA](../../inventories/analysis/2026-11-07/v4_communication_author_v1/method_freeze.json)。

准备/输入审核/响应/汇总/总审计使用`runpy.run_path`并注入`PROJECT_ROOT: Path`；CPDB和其审核接收`RUN_ROOT TAG`，R接受`RUN_ROOT`。有限控制器由外部环境变量`RUN_ROOT`/`PROJECT_ROOT`定位，不含实体服务器路径和连接凭据。

原始计数准备保持全部既有V4细胞/标签/CMT；13来源共同实测基因防止结构性缺测当0。独立输入审计抽原CSR重算标准化；技术小试不可冒充生产结果，正式均1,000置换，正式High/Low不截断到每标签1,000。完整25组和所有14CM保留，空网络无内部结果，缺失细胞对NA。

CellChat独立算法和直接Scanpy/Immune Dictionary响应同时保留；超几何为等价重建、签名映射/方向/背景与CM最小q的含义分开。主空间结果尚未全范围完成，通讯执行是预登记探索，不自动获得高共定位CM主结论资格。

模型、h5ad、mtx、个体行、患者/样本桥与完整日志留私有。此目录是可审阅科学代码，不包含私有输入，不能仅凭本代码宣称Cloud已复算真实结果。

实际使用`audit_V4_CPDB_v2.py`：先核对唯一互作ID集合完全相同，再将means/significant_means按pvalues的ID顺序重排；原版同序assert失败保留。`finite_V4_communication_controller_v2.py`复用成功输入审核和技术小试，重新审核既有小试输出后从正式global接续，不重跑已成功科学阶段、不改参数；见检查点的`execution_code_recovery_v2.json`。初始冻结code_manifest保留，恢复另行绑定SHA。

`audit_completed_global_response.py`独立审核真实全局资源覆盖、全部可用身份的直接Scanpy DE/超几何/BH和生产表达；`review_cytokine_symbol_aliases.py`只做IL17E/IL25 canonical符号审阅、不重跑DE；两者部署SHA另在聚合审核目录绑定。`render_completed_communication.py`在本机读取公开聚合目录，设置`CCI_REVIEW_ROOT`后生成两张静态科学图，不访问私有数据；输出图已独立视觉检查，初版图例遮挡修正后重新渲染。
