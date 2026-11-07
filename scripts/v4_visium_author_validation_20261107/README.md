# 同V4标签的Visium验证脚本

见[方法](../../methods/V4-Visium-CoVarNet-author-validation-20261107.md)及[实际run](../../runs/20261107_v4_Visium_author_validation_v1.md)，`code_manifest.json`绑定公开字节。没有附私有输入，Cloud不能假定真实分析已重跑。

- `prepare_visium_reference.py`/`audit_reference_V4.py`：已执行并独立通过的同V4参考及原始计数/源覆盖准备和核验。
- `cell2location_V4.py`：保留的训练前devices重复失败版本；`cell2location_V4_v2.py`是当前GPU生产版本，不改科学参数。
- `measure_reference_separability.py`：基因面板上的同标签留出重放，仅诊断，不重定义注释。
- `score_author_Visium.R`：直接提取固定作者共定位/60近邻有向Moran函数；全部CM及两posterior完整保留。
- `audit_score_Visium.py`：独立矩阵分数、距离权重、完整epoch/丰度审核，同谱系和参考谱诊断，以及成员W活动；原文指标和补充控制分开。
- `summarize_Visium_validation.py`/`audit_Visium_cohort.py`：全部请求来源/14CM的汇总与独立读回；不推断患者独立、不汇合spot做P值，原空CM为NA。
- `spatial_math_fixture.R`/`test_author_spatial_math.py`：明确合成的100spot公式重放；实际恢复版本补入Path导入，仅修测试器，不把此测试称真实空间结果。
- `render_Visium_report.py`/`finite_report_controller.py`：仅在真实全家族独立审核通过后生成全部CM汇总图与按输入顺序选定的CM03投射图；私有输出，需视觉/科学审阅后再公开，不自动push。

准备/审核/资格/汇总通过`runpy.run_path`注入`PROJECT_ROOT: pathlib.Path`；GPU脚本参数为`RUN_ROOT reference`或`RUN_ROOT section SECTION_CODE`；R指标及独立分数审核接收`RUN_ROOT BASE_CM_RUN SECTION_CODE`。私有`section_jobs_PRIVATE.json`绑定已资格化的实际输入，不能擅自换原始矩阵仍称同一运行。

运行版本v2保留原50标签（5个上皮仅背景参考）、原45 CM背景节点。输入参考已通过不等于模型/评分/总体审核已结束。任务PID、连接、物理服务器路径、矩阵/模型、逐细胞/spot/患者行和完整日志不公开。有限controller仅串接科学阶段，没有自动发消息、自动公开推送或长期计划任务。

后续：`audit_fitted_reference_v2.py`已对真实250轮/50标签签名独立读回通过；首版读本任务自产pandas object类型npy失败，v2明确转float，原文件与失败代码保留。逐切片正式审核改用`audit_score_Visium_v2.py`同一读取兼容修复，主模型/科学参数不改。

`fit_CellCharter_Visium.py`、`audit_CellCharter_Visium.py`和有限controller为独立表达—位置niche分支，不用CM效应训练。v2 fitter/audit仅为首版训练前零HVG spot失败建立合格niche视图；低信息spot仍保留原轴/主模型。有限recovery复用成功模型、只重试已识别的该类训练前失败；`summarize_CellCharter_Visium_v2.py`/`audit_CellCharter_cohort_v2.py`以实际接受分支和合格spot对齐，全范围新汇总写入`niches_final_v2`。原分支/失败/代码和原53请求均保留。模型/指标/审核的执行状态以较晚检查点为准。
