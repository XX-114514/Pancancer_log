# V4全部14CM的评估与选择

用户要求比较所有新CM并选择最适合研究者。当前分析以2026-11-06 V4主v2为基准，历史V5/K9与新TME图谱不混入。见[方法](../methods/V4-CM-selection-source-stability-20261107.md)。

2026-11-06 22:12:09 UTC／2026-11-07 06:12:09 UTC+8开始取证；客户端日期2026-11-07。独立逻辑run`${PROJECT_ROOT}/runs/20261107_v4_CM_selection_v3`；v1/v2因A导出轴名检查失败保留，未启动NMF。v3在验证A=W·S后显式对齐V1..V14与原CM顺序，22:15:09 UTC描述诊断完成并启动有限队列。五额外seed＋13逐来源留出，各30次固定K14，共540次，6并行、BLAS单线程；摘要和独立审核由有期限的接续执行。

输入绑定沿用原主v2的频率、分母、来源准入、W/S/H和注释；不重新载入全RNA、不修改原CM、不因某模块有趣改变来源范围。新代码及实际执行SHA随结果收据保存，研究环境R4.3.3/NMF0.28、Python扫描环境及相关numpy/pandas/scipy。公开不含连接/PID/物理路径或个体行。

全部18模型22:18:43 UTC退出0，22:18:49完成252个因子对应；22:20:01 UTC（2026-11-07 06:20:01 UTC+8）R/Python独立核验通过。[全14结果与图](../inventories/analysis/2026-11-07/v4_CM_selection_v1/README.md)。CM03网络质量65.7%，seed W余弦中位0.990、原边5/5保持；LOSO中位0.984、原边10/13，严格成员5/13一致。CM06、CM09有来源敏感的边；CM09额外seed W谱中位仅0.196；CM11单节点Glial父类造成条件频率常量。CM07网络质量仅9.7%，与其Melanocytic最高载荷不一致。CM03定为优先验证候选，CM06备选，不称当前已发表级模块。

R首次核验因producer JSON的CCC默认四位小数与完整精度严比失败；保留失败收据，v2直接按模型完整精度独立计算CCC，以1e-10对模型原生值、5e-5舍入界对原JSON核验，18个最大序列化误差4.85e-5，WSH/导出及每模型30次全部通过。未重跑科学模型。全NA空图指标引发numpy空切片警告，NA语义保留。原CM、原注释和既有Git修改未变。

复现合同：`diagnose_all_CM.py`通过`runpy.run_path`注入`PROJECT_ROOT: Path`；`refit_CM.R RUN_ROOT BASE_ROOT JOB_ID`按`jobs.json`执行；`summarize_refits.py`、独立`audit_selection_v2.R RUN_ROOT`及注入PROJECT_ROOT的`audit_selection.py`后验收。代码见结果内`scripts/`，不含连接或物理路径；`audit_selection.R`保留为首次舍入比较失败版本。当前所有结论是发现后研究优先级；固定K14没有给各留出重新选秩，也没有独立空间/临床/机制或患者层确认。
