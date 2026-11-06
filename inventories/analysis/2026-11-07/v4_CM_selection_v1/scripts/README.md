# 代码与输入合同

`code_manifest.json`保留代码SHA及执行角色。私有输入沿用2026-11-06 V4 CoVarNet主v2；公开结果中不含模型/频率/W/H矩阵或个体行。

- `diagnose_all_CM.py`：全14原CM、全部18全局边、来源中心化/分母/活动诊断；以`PROJECT_ROOT: Path`注入执行。内部验证A=W·S后处理R的V1..V14序列化列名。
- `refit_CM.R RUN_ROOT BASE_ROOT JOB_ID`：按`jobs.json`执行5seed和13LOSO，各30次固定K14；不重新选择各视图最优K。
- `summarize_refits.py`：W/A一对一匹配、独立精确DP检查、全部252对应与网络恢复；注入PROJECT_ROOT。
- `audit_selection.R`：保留的首次CCC JSON精度审核失败代码；不是接受依据。
- `audit_selection_v2.R RUN_ROOT`：成功的保存模型/完整CCC/WSH导出独立核验；以1e-10比较完整模型值，对默认4位JSON仅检查其实际舍入界。
- `audit_selection.py`：独立数组级来源/分母/活动/全部匹配核验；注入PROJECT_ROOT。
- `export_edge_LOSO.py`及`trace_CM03_edges.py`：补充来源删去边明细与全部26 CM03原边筛选逐步重建；注入PROJECT_ROOT。
- `render_selection_report.py RESULT_ROOT ORIGINAL_AGGREGATE_ROOT`：只用已接受的聚合表生成全部14取舍及图；不改模型/成员/阈值。报告文字给出科学判断，不能用模板替代换输入后的重新审视；如已接收CM03逐步追踪表，则验证其26行和实际丢边原因后加入报告。

实际科学逻辑run是`runs/20261107_v4_CM_selection_v3`。v1/v2仅有因子轴检查失败历史，未启动NMF；所有原模型及原发现输出不覆盖。SSH、PID及物理根目录不公开。Cloud缺私有输入不能宣称已经重新执行本轮科学分析。
