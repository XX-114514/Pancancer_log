# HNSCC亚型内程序与独立验证入口

## 请求、起点与计划

继续完整细胞互作生物学/临床目标。上一goal turn为progress：HNSCC预设临床/配对家族及独立复算完成，仍不足以宣称两端联合或直接互作。此次按全部巨噬候选亚型分解、全部50独立功能读出、元数据/外部覆盖与方法审阅推进；原9CM/172身份/CRC模型均冻结，保留并发atlas/state工作。

## 实际时间线（事件UTC；括号UTC+8）

- 04:11:28（12:11:28）：所有现有ICB补充来源原字段/细胞谱系/全局基因覆盖核查完成；没有程序或结局检验。
- 04:13:44（12:13:44）：全部Macro_多基因候选亚型内固定程序、捕获及3描述变体/LOO完成；细胞选择独立于结局，原绝对权重95%门槛和质量系数不变。
- 04:14:29（12:14:29）：ICB精确稳定ID桥接和来源患者/时点/结局/端点核查v1完成。审查发现默认字典创建不存在ID键，影响库存展示；04:18:47（12:18:47）v2改get读取，原映射/覆盖/分母数值逐项保持，未运行效果检验，旧版保留。
- 04:18:27（12:18:27）：全部50 Hallmarks×3亚型×2描述变体真实完成；训练基因从签名和排名背景同时移除，源可测轴/固定细胞与平均并列秩maxrank1500保持，独立原生计数排序标量核对PASS。
- 04:19:51（12:19:51）：完整功能结果触发性别/年龄诊断，全部临床端点/亚型/50功能保留；仅描述，不新增P/BH、不替代原主分析。人口组成不均衡限制保留。
- 04:23:11（12:23:11）：全部471描述性效应及全部留一患者范围独立manual-rank/QR复算PASS，包含正负及低df不可估计情形。
- 04:29:16（12:29:16）私有完整解释报告v1完成；04:30:52（12:30:52）通用展示代码v2完成，数据依赖的实际结果由服务器表填充，科学输入/效应未改。

以上为2026-10-14真实终态；日志同步时间见本轮节点标题。每次实际production/audit均exit0；观察等待只重读原句柄，没有重启科学任务。

## 已完成证据与边界

亚型内变化、功能与人口/质量诊断已审，完整正负/覆盖不足保留。没有新增显著性检验、正式状态标注或通讯机制宣称。学习权重与患者/细胞行、派生科学统计只留获准对话和服务器私有报告；公开方法无实际统计或学习基因列表。

外部全局覆盖只是来源测量上界，原95%门槛不降；旧ICB补充的粗T/髓系不能充当原细节点映射，最低细胞数保留过滤不能当总体组成。Bassez原结局为T克隆扩增而非病理退缩，不能跨研究把Combined_outcome统一为疗效。隔离来源不进入验证。

实际访问作者新版网站与公开HTML，未发现处理后count直链；旧Biokey地址web工具不可达，普通作者页面HTML成功。另读Franken公开空间仓库，确认受控原始FASTQ及公开说明/清单，并非已获取可运行空间矩阵；没有访问申请、作者联系、权限绕过或新表达下载。

已读scIVCCC原论文：治疗工具变量/条件独立的组织可借鉴，但当前没有保证治疗仅经发送者影响接收者的证据，活检/切除差异仍在；不因条件独立零假设未拒绝就宣称中介。详见[完整方法与未来入口](../methods/HNSCC-programme-subtype-and-external-validation-v1.md)。

## 产物与下一步

scanpy Python运行如下公开方法代码：diagnose_HNSCC_within_macro_subtypes_v1.py、diagnose_HNSCC_full50_orthogonal_functions_v1.py、diagnose_HNSCC_sex_function_sensitivity_v1.py、audit_HNSCC_subtype_function_followup_v1.py、probe_ICB_programme_validation_metadata_v1.py、audit_ICB_clinical_programme_scope_v1/v2.py、report_HNSCC_subtype_functions_private_v2.py。环境变量PANCANCER_PROJECT_ROOT/PANCANCER_PREFLIGHT_ROOT提供服务器根，所有输出拒绝覆盖。

私有逻辑目录${PREFLIGHT_ROOT}/HNSCC_within_macro_subtype_diagnostics_20261014_v1、HNSCC_full50_orthogonal_function_diagnostics_20261014_v1、HNSCC_sex_function_descriptive_20261014_v1、HNSCC_subtype_function_independent_audit_20261014_v1、ICB_clinical_programme_scope_20261014_v2、HNSCC_subtype_function_private_evidence_20261014_v2。各版本记录计划/终态、代码/参考SHA，原报告v1保留且不公开实际解释载荷。

下一先完成T端及跨来源的完整程序解释，避免只有接收者的结论；核合法、完整RNA的独立生物学/临床/空间入口。若采用新可转移程序，需另开结局盲发现/患者留出版本和完整固定来源/槽位家族；不能改旧权重使原测量失败通过。原CM、程序协同、单端临床线索和直接互作仍分别登记，goal保持active。

轻量验收包括代码compile、仓库隐私/链接/TSV检查及diff；普通commit/push独立核云端，具体发布时刻保存在receipt。建议commit：analysis: audit HNSCC within-subtype functions and external validation gates。
