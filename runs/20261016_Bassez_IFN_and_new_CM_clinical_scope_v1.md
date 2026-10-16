# 独立乳腺IFN探索检验及新CM临床范围

## 要求、起点与进展审计

继续完整细胞互作生物学/临床goal。上一goal turn完成实际分母审计、Bassez来源测量和患者资格，已改变后续动作，记progress。当前先重读Git/服务器状态：尚无新的Bassez效果目录，不把旧计划或终态当实时进程；保留另一对话正在整理的原2026年10月CM（v20261010）状态结果、daily及未跟踪文件，不消息其他任务或生成子代理。

实际按上一节点先计划再实施：固定三主要问题、完整家族/控制项和缺测规则，然后全部原生基线细胞评分、效果、独立复算和下一模块临床入口核查。原172/9CM、CRC权重、原HNSCC临床结果和权威Zarr均不改。

## 2026-10-16实际事件时间线

事件统一UTC；括号为UTC+8。用户本轮Melbourne时间另见日志节点标题，不把历史+8时刻重标为用户现时。

- 04:04:03（12:04:03）：Bassez固定公共IFN部分读出准备exit0；效果前plan已落盘，全部合格原生基线T/髓系细胞使用，无抽样/效果选群。接口从签名/rank背景剔除，原始计数/分母和来源确认列保持。
- 04:05:58（12:05:58）：三主问题BH3、未调整BH3和治疗分层Wilcoxon BH2全部8记录生产exit0。全部正负、固定类别/低df/缺测及bootstrap有效数保留；未以敏感性替代主结果，未重开有利家族。
- 04:07:31（12:07:31）：全部三问题24种控制组合/全LOO描述及跨访视标签/治疗/亚型一致性完成，无新P/BH。稀有固定类别缺失可使bootstrap/LOO不可识别；不删除该类别或患者来补出区间。
- 04:08:04（12:08:04）：全部8真实系数/t/20000患者置换/2000重采样/三家族BH独立manual-rank/full-QR复算PASS；Wilcoxon另用成对大小和SciPy U复核。
- 04:09:17（12:09:17）：全部15新图谱CM范围的原作者HNSCC精确临床基线、W/频率节点及H键审计exit0。15是scope数而非独立模块总数；四范围有基线入口，已确认原患者H多为治疗前后平均，不接受其作为基线预测量。
- 04:10:26（12:10:26）：24控制组合/完整LOO再经独立manual-rank/QR核对，完整私有证据报告exit0。报告合并全部效果、反例、区间缺失、接口RNA和下一临床入口，不公开科学数据。

## 科学决策与限制

未调整、Wilcoxon和完整主模型分别解释；本轮不足以晋级稳健两端IFNG通讯候选。控制项分解有助解释敏感性，但不能证明RNA复杂度只反映技术因素。原HNSCC巨噬单端线索保留，未命名耐药模块，未把克隆扩增当病理疗效。粗T/髓系、部分公共签名、来源保留过滤和项目独立验证资格限制都保持。

停止在同一IFNG家族追加阈值；转模块临床/队列意义。四合格新图谱scope及全components在效果前固定，保留与2026年10月CM（v20261010）原9模块的区别。下一只用基线频率及精确原minmax/W投影；原发现使用过治疗后样本，故后续只能称内部探索关联，不能称独立预测。T-only与跨谱系范围分开。详见[方法与下一轮合同](../methods/Bassez-IFN-and-new-CM-clinical-followup-v1.md)。临床新模块效果尚未运行。

## 失败、验证及复现

本轮科学生产和审计均首版exit0，无重启或删除。Nature正文web取读500/跳转不可达，改用原论文PubMed和搜索索引确认克隆扩增定义，没有绕权限。一次本机只读glob未匹配，后用rg --files定位；不影响科学输入或计算。没有新表达下载或外部消息。

环境：服务器scanpy Python；PANCANCER_PROJECT_ROOT/PANCANCER_PREFLIGHT_ROOT给真实路径。代码：prepare_Bassez_IFN_validation_v1.py、test_Bassez_IFN_validation_v1.py、audit_Bassez_IFN_validation_v1.py、diagnose_Bassez_IFN_controls_v1.py、preflight_HNSCC_new_CM_clinical_scope_v1.py、report_Bassez_IFN_validation_private_v1.py。每版新目录拒绝覆盖；plan/terminal记录时间、参数与输入/代码SHA。

私有逻辑产物为${PREFLIGHT_ROOT}/Bassez_IFN_validation_preparation_20261016_v1、Bassez_IFN_validation_results_20261016_v1、Bassez_IFN_control_diagnostics_20261016_v1、HNSCC_new_CM_clinical_scope_20261016_v1、Bassez_IFN_private_evidence_20261016_v1。主要结果SHA d1404539e214c7bf3e340a0ea155dbaa75ba089508ed00b440dc93de53078deb；准备表SHA d0cb84207f91fb44415b979b13dae6e9ecbbe5ecc40d2762fa10d3a3e1383402；完整私有报告SHA f987d0adafff6d218bde43ce46278bd1aea8472f927cce7048c2e5a68759a3aa。源患者/细胞行、科学统计、学习签名仅server/获准对话，不写公共GitHub。

本机源码compile、仓库隐私/链接/TSV检查与diff在发布前执行；显式stage本任务代码/日志，保留其他修改。普通commit/push及独立remote SHA核验后才称云端已同步；收据单列。goal仍active。建议commit：analysis: complete fixed breast IFN validation and audit baseline module clinical scope。
