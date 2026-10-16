# 头颈癌作者亚型/状态计数临床研究

## 起点、执行和输入验收

接续[完整功能临床阶段](20261016_HNSCC_longitudinal_function_clinical_v1.md)及[用户Table S1](20261014_HNSCC_author_Table_S1_input_review_v1.md)。原作者状态层在17:07 UTC完成初始构建，当时尚未审阅新功能临床效应；当前goal继续后先核完整结果并收束原家族。新科学问题与范围见[方法合同](../methods/HNSCC-author-state-count-clinical-v1.md)，不是旧家族改阈值重跑。

2026-10-16 17:14:32.502428 UTC /2026-10-17 01:14:32.502428 UTC+8启动独立状态计数审计，代码166f119cf4a6ebe851d35d963c6991d5fe8d713071ba55d374d15aa9a833a80d；17:14:33.485228 UTC全部原作者计数、35定义/1995比例/285父分母/35单端/241组合机会通过。输入producer终态SHA fcb2660b37299e840f947813f5acd3c94054a5903b8dd2a439b45adb47a3db4a；不含逐细胞映射或新临床效应。

17:16:09.436141 UTC /01:16:09.436141 UTC+8启动全部276临床槽预检，源55c02b64eb169a052fc5961d3900ff3e85467f982537723b299e7e049a5cda8a；17:16:12.087733完成，合同4584367070fbab0e01de6cf03996e68f892cd7e5ecfa7cae144b0f0cf129871c。17:16:52.823819开始原计数独立重建，17:16:55.724278 UTC全部58数值设计/276槽通过，证明cbb2a7498d16f17efd91e1da88b1d7daff1477f68a98409ed590a2d5a4fd3fcd。全部零计数、患者交集和七列协变量保留，无新真实效应。

## 校准失败与限定修正

17:17:22.214019 UTC启动v1校准，17:17:36.282781实际FAIL，终态5c0aae6fdf24ab344f119ecee83968cb98026239626b33c453b7cdcd8812e706。17:17:56.258904派发的门控生产因该FAIL明确退出，未创建真实结果目录；17:18:07.106574伴随审计因原生产进程缺失退出，原错误/目录保留。17:18 UTC实际/proc及日志取证，不是工具观察超时，没有重启这些句柄。

新效果前只评估一次HC3受限零残差修正；v2内核38b08ccf9af8144a14bee222400ad3f44c161e8771f5b6c723356be6810112f9，17:19:16.109641 UTC /01:19:16.109641 UTC+8校准v2启动，代码5d060d6b113d91ee9d7f45084cac6ca5b7058e5b54613d0dd617ebd395a414b4。每情景4096独立新随机流、原三代表设计和.075门槛固定；尚未核最终校准，无新真实临床效果声明。失败则只描述本计数层，不第三次改法；原范围/协变量/身份与历史不变。

## 产物、命令及边界

scanpy Python，OMP/BLAS单线程，使用PANCANCER_PROJECT_ROOT和PANCANCER_PREFLIGHT_ROOT。脚本依次为preflight_HNSCC_author_state_count_layer_v1.py、audit_HNSCC_author_state_count_layer_v1.py、preflight_HNSCC_author_state_clinical_v1.py、audit_HNSCC_author_state_clinical_inputs_v1.py、calibrate_HNSCC_author_state_clinical_null_v1.py/v2.py，后续真实run/audit必须服从校准。版本目录HNSCC_author_state_count_layer_preflight_20261016_v1、HNSCC_author_state_clinical_preflight_20261016_v1、HNSCC_author_state_clinical_calibration_20261016_v1/v2。收据绑定源码SHA及PID/startticks，原表只读。

无新下载/删除、无逐细胞正式状态或原172/2026年10月CM（v20261010）变动。病例、科学表/效果/图只server及获准当前对话；公共只代码、方法、轨迹和SHA。技术审计不等于发现，整体goal active。普通commit/push必须在diff和仓库验证后，云SHA另行核验。

## 17:20—17:22 UTC校准及完整推断验收

v2原句柄17:20:07.298650 UTC /01:20:07.298650 UTC+8通过所有预定工作模型诊断，终态12b9eb353a9b2bdb89000ab100e29504af8372cc867b1b1bad56d60c7ccf41a5；正式修正合同e22b5291f9adab71675493ed761babd3119519baba4132f74a735012d590a420。原v1失败未删除或改写，真实临床效应在v2PASS之后才允许计算。

17:20:34.960579 UTC生产v2启动，源cb8b919240d0e7823527523d05ad6376ddb103aa0fd571a993ca47b19314d893；17:20:42.191333 UTC完成，终态a6fdb8774aae660db39605bcdc13731b66e1d4f0af95bf93a59962df57d1631c。独立QR审计17:20:42.162065开始，17:20:50.869380 UTC /01:20:50.869380 UTC+8通过全部58模型、11599942wild、116000bootstrap、全部留一及276逻辑槽/四家族校正，证明2fd65c39155cc8812bce11eb256a2d5541027567f626f6a4753579bec7f735cb。最大bootstrap系数复算差1.63e-14以内，是数值精度证据而非生物效应。

全部科学结果已在授权对话审阅，本完整作者状态变化—疗效家族收束，不追加阈值、分母、状态或检验。尚未得到可晋级直接通讯/临床预测的结论。下一核查肺癌队列原生多细胞程序的发现资格，区别于既有跨情境程序投影；不将已反复探索队列的内部留出称完全独立验证。
