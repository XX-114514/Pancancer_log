# 已缓存bulk固定接口：完整资格与临床家族

[事前合同](../methods/cached-bulk-fixed-interface-clinical-family-20261020.md)，[技术收据](../inventories/analysis/2026-10-20/cached_bulk_fixed_interface_v1/technical_metadata_accepted_clinical_family_start_v1.json)。原CM功能代理及其它已闭合家族保持。

## 实际执行

- 物理资格v1：2026-10-19T19:14:58.043667 UTC /次日03:14:58.043667 UTC+8启动，PID157368/startticks325239672；19:15:24.674731 UTC /次日03:15:24.674731 UTC+8完成，终态faa864add9abe72c033994e76c3bc96fe6e8b0493474e983ab2c742f3b6538a0。
- 另解析器全范围审核：19:16:51.802182 UTC /次日03:16:51.802182 UTC+8启动，PID167368/startticks325251047；19:17:37.982643 UTC /次日03:17:37.982643 UTC+8接受，证明ebce4259247d239ed258b3dfb92b312bda32a1047edbef4577c466008c5368a1。全部原文件SHA、完整物理轴、全部临床来源计数另实现相等。只读元数据/来源分母，不解析RNA数值或计算效应。
- 全84临床/39描述问题：19:21:57.454040 UTC /次日03:21:57.454040 UTC+8真实启动，PID182973/startticks325281613。前置合成例子通过，尚待原句柄终态与完整数值/科学审阅，不因观察超时重跑。

运行解释器沿原scanpy Python3.10、NumPy/Pandas/SciPy/statsmodels，BLAS1，无新安装。完整元数据输出/配置仅任务tmpfs，原项目只读；本机只接来源级/技术摘要，不接个体行/表达/模型数组。环境变量PANCANCER_BULK_ROOT、PANCANCER_OLD_BULK_RUN、PANCANCER_BULK_QUALIFICATION、PANCANCER_RAM_OUTPUT由外部私有运行器传入，配置经stdin，不提交私有候选基因表。

```text
python -B scripts/qualify_cached_bulk_fixed_interface_metadata_v1.py < PRIVATE_FIXED_GENE_CONFIG
python -B scripts/audit_cached_bulk_fixed_interface_metadata_v1.py
python -B scripts/run_cached_bulk_fixed_interface_clinical_family_v1.py < PRIVATE_FROZEN_FAMILY_PLAN
```

三步使用各自新任务输出/已接受前步根；生产输出根不覆盖旧版。独立科学审核尚待。原作者候选代码入口不能代替患者真值/随访起点。新临床表及公开推送的两个许可请求均待答，没有绕过；本地技术/科学进度不依赖它们。


## 2026-10-19T19:25:38.973472+00:00 UTC /2026-10-20T03:25:38.973472+08:00 UTC+8 原文本键读取修正

v1在首来源旧分数键集合核对处AssertionError停止（来源行192/176/152），尚无完整终态或任何来源效果接受。只读原始CSV键/当前临床键的集合确认完全对应；Pandas index_col路径即便dtype=str仍推断数值索引，先文本列读取再set_index则逐字保持原键。未展示任何ID。v2只修此读取及已转文本背景列的数值转换，另用csv核每个原始文本键；所有模型/评分函数内容、固定计划、病例、背景及门槛保持，v1源码/失败保留。v2于2026-10-19T19:24:48.348988 UTC /次日03:24:48.348988 UTC+8启动，PID196977/startticks325298702，同计划SHA595950347ed10e52e02826a8ca7dbc4251cc194df06c5298d27ad0e137aa61c0。科学结果仍待原完整运行/审阅，不能以合成通过代替。


## 2026-10-19T19:30:19.290962+00:00 UTC /2026-10-20T03:30:19.290962+08:00 UTC+8 bulk全范围/私有科学图接受

全7既有bulk轴/原来源临床分母、全84临床/39来源描述的实际RNA与完整模型、另一实现全/空模型和校正、私有整批科学及实际两图审阅结束。v1文本键推断故障保留，v2只修读取，原科学计划与核心数学保持。有限问题闭合；下一回到已完成CM整体临床/队列证据，区分总程序关联与给定RNA背景后的增量。正式身份/状态及两个CM版本保持，整体goal active。公开推送及新CSV仍待各自特定许可，无绕过。

## 全家族实际结束与科学交接

v2于2026-10-19T19:25:15.722637 UTC /次日03:25:15.722637 UTC+8结束，终态2350877ef5b82d341d9de06be02a8230612795ab56a72d0492410ef10d3e12e5。全84临床与39来源描述逐位置保留；所有实际可拟合的全/空模型用另一Newton及显式Efron风险集复算，分数另平均秩、旧临床字段/共享背景及原始文本键全恢复。整个家族BH/BY再由另一实现核。共享I/O与既有病例，不称独立采集或新患者身份。

四份来源汇总/证明63194字节已本机私有核SHA及字段白名单，无患者行、基因配置、表达或设计。另完成完整覆盖/来源描述和全部可估计模型森林图，两实际PNG已逐一查看，标签/图例/范围无截断，所有四PNG/SVG及终态513649字节另核。图和科学分级仅私有，不在公共Git。科学全家族接受不等于互作、临床预测或正式身份确认。

原生数组、真实分数、临床设计及拟合数组只RAM；所有新小型科学文件专属tmpfs，原项目只读，未下载新RNA或处理待许可CSV。该有限问题关闭，不改背景/门槛/端点以取有利结果。后续先审已有CM整体证据，明确总程序关联和给定背景后的增量是不同问题；某些背景反映真实状态，不能将条件模型简单解释为去除了全部混杂。原NSCLC来源限定候选与γδ证据仍分层保留，原正式身份/状态/旧CM不改。

公开同步：上一次37763d4的普通公开push被自动审批拒绝于执行前，已发该批具体范围确认请求；本批没有尝试任何绕过或新公开推送。云端最新已核仍6eb0059，本地内容会另commit。原外部临床CSV许可也待定，其它既有数据处理未依赖该许可。

本goal turn为progress：新的完整物理资格与实际有限临床家族改变后续路径；整体目标未完成，也非全局阻塞。建议commit：analysis: qualify cached bulk axes and accept complete fixed interface family。
