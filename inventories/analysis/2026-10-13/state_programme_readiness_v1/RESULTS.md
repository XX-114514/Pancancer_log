# 细胞状态准备诊断：真实结果与研究方向

## 运行与验证范围

2026-10-13 UTC12:46:42.997136—12:47:30.981673生产，worker到12:47:31.101092 exit0；12:47:46核实worker/child均退出。11来源、52source-major组合、2028已审major-library单位，固定181参考程序×52=9412记录和16290参考程序对。所有原输入SHA、轴和完整键检查；52组合固定Mac_IFNG程序按直接标量逐单位均值/检测循环复算通过。本机取回后核所有文件SHA、全部9412键/分母/汇总。**没有独立重演每个程序分数，没有新p或单细胞状态标签。**

[全表](all_source_major_program_readiness.tsv.gz)、[字典统计](all181_reference_program_catalogue.tsv)、[全部重叠](all16290_reference_overlap.tsv.gz)、[52标量检查](direct_scalar_diagnostic_checks.tsv)、[执行](cell_state_programme_readiness_20261013_v1.execution.json)、[登记/刷新](acceptance_receipt.json)。自己的诊断证据登记完成/0冲突；原CM证据v17、身份/9CM/权威库未改。

## 三项实际发现

| 诊断 | 结果 | 含义及边界 |
|---|---:|---|
| 任意预声明覆盖门槛通过 |8217/9412，87.3%|支持多数reference继续检查；门槛是准备诊断，不代表这些状态存在或fine身份适用 |
| 均值程序与非LR检测基因复杂度abs Spearman≥0.5 |6386/9412，67.8%|当前mean-log程序与复杂度高度相关；可能生物或技术，不因相关就回归删除或宣称伪影；没有患者独立性或p推断 |
| 全局剔除LR损失至少50%已匹配签名基因 |1036/9412，11.0%|全资源剔除会削弱若干真实状态读出；需保留完整状态字典，通讯关联时做候选接口/共同基因针对性诊断 |

这些数目跨source-major-reference记录，不是9412独立状态或统计假设。参考细胞类型程序全部检查覆盖，但CD4/CD8/Treg适用性需身份审阅，不能直接套 broad T/NK。

### 全局LR删除具体损失

| 官方程序 |参考基因|在全LR字典内|全局剔除后剩余|
|---|---:|---:|---:|
|Mac_IFNG_response|21|16|5|
|all_type-II-ifn-response|17|15|2|
|CD8-T_tumor-reactive-like_UP|8|7|1|
|T_tcr-activation|66|35|31|
|CD8-T_terminal-exhaustion|13|5|8|
|TNK_cytotoxicity-effectors|14|6|8|
|all_type-I-ifn-response|17|1|16|

这是固定外部字典与原4624LR字典的重叠，不是按本数据p选择的基因。7个已有Macrophage来源Mac_IFNG均匹配21基因，剔除后确实仅5；没有用这些结果宣布IFNG轴成立。参考typeI/typeII仅共享1基因，与Hallmark gamma/alpha共享73的字典不同；低基因重叠也不保证生物刺激特异性。

## 最新通讯结果如何指导状态问题

此前原CM精度阶段确实有178q记录/147区间支持，主要CRC，不能统称所有阶段失败。巨噬–T从8主BH到共同程序控制后0，提示共同状态/自身RNA值得检查，不单独证明技术伪影；亚型控制也可能控制生物中介。

并行10/13外部接收程序阶段已由主线登记[147968模型结果](../mac_T_receiver_programs_v1/RESULTS.md)。IFNG→巨噬gamma主β.329/BH.918，亚型β.131/BH1；α/TNF也有名义关联。其他3条BH线索BY均不通过、亚型敏感。它们支持先分清共享炎症、身份混合与状态程序，不支持以某个最有利变体发布轴。完整IFNGR1/2同细胞检查由并行支线推进，本任务不重复或覆盖。

## 下一步可执行研究

**方法**：Spectra为优先候选，cNMF作为无先验对照，固定模块后UCell类rank评分。当前三者未安装/训练，本实验只是既有major均值缓存诊断，不能冒称单细胞factor或UCell。

**首批对象**：T/TNK和髓系；从现有已审来源冻结源/患者/文库均衡的单细胞试点与留出范围，保留各数字subcluster。程序至少分开细胞毒、TCR/活化、耗竭相关、IFN、增殖、炎症、代谢和应激。重点检查同一原编号内的连续变化、来源/复杂度和共同应激；有多基因/稳定性/投影证据才审阅为state候选。

**标签**：原major/subcluster不动，另存programme scores、多值候选标签、置信度、缺失/待审状态、基因覆盖与模型SHA。例子`T_cell_CCL5_38 | cytotoxicity_high; IFN_response_high`仅说明格式，不是本轮赋值；不以最高分强制单状态、不按每库四分位制造固定比例。

**回到生态/通讯**：用样本级、身份内程序及其协同描述原9CM；另建state-programme网络用SP编号。总体关系、身份混合、两端共同状态和接收自身表达分步建模；状态可能是混杂、中介或结果。检验固定候选LR时删除其直接重叠基因，并列全局LR剔除敏感性，不根据p挑最优版本。状态/生态可来源限定；分选、患者、癌种混杂及原发现重叠持续保留。

详细方法和全部阶段复盘见[方案](../../../../methods/cell-state-programmes-and-communication-review-v1.md)。原历史阴性、来源限定阳性、技术失败及新模型敏感线索全部保留。尚无新正式状态、独立通讯机制或耐药结论。
