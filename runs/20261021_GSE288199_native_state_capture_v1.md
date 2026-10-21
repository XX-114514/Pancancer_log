# 独立HNSCC全部处理材料的原生RNA状态捕获

## 目标与范围

继续完整生物学/临床互作goal，优先2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）。[上一轮](20261021_GSE288199_external_CM02_measurement_qualification_v1.md)完成独立来源技术资格并改变下一动作，分类progress。当前推进实际表达，不将基因在轴上当状态已经测得。

全部官方声明处理矩阵先做完整header/预算，再逐份RAM流式测量。原22基因/双parent/六观察状态不改，CXCL13仅描述两个parent的检测情况，不新增正式身份。全pre/post、三分选范围、两QC、双parent、四捕获深度，共96来源位置，所有不足/零值/双谱系交叠保留。来源case代号不是已审阅独立患者；临床药物/病理表和空间联合资料均未接受，本轮无新临床、状态对相关或空间检验。

## 效果前读取合同

- 固定原SOFT、全部feature摘要、原生产/独立审核终态。每矩阵绑定自己的feature物理轴；稳定人类Gene Expression、其他Gene Expression及Antibody Capture分开。
- 64KiB压缩header前缀可能含未解析的count字节，只读取/解析MM格式和维度，不提前查看状态。全部header与原feature SHA完成后才接受预算；total压缩≤16GiB/单份≤4GiB/列≤100万，实际RAM估算加3GiB余量，正式流并发2，无矩阵落盘。
- 全compressed字节数/LastModified/header/prefix固定，完整流计算SHA并验证gzip CRC、声明nnz及非负整数、唯一单调坐标。兼容feature优先或cell优先两种次序；不是按效果调整。每chunk的marker、三类UMI、线粒体和检测人类基因数由原坐标加法与另一稀疏矩阵实现全部复算。
- 所有列及human RNA≥500 UMI/≥250基因/线粒体≤20%的RNA重分析两口径并报。仅借鉴原论文RNA过滤，不声称复现作者双细胞/TCR/聚类细胞集合，也不把公开材料直接当论文最终子集。
- 原marker逐细胞vector/scalar逻辑一致，R/CD8双阳性交叠保持未决；未达到QC的列不是阴性细胞。每文库/QC在人类RNA值三分位分捕获深度，ties一起移动；不强行等量。
- 来源机会≥5 coded单位才给池化观察数；每parent≥20观察且≥5支持单位才给病例平衡分位。不足5分布抑制。全96、原parent/marker/质量和可供下一状态对问题的覆盖全部保留，无P/FDR。
- per-cell、coded-case/phase、表达及中间数组只在RAM，输出仅来源聚合/技术证明，无barcode/VDJ/新临床读取。原CM/W/S/H、权威标签及所有旧家族保持。

## 实际执行

```text
python -B scripts/measure_GSE288199_native_state_capture_v1.py
```

环境：既有服务器scanpy/Python、NumPy/SciPy，线程1。新版本控制目录由PANCANCER_RAM_OUTPUT指定；PANCANCER_EXTERNAL_MEASUREMENT_OUTPUT/PANCANCER_EXTERNAL_FEATURE_AUDIT_OUTPUT固定前置终态。网络仅direct官方HTTPS，有界3次同URL重试，无继承代理/备用镜像。

主脚本SHA7b0ff86d68a5868d4213ad27a88e6e25b239a9837463b1cab6669b162e88b1c7；reader SHA c62ef83dbd9810e83a409e236cee7243f2189b113398704c461bd59ae02f678c。2026-10-21T03:30:11.320591+00:00（UTC+8 11:30:11.320591）实际启动，真实PID/startticks留私有。

03:30:29.586746 UTC取证原PID/startticks仍活跃：合成MM全读取/mmread、两坐标次序/小chunk跨行、显式零、人类/抗体/线粒体边界、重复坐标拒绝已通过；全部矩阵header census进行中。完整预算、RNA结果与接受仍待，不能由启动/部分进度声称完成。

## 验收及后续

实际全60流数值通过后，再核来源完整96位置、QC/深度/parent/命名空间分母与不足分布、marker捕获和双端共同支持。执行与数值检查、来源科学解释分别验收；source summaries并不等于独立患者身份或稳定状态组合验证。病例平衡分位若只核界限，不冒称第二次逐case数组重放。

后续先依真实覆盖冻结源内状态关系及采样/质量背景，临床分支待原Table S1精确桥，空间分支待坐标/共面板/病例。TCF7/IL7R模式不直接改名为作者TEM/TRM。当前尚无新的RNA结果或临床/空间增量，整体goal active；原受拒公开push/新随访表未获特定许可，不重试。


## 2026-10-21T04:02:31.879243+00:00 UTC /2026-10-21T12:02:31.879243+08:00 UTC+8 用户确认临床附件待补与测量继续

用户明确回复：先继续，记录GSE288199 TableS1缺失，稍后补充。临床表状态设为待用户补充，不重复请求，不从图或来源编号补出病例药物/病理值；不因此暂停结局盲RNA测量或整体goal。NeoTRIP许可问题仍未答，不访问受限文件。

全矩阵header/预算于03:34:11.743244 UTC通过；直连服务器完整读取较慢。相同官方1MiB片段服务器20.57s、本机2.25s，prefix SHA完全相同；SSH Unix临时通道保持官方TLS但26.06s，无速度收益，未切换科学作业，已请求关闭本任务通道，不更改既有代理或TCP端口。

既有本机NumPy2.4.1/Python3.13.9与R4.5.2/Matrix1.7.4可用，未安装。最初R输出别名fixture生成仅4byte握手文件，未触真实表达，原失败代码/fixture保留；显式/dev标准管道在默认沙箱受阻，使用与实际任务相同执行权限后原代码合成通过。保持原计数/坐标/CRC门槛，独立reference改用每chunkR Matrix binary pipes，全部同scope、不输出个体数组。

03:58:21.694252 UTC一份完整官方真实矩阵NumPy/R Matrix所有marker/命名空间/检测/线粒体逐块精确一致，完整字节/CRC/LastModified/prefix通过，约5.16s；未审阅状态或疗效，原raw传输/header/parser方法AST与原reader一致。0400UTC按原PID/startticks受控停止旧v1，保留全部历史；当时5/60已校验，未接受完整结果。此次停止基于验证过的性能替代，不因观察超时判失败。旧作业停止证明已转存，原完整header摘要固定。

04:01:12.655285 UTC本机private RAM全量v3真实启动，先重放所有60 header并要求与原摘要严格相同，再计数；原96机会、六状态、QC和全部材料不变。并发2，R Matrix逐块reference与全部cell scalar规则；数据和病例数组仅RAM，最后回传来源汇总。新任务实际终态仍待。现有临床/空间旧家族、CM02/W/S/H和注释保持，整体goal active，原受拒公开push/新随访表不重试。

## 实施版本对照

|组件|实际角色|验收边界|
|---|---|---|
|main v1 + reader v1|原服务器SciPy作业，已有完整预算后受控停止|部分材料，不接受完整科学结果|
|main v2 + reader v2|仅准备SSH Unix relay替代|速度探针无增益，未运行科学全量|
|R reference v1 + reader v3|本机合成fixture失败历史|标准管道别名问题，无真实表达|
|main v3 + reader v4 + R reference v2|当前实际本机全量版本|既有NumPy/R Matrix；全量终态与来源审核待|
|source auditor v1|独立来源汇总守恒审计|完整96位置/边界，不冒称独立逐case分位重放|

原reader v1与实际v4的opener/feature_contract/header/census/chunks/Tracked/read_full函数AST逐项相同；独立数值reference由SciPy换成现有R Matrix，需单独合成、真实完整文件及全部chunk一致，不能仅由传输函数相同推定计数正确。


## 2026-10-21T04:12:38.248152+00:00 UTC /2026-10-21T12:12:38.248152+08:00 UTC+8 原生RNA全材料与完整来源测量接受

本机实际v3于2026-10-21T04:09:44.092091+00:00 UTC结束，全部官方材料完整流/gzipCRC/长度/LastModified/header/nnz与全部chunk跨语言数值核对通过；实际代码副本和终态清单逐SHA相同。source auditor v1全部96位置、父类/交叠/marker/深度分区与QC子集和小分布规则通过，原生RNA质量与全部预设来源范围已私有审阅。逐case分位仅核支持与边界，不冒称独立逐case数组重放；本轮没有新的状态对/临床/空间效果检验或原CM投影，不新增正式身份。

原v1服务器慢速作业按同PID/startticks受控停止，旧部分产物保留；SSH Unix relay无速度增益、科学main v2未启动。实际本机main v3/reader v4/R reference v2完成，所有版本/fixture失败历史保留。已确认本任务临时relay及实际RNA worker均结束，原任务和代理配置未改。

用户确认临床Table S1稍后补充，当前为deferred，不重复请求。旧mmc2输入属于原Cha/GSE286827，并非当前Li/GSE288199；临床材料/药物/病理桥和联合空间仍未接受。下一先根据完整测量支持固定全部R×E来源状态问题与捕获/两端丰度背景，不能由可测直接声称稳定互作。具体观察数量/状态比例/捕获结果留私有；公共仅技术轨迹、通用代码和SHA。

自动审批新拒绝将本轮非公开来源汇总/科学审阅写回具体Project_v5目录，理由是该批派生结果向该目的地出站缺明确授权；拒绝命令未执行，未重试或绕过，已询问只含来源汇总/通用代码/技术证明的具体目录许可。本机完成验收及日志整理继续。此服务器回传问题与原公共push/NSCLC随访表拒绝分别记录；原两项仍无特定新许可，不重试。

技术验收：[receipt](../inventories/analysis/2026-10-21/GSE288199_native_capture_v1/technical_complete_acceptance.json)，SHA f11767f9e3f7ae6f6cbc482ce3d4b089f4007d80b5daef36b479f38a7e8ec982。执行与数值和完整测量范围接受；状态组合、空间组织与同终点临床增量尚未完成，整体goal active。本轮分类progress，建议commit: Measure independent HNSCC native RNA state capture and defer clinical bridge。
