# NSCLC两端协调及完整程序描述：实际执行与验收

## 时间和环境

2026-10-16 UTC / UTC+8：预检10:07:02—10:07:08 /18:07:02—18:07:08；效果10:10:31—10:17:06 /18:10:31—18:17:06；独立效果审计10:15:59—10:18:19 /18:15:59—18:18:19。程序参考描述10:12:45—10:12:54 /18:12:45—18:12:54，独立描述审计10:18:31—10:19:01 /18:18:31—18:19:01。上述是原收据/终态事件，本记录在10:30 UTC之后补录，非早期已发布。10:30后重新只读核全部终态。

scanpy Python3.10.19、NumPy1.26.4、pandas2.3.3、SciPy1.11.4、Statsmodels0.14.6；BLAS线程固定1。逻辑根`${RUN_ROOT}`为既有Project_v5 preflight，程序由独立会话运行，收据在`${RUN_ROOT}/HNSCC_receiver_followup_job_control_20261016_v2/<script_stem>/launch_receipt.json`。基于6c13b1c后的原生测量/临床验收，实际源码SHA为执行依据，不以之后提交假称运行前commit。

## 方法和实际范围

[冻结方法](../methods/NSCLC-two-end-coordination-and-redundancy-v1.md)：全部572计划位置、360可估计；其余保留覆盖或杠杆原因。24签名全276权重对、1106同亚群描述、2400参考槽均完成。科学效应、临床交集、P/BH及基因内容仅服务器私有目录和获准对话审阅；公共仓库只保存代码/方法/执行轨迹和指纹。

主要命令由服务器环境提供根目录，按依赖顺序运行：

```bash
python scripts/prepare_NSCLC_two_end_coordination_v1.py
python scripts/test_NSCLC_two_end_coordination_v1.py
python scripts/audit_NSCLC_two_end_coordination_v1.py
python scripts/review_NSCLC_signature_redundancy_function_v1.py
python scripts/audit_NSCLC_signature_redundancy_function_v1.py
```

后三个作业的实际启动有依赖等待，不能将启动时间当完整输入已就绪。原始生产终态的audit_pending保留，后续独立审计覆盖该阶段状态。

## 指纹与验收

| 证据 | SHA256 |
|---|---|
| 预检代码 | b224fef2ff4146476631ed46362954322d10392a8846e4f3956b4a8e33172572 |
| 效果代码 | a6dea51b1621b8c8e545e6be1db89ad70fc45983aab3a5aca1e768ab0166b464 |
| 效果审计代码 | a20015c1d55e020bdffe544d8482ca4cd7d65963f8c6b59c0167f3c1a9732f9f |
| 描述代码 | 4b1d3fb89c480843588004076ab519fe6d6f88ca7ea6fcc7d5f8e1128649efb5 |
| 描述审计代码 | 00eba5c36c4574d4a94344c60a748e3f963aac1ee72c9e6d4f19d5bca12cf9be |
| 预冻结合同 | 7a0ed483a47dfb35873d15f3d2e424888a46b8c67712d0bfeffbf3c57bc2bb41 |
| 全572范围 | a5694dd31a52fde457a7b5f743f6b5459ea21ca7341607686a6e29ad6d8f75ae |
| 全572结果 | d05bb3d85a5985f01ba539669f6f2e2b348762cd7787c01f92df5cb4474a21ac |
| 效果独立审计 | f77312a22226752ba360afb20a6e0f45ffa0a9397bc99d4a2cf89f039d0261bc |
| 描述生产终态 | 3f1d3ff10879c83a9bee07a6e573cd66495476c936cc9ef1b8779779b43610ca |
| 描述独立审计 | 202668a01cb9685a3e6a6afade1ccc535da8d2c3f39ebdc4eb9872675aa1c505 |

对应私有目录：`${RUN_ROOT}/NSCLC_two_end_coordination_preflight_20261016_v1`、`NSCLC_two_end_coordination_results_20261016_v1`、`NSCLC_signature_redundancy_function_20261016_v1`。后两者的`private_independent_audit.json`分别PASS完整范围/点/HC3/敏感性/BH/留一/三例重放，以及全276标量余弦/24聚类规则/96权重描述/1106分数/2400超几何记录；不称全部模型随机流均独立重放。

## 解释边界及下一

当前为来源内固定程序协调和来源内临床线索，不是直接LR、预测、正式细胞状态或机制确认。病理条件化可能有collider偏倚，不控制病理的敏感性单列。参考富集不替代去基因独立响应。原MCP编号不互换，原MIF/TF/参考不支持证据保留。全部目标保持active。

下一自动先验排序、冻结原生表达/受体覆盖与后续完整关联家族；不重复原协调或临床家族。只读汇总一度因重复列名JSON失败，唯一列查询恢复；科学生产未重启。一次广泛web加远程输出过长，后续改为有界分批读取，未改变结果。没有新下载、删除或覆盖并发任务。

本轮另将完整临床—协调候选交集、原反例、条件化/来源限制和文献增量判断写入服务器私有 `private_scientific_review_20261016_v1.md`，SHA `5cbe3fdd50b9d6b08ad457d6049eb70268fd29498562479ac289191da280430f`。该报告包含科学值，未复制本机或公共仓库；交集不是联合FDR，后续输入未据此缩减。
