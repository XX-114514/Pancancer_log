# Run：CM定义与免疫相关来源内部复现诊断

- 输入基线v20261010/172节点/9CM/1290样本；本任务4源131文库/6情境。不修改原身份、分母列、W/H或权威库。
- 合同/代码commit：v1 4641ec216eea991255410fc3738609c3624c0a51；v2 ecb922f42f2a980b7e3763c3e29c52f84bcb279e。
- 环境：既有scanpy Python3.10.19、NumPy1.26.4、pandas2.3.3、SciPy1.11.4及matplotlib；本地独立设计复算bundled Python3.12、NumPy2.3.5、pandas2.2.3。未改变全局依赖。
- 状态：首轮/敏感性真实计算、选定投影数值复核及全部设计/观测复算完成；独立生物、完整CM复现、通讯轴/临床机制未验证。
- 论文：原文完整35页PDF SHA dea836004cba9806ef6f61c1448715cb9bf2065a6cbf4ae59f3a928945584aad；公开DOI10.1038/s41586-025-09053-4，来源d-nb开放镜像。PDF只用于阅读，未放入public仓库。

## 真实执行时间

| 阶段 | UTC | UTC+8 | 状态 |
|---|---|---|---|
| v1预声明推送核验 |14:48:19 /14:48:27|22:48:19 /22:48:27|普通push及独立SHA一致|
| v1生产 |14:49:44.886755—14:49:52.147499|22:49:44—22:49:52|exit0|
| v2预声明推送核验 |15:02:24 /15:02:32|23:02:24 /23:02:32|普通push及独立SHA一致|
| 独立投影deviance复核 |15:03:12.111565—15:03:14.085238|23:03:12—23:03:14|exit0，36fixtures通过|
| 完整重连/分母敏感性 |15:03:12.116247—15:04:05.311894|23:03:12—23:04:05|exit0，432配置|
| 本地独立设计/观测复算 |准确执行时间未单独记录，15:15取证已完成|23:15取证|2670检查通过|

正式输入文件SHA与v1代码SHA在[摘要](../inventories/analysis/2026-10-12/CM_immune_recovery_v1/summary.json)，v2方案/代码SHA在[新摘要](../inventories/analysis/2026-10-12/CM_immune_recovery_v1/rewiring_v2/summary_v2.json)。结果根为`${PREFLIGHT_RUN_ROOT}/CM_immune_recovery_diagnostic_v1`及`CM_immune_recovery_rewiring_v2`，public只有情境汇总；样本活动/置换向量留private。

## 可复现命令结构

```bash
python cm_immune_recovery_diagnostic_v1.py \
  --frozen-root "${FROZEN_DOWNSTREAM_ROOT}" --review-root "${PREFLIGHT_RUN_ROOT}" \
  --output "${PREFLIGHT_RUN_ROOT}/CM_immune_recovery_diagnostic_v1" --plan "${CODE_ROOT_V1}/plan.json"
python audit_cm_recovery_KL_deviance_v2.py \
  --run "${PREFLIGHT_RUN_ROOT}/CM_immune_recovery_diagnostic_v1"
python cm_immune_rewiring_sensitivity_v2.py \
  --frozen-root "${FROZEN_DOWNSTREAM_ROOT}" --review-root "${PREFLIGHT_RUN_ROOT}" \
  --output "${PREFLIGHT_RUN_ROOT}/CM_immune_recovery_rewiring_v2" --plan "${CODE_ROOT_V2}/plan.json"
python recount_CM_immune_design_from_local_counts_v3.py --workspace "${LOCAL_WORKSPACE}"
```

所有真实运行在新目录，旧失败、图v1和投影检查保留。SSH driver与子进程仅属于本次有界诊断，不宣称已枚举全服务器或重算全量Zarr。

## 结果、失败及边界

首轮54中2过捕获/0可估计旧边共变；51在不满足捕获时仍有非零活动。12投影初次2个SLSQP目标偏高，另有界deviance算法从独立初值复核36全部一致，原H不变。2万置换/2千bootstrap的v2 432配置中120可估计/94唯一问题，14正向BH记录/11唯一问题，BY3全是CM09原发前列腺的重叠读出；不算独立CM或轴。

本地科学依赖缺失改用服务器既有环境；公开论文首次下载超时/续传不支持，另存完整副本；初始患者路径漏family子层已纠正。首轮图文字拥挤，聚合数值不变另绘v2并图像审阅。没有降低v1门槛宣称复现；v2是随后预声明的完整敏感性。

来源使用历史、患者替代块、富集、组织标签和分母影响保留。数值和输入审计不是外部生物复制，原172中未决身份不能升级。详见[全部结果](../inventories/analysis/2026-10-12/CM_immune_recovery_v1/RESULTS.md)与[定义合同](../methods/CM-definition-and-recovery-v1.md)。
