# 20261018 priority immune CM bulk v1

- 范围：2026-10-14免疫图谱CM v1、T_whole_TME_primary_tumor、CM02/CM04固定功能代理。
- 状态：312预设主要槽全部登记；162可拟合设计全部独立数值重拟与完整家族校正接受；其余150不可估计保留。4幅v2图已全部实际查看。临床结论不自动晋级。
- 代码：计划先行commit2039bcc；生产d3b6588、38a3e5c、9379f89；最终补充与渲染代码SHA见技术收据。远程部署快照和失败版本保留。
- 环境：scanpy Python3.10.19 /NumPy1.26.4 /pandas2.3.3 /SciPy1.11.4 /statsmodels0.14.6；未安装新环境。
- 输入：`${PROJECT_ROOT}/Project_v5/bulk_cohorts_20260904`既有TCGA与5ICI；全native充分统计只读；作者公开SU2C-MARK v3处理数据。
- 主根：`${PROJECT_ROOT}/Project_v5/runs/20261018_priority_immune_CM_bulk_v1`。

## 实际时间线（UTC；括号UTC+8）

- 08:55:38（16:55:38）v1启动，08:56:05因Rose非精确临床key断言退出；不接受不完整结果。
- 08:57:49（16:57:49）v2按精确交集修复，08:58:53生产完成；08:59:14全304槽/154设计独立接受。
- 09:02:30（17:02:30）本机公开三成员CRC/SHA验收；服务器直连Range失败及全TPM慢传输部分文件均保留不用。
- 09:08:29（17:08:29）SU2C外部链启动，09:08:31全8槽完成，09:08:32全部8设计独立接受。
- 09:08:45/46（17:08:45/46）原来源与外部描述性诊断完成；09:11:13（17:11:13）外部8临床控制槽完成。
- 09:15:06（17:15:06；Melbourne20:15:06）完整v2四图实际复核登记，原312槽效果表与v1字节一致。更精确事件ISO时间见各收据。

## 可复现执行结构

```text
explore_priority_immune_CM_bulk_v1.py --base ${PROJECT_ROOT}/Project_v5/bulk_cohorts_20260904 --atlas ${ACCEPTED_ATLAS_ROOT} --stats ${ACCEPTED_NATIVE_STATS_ROOT} --plan FROZEN_PLAN --output ${NEW_RUN_ROOT}/analysis_v2
# ACCEPTED_ATLAS_ROOT需含tables/subcluster_catalog_refined_private.tsv；STATS需含genes_private.tsv、sum_CP10K.npy和measured_counts.npy。
# NEW_RUN_ROOT必须是新版本目录，原已接受目录拒绝覆盖；其余审核/渲染命令相应使用新根。
audit_priority_immune_CM_bulk_v1.py --run ${RUN_ROOT}/analysis_v2 --output ${RUN_ROOT}/audit_v2
explore_SU2C_priority_immune_CM_v1.py --source ${RUN_ROOT}/SU2C_MARK_panel_ready_v1 --plan FROZEN_PLAN --external-plan EXTERNAL_PLAN --output ${RUN_ROOT}/SU2C_analysis_v1
audit_priority_immune_CM_bulk_v1.py --run ${RUN_ROOT}/SU2C_analysis_v1 --output ${RUN_ROOT}/SU2C_audit_v1 --expected-slots 8
diagnose_priority_immune_CM_bulk_v1.py --run ANALYSIS --audit AUDIT --output DIAGNOSTICS --mode existing
# 外部诊断相同参数改为--mode external。
describe_SU2C_clinical_controls_v1.py --run ${RUN_ROOT}/SU2C_analysis_v1 --audit ${RUN_ROOT}/SU2C_audit_v1 --output ${RUN_ROOT}/SU2C_clinical_controls_v1
render_priority_immune_CM_bulk_v1.py --run ${RUN_ROOT}/analysis_v2 --audit ${RUN_ROOT}/audit_v2 --external ${RUN_ROOT}/SU2C_analysis_v1 --external-audit ${RUN_ROOT}/SU2C_audit_v1 --output ${RUN_ROOT}/research_results_v2
```

## 收据和产物

公开技术入口：[最终技术验收](../inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/final_technical_acceptance_public.json)、[304槽独立审核](../inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/independent_main304_audit_public.json)、[8槽独立审核](../inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/independent_external8_audit_public.json)、[公开下载](../inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/SU2C_public_members_download_receipt.json)。原生产terminal中的PENDING保留，后续独立审核是当前终态。

私有/授权产物：analysis_v2、audit_v2、SU2C_analysis_v1、SU2C_audit_v1、diagnostics_v1、SU2C_diagnostics_v1、SU2C_clinical_controls_v1、research_results_v2/RESEARCH_INTERPRETATION_CN.html与index.html及4组PDF/SVG/400dpiPNG。全部原失败/partial文件标为不供分析，不删除。

科学数值、患者分数与设计不进public Git。数值审核不是独立采集或治疗验证；部分来源缺固定基因、时点或临床控制，全部保留。全33癌种计划中仅32癌种有01材料。样本/疾病/终点不同，不能拼接正向结果宣称泛癌预测；跨原生与bulk研究患者完全不重叠仍未确证。

方法与下一阶段：[合同](../methods/priority-immune-CM-TCGA-ICI-exploration-20261018.md)。任务细节、错误和同步记录见[日志](../logs/codex/2026/2026-10-18_priority-immune-CM-TCGA-ICI-exploration.md)。
