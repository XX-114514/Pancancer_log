# 原生UCell真实评分与独立核查 v2

## 已执行

固定98signature/65严格稳定候选（跨K重复不是独立状态），03:42:30—03:52:16 UTC worker/child exit0；真实审计03:59:04—03:59:44 exit0。5source全部240533cells，每源完整原生基因轴，CPU平均并列秩，主maxrank1500、1000/3000敏感性及90CN的TPM signature敏感性；不改变type/node，也不按每库百分位强造状态。

[固定字典](../frozen_state_programme_dictionary_v2/frozen_programme_dictionary.json)、[490覆盖](all490_native_gene_coverage.tsv)、[source-major分布](all_source_major_programme_distribution.tsv)、[原编号内分布](all_source_node_programme_distribution.tsv)、[审计](independent_audit.json)。全部490 gene-axis记录50/50genes实测；模型scope外保留NaN，实测零就是0。完整原cell IDs、gene检测/fraction、全部score和91599个sample-node-programme均值留server，不能把这些unit称独立患者。

## 独立核查的范围与限制

- 全240533cells/98programmes的47144468个检测/fraction值独立稀疏incidence复算，所有scope/coverage mask、cell IDs及原identity metadata精确。
- 200个固定随机真实cells（每source/major20）手工分组排序处理并列秩，三rank与TPM共38700scalar，最大差4.47e-8；**不是所有cell score再次全算**。
- 全91599私有sample-node-programme均值/valid及informative计数独立incidence聚合复算，最大均值差2.02e-7。
- 公开分布及相关字段未逐项全独立重演，生物/机制和正式state命名尚未通过，不能由numeric PASS直接升级。

## 实际发现

330个main稳定candidate×source×允许major记录中，92个（27.9%）与完整RNA abs rho≥.5，median abs rho.328。它们不是330独立状态；与此前mean-log准备诊断的变量/分母不同，不直接宣称降低多少技术偏差。

maxrank1000对1500的rho median.952、min.656，说明部分programme对cap仍敏感；TPM signature对正coef signature的rho median.512，两种gene选法不能当相同读出，更不能追最有利LR p选读出。

main programme的>=5原生gene检测信息支持比例差异很大：Mac median约24.5%，T约97.4%，有programme/来源近零。信息支持阈值不是阳性状态或celltype判定。外来谱系RNA、基础ribosomal/mitochondrial、source/原node依赖必须单列，再讨论是否适合状态命名。

GSE183904是此前已使用来源的内部投影，不称独立验证。当前pilot主要提供连续、多程序候选状态层，尚无新的正式type、耐药或免疫抑制结论。

## 完整目标

[全部172/1290预检](../all_CM_state_input_preflight_v1/RESULTS.md)已覆盖52sources/4532933cells。补其余18major的发现/原生评分，分清identity-associated与activity-associated，再把状态映射原9CM；固定状态LR/接收family，完成状态后DIALOGUE及新空间。原172/9CM/Zarr及历史阴性/来源限定阳性保持。
