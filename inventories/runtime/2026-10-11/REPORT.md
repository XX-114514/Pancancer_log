# 下载与Project_v3运行快照

采集时间（UTC）：**2026-10-11T03:22:53.508526+00:00**。这是一次实测快照，不是GitHub实时监控。

## 新公共数据收集

- 已验证下载/源矩阵库存：**5,341,204源细胞**。
- QC分片：**3份，608,684细胞**，仍待最终来源验收。
- 正式验收新增：**0**；目标为排重、QC后450万—550万。
- 用户已授权下载与逐来源前处理；跨研究整合等待单独决定。

| 队列 | 源细胞 | 已完成的处理 |
|---|---:|---|
| CAR-T治疗产品 | 417,167 | QC保留416,800；待最终来源验收 |
| 卵巢癌SPECTRUM | 927,205 | 原始计数Zarr读回验证完成；待QC/来源验收；完整计数对基线精确匹配0个；不能排除不同基因覆盖下的重叠 |
| 大B细胞淋巴瘤 | 598,289 | 原始计数Zarr读回验证完成；待QC/来源验收；完整计数对基线精确匹配0个；不能排除不同基因覆盖下的重叠 |
| 儿童高级别胶质瘤 | 401,253 | 原始计数Zarr读回验证完成；待QC/来源验收；完整计数对基线精确匹配0个；不能排除不同基因覆盖下的重叠 |
| 转移性乳腺癌 | 394,534 | 原始计数Zarr读回验证完成；待QC/来源验收；完整计数对基线精确匹配0个；不能排除不同基因覆盖下的重叠 |
| 未治疗TNBC | 427,823 | 已下载；来源/文库/QC策略审计中 |
| CTCL图谱 | 1,533,956 | 已下载；来源/文库/QC策略审计中 |
| 头颈癌免疫治疗 | 416,211 | 已下载；来源/文库/QC策略审计中 |
| 肾癌治疗 | 145,316 | QC保留137,728；待最终来源验收 |
| 胶质瘤疫苗 | 55,570 | QC保留54,156；待最终来源验收 |
| 肾髓质癌 | 23,880 | 原始计数Zarr读回验证完成；待QC/来源验收 |

### 第二批文件传输

| 队列 | 已校验文件 | 未完成文件进度 |
|---|---:|---|
| GSE315928 | 1/2 | GSE315928_RAW.tar：12%，0B/s（日志超过10分钟未更新，需检查作业） |
| GSE294273 | 0/1 | GSE294273_RAW.tar：68%，0B/s |
| GSE256326 | 4/4 | 文件校验完成，计数与临床验收另计 |
| GSE236581 | 3/4 | GSE236581_counts.mtx.gz：10%，0B/s |
| GSE246613 | 0/2 | GSE246613_PembroRT_immune_R100_final.h5ad.gz：等待/未观测到传输<br>GSE246613_PembroRT_non_immune_cells.h5ad.gz：等待/未观测到传输 |
| GSE332846 | 0/1 | GSE332846_RAW.tar：等待/未观测到传输 |

Zenodo 19586163：36%，24KiB/s；官方校验完成标记：False。传输进度不作为完成证明。

### 当前审计限制

- TNBC的7个orig.ident是分组，不能作为原始文库；472种barcode前缀也尚未获得作者文库映射。
- CTCL的13个来源需追溯原研究别名、文库及基因覆盖；570,282个非primary标记细胞处于来源隔离审计。
- CAR-T为治疗产品，原位肿瘤生态资格单列。临床缺失不从表达矩阵推断。
- CRC元数据有975,275条细胞记录，但计数文件未完整校验时不计入上面的完成下载规模。

## Project_v3当前动态

后续独立审计和最终交接登记优先于生产者早期manifest的pending检查点；保留早期状态用于追溯。

源状态表生成时间：`2026-10-11T03:03:12.990848+00:00`；本次重新读取了当前登记和轻量审计文件，并保存SHA-256。报告中的PASS是服务器已有审计的记录，本次没有重算数百万细胞。

- 核心范围：4,676,787细胞；登记细身份非空926,810。
- 最终科学身份发布：`NOT_READY_FOR_FINAL_IDENTITY_RELEASE`；GSE299623仍有40群/54,099细胞未获最终身份决定。
- 表达注释与新下游交接：`COMPLETE_USER_REQUESTED_EXPRESSION_ANNOTATION_ZARR_UMAP_COVARNET_LIANA_AND_DOWNSTREAM_REVIEW`。描述性单基因命名保留既有身份，不自动成为新细胞亚型。

| 运行/证据 | 已记录状态或规模 |
|---|---|
| 20261010_uniform_expression_annotation_v1/input_contract.json | 已登记 |
| 20261010_uniform_expression_annotation_v1/run_register.json | COMPLETE_USER_REQUESTED_EXPRESSION_ANNOTATION_ZARR_UMAP_COVARNET_LIANA_AND_DOWNSTREAM_REVIEW |
| 20261010_uniform_expression_annotation_v1/manifest.json | COMPLETE_ALL_GSE_EXPRESSION_PENDING_INDEPENDENT_AUDIT |
| 20261010_uniform_expression_annotation_v1/independent_audit.json | PASS_FULL_SCOPE_EXPRESSION_READBACK_AND_INDEPENDENT_PILOT |
| 20261010_uniform_expression_annotation_v1/final_goal_scope_audit.json | PASS_FULL_USER_GOAL_ARTIFACT_SCOPE_CANONICAL_FULLROW_AND_ALL_DOWNSTREAM_GATES_PENDING_FINAL_HANDOFF_REGISTRATION |
| 20261010_uniform_expression_annotation_v1/final_goal_handoff_registration.json | COMPLETE_USER_REQUESTED_EXPRESSION_ANNOTATION_ZARR_UMAP_COVARNET_LIANA_AND_DOWNSTREAM_REVIEW |
| 20261010_expression_annotation_downstream_v1/run_contract.json | PREPARED_NEW_LABELS_EXACT_ROSTER_SOURCE_HOLDS_PENDING_INDEPENDENT_AUDIT |
| 20261010_expression_annotation_downstream_v1/run_register.json | COMPLETE_USER_REQUESTED_EXPRESSION_ANNOTATION_ZARR_UMAP_COVARNET_LIANA_AND_DOWNSTREAM_REVIEW |
| 20261010_expression_annotation_downstream_v1/input_independent_audit.json | PASS_NEW_LABEL_DOWNSTREAM_INPUTS_EXACT_SCOPE_COUNTS_AND_SOURCE_HOLDS；included_rows=4534433；samples_min100=1290；nodes=172 |
| 20261010_expression_annotation_downstream_v1/covarnet_audit_v1/manifest.json | PASS_COVARNET_EXACT_FREQUENCY_NORMALIZATION_NEW_LABEL_FACTOR_DIMENSIONS_AND_SOURCE_DIAGNOSTICS；samples=1290；nodes=172；K=9 |
| 20261010_expression_annotation_downstream_v1/liana/independent_audit.json | PASS_LIANA_EXACT_SCOPE_NEW_LABELS_NATIVE_GENES_INPUT_SIGNATURES_NORMALIZATION_AND_OUTPUT_READBACK；samples=1312 |
| 20261010_expression_annotation_downstream_v1/cm_liana_summary_v1/independent_audit.json | PASS_NEW_CM_LIANA_FULL_STREAM_PARITY_NEW_WEIGHTS_NATIVE_COMPLEX_DENOMINATORS_AND_SUMMARIES；all_LIANA_rows=31150424；within_CM_rows=13729410 |
| 20261010_expression_annotation_downstream_v1/covarnet_network_corrected_v1/network_independent_audit.json | PASS_CORRECTED_COVARNET_ALL_GLOBAL_AND_CM_EDGES_METADATA_WEIGHTS_NO_MODEL_RECOMPUTATION |
| 20261010_expression_annotation_downstream_v1/downstream_review_v2/table_and_render_audit.json | PASS_DOWNSTREAM_REVIEW_FULL_TABLES_SOURCE_CONFIDENCE_EDGE_JOIN_LR_SELECTION_AND_RENDER_FILES_PENDING_AGENT_PIXEL_REVIEW；nodes=172 |
| 20261010_expression_annotation_downstream_v1/downstream_review_v2/agent_visual_review.json | PASS_AGENT_PIXEL_REVIEW_ALL_FIVE_DOWNSTREAM_FIGURES_AND_NINE_CM_NETWORK_COMMUNICATION_PANELS |
| 20261010_CM_source_replication_v1/run_register.json | COMPLETE_FIRST_ROUND_SOURCE_REPLICATION_AND_SAMPLE_LR_AUDITED |
| 20261010_CM_source_replication_v1/results_summary.json | COMPLETE_FIRST_ROUND_SOURCE_REPLICATION_AND_SAMPLE_LR_AUDITED；samples=1290；reviewed_samples=121；reviewed_patient_units=45；primary_model_BH_significant_strata=0 |
| 20261010_CM_source_replication_v1/execution_final.json | COMPLETE_AUDITED |
| 20261010_CM_source_replication_v1/audits/independent_audit.json | PASS_FROZEN_CM_SOURCE_REPLICATION_AND_SAMPLE_LR_INDEPENDENT_AUDIT |
| 20261010_CM_source_replication_v1/audits/meta_readback.json | PASS_ALL_COMPOSITION_META_SCALARS_AND_EMPTY_COMMUNICATION_META_SCHEMA |
| 20261010_CM_source_replication_v1/audits/visual_review.json | PASS_AGENT_PIXEL_REVIEW |

本次观察到Project_v3科学/流水线进程：**0**。代理、日志刷新等基础进程不计为科学计算。进程存在不等于结果验收；观察范围为当前用户可读进程与所登记路径。

### 下一步

1. 补齐CM候选来源的患者、癌种、组织、分选和终点映射；六个可估计主模型分层没有BH显著结果，不升级为已验证机制。
2. 复核K与来源/分母敏感性；9个CM是当前模型，不声称已证明全局最优。
3. 新临床队列完成来源、文库、QC和临床端点验收后，冻结独立验证集；不在验证集中重新选择marker或重拟合发现模型。
4. GSE299623身份审阅单独闭环；表达状态、细胞身份、CNV恶性和生态/临床验证分别登记。

完整机器快照与文件哈希：[runtime_snapshot_v2.json](runtime_snapshot_v2.json)。
