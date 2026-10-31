# 2026-10-31 数据库存核查：WEHI目标端复核与缺口补证

## 用户要求

用户把源服务器复核产物上传到R2（文件清单、JSON manifest、7份证据JSON、两份核查日志、复核结果与脚本），要求继续推进。源端结论：16个已记录SHA一致、14个补算SHA、50个依赖引用存在；WEHI目标端未核。

## 开始状态

- 仓库与origin/main一致（`0386572`），工作区干净。
- 本会话此前已在WEHI读过三原生zarr的obs/var与旧sample_metadata（见对话方案，未入库）。源端报告纠正了其中几处表述：联合层只读immune一份五重叠来源；旧sample_metadata不代表当前资格；v20261010 CM为1,290样本（1,303属于0825 CoVarNet冻结版本，需带run名引用）；不另建平行登记。

## 计划

1. 从R2取回上传产物，核对证据JSON的SHA。
2. 按`PATH_MAP.tsv`把清单32项与50个依赖映射到WEHI，逐项算SHA（目录与>2 GiB文件只核存在/大小，并对照迁移逐文件SHA）。
3. 处理清单"仍需补证"中能在WEHI独立完成的项：GSE299623/GSE169609范围依据、GSE278694理由、ADT精确分类、空间与临床扩展元数据的目标端状态。

## 执行与发现（取证 2026-10-31 11:00—11:23 UTC / 19:00—19:23 UTC+8）

### 1 上传产物

- R2 `server/XX/work/data_claim_audit_20261031/` 13个文件取回到`${PC_ROOT}/migration/local_upload/data_claim_audit_20261031/`。
- `evidence_hashes.json`列7份证据，6份SHA一致；`followup_probe.json`列在其中但未上传（已观察）。

### 2 WEHI目标端复核（已观察）

工具：[wehi_verify_handoff.py](../../../scripts/wehi_migration/wehi_verify_handoff.py)。回执`WEHI_target_verify_20261031.json`（SHA256 `212651a2…a15b`，含路径，留服务器）。

| 结果 | 数 | 说明 |
|---|---|---|
| SHA与源端一致 | 52 | 清单项及依赖小文件 |
| 目录存在 | 3 | core、public union、immune三zarr；迁移时已逐文件核对（83,409 / 38,296 / 24,431文件，0不一致，`migration/verify_relocated.tsv`） |
| 仅核大小 | 1 | ICB扩展的3.06 GB h5ad，源端未算SHA；大小一致，且在M2逐文件SHA清单中已核 |
| 不在WEHI | 3 | GSE162498 raw/filtered（Project_v1未迁移，无PATH_MAP项）、GSE278694 `sample_qc_filtered`（Project_v3旧run未迁移）；源端清单本就注明"无需回传矩阵" |

结论边界：清单所列证据与依赖在WEHI可用且与源端一致；不代表WEHI上所有运行都已与源端当前状态重新比对。

### 3 GSE299623与GSE169609不在选定范围（已观察）

`native_archive_collection_v1.json`的`outside_selected_source_scope`为GSE299623 94,028、GSE169609 17,478。继承自冻结52-GSE下游输入（`legacy_bank_descriptor_v1.json`："preserve GSE169609/GSE299623 exclusions"）：

- GSE299623：最终身份门槛未闭合（40群/54,099细胞，TODO ID-001）。
- GSE169609：`20261010_expression_annotation_downstream_v1`登记为`GSE169609_duplicate_source_hold`（17,454细胞）。WEHI记录中没有找到它与哪个来源重复——**未解决**。

### 4 GSE278694排除理由

WEHI上仅有预检（81,846×36,601，`pending_modality_and_provenance_review`，accepted 0）。仓库此前已登记为治理缺口（TODO GOV-002，2026-07-31日志）。WEHI证据无法闭合，需要原执行者的决定。

### 5 core中的非RNA特征：精确分类（已观察）

在WEHI重新下载GSE274229 GEO上44个人类样本（GPL24676）的`*_GEX_features.tsv.gz`。系列顶层附件是小鼠数据（ENSMUSG），不能使用。NCBI限速导致25个文件截断，重试后44/44通过`gzip -t`，SHA见`human_features.sha256`。工具：[core_feature_type_audit.py](../../../scripts/wehi_migration/core_feature_type_audit.py)，输出`feature_type_audit.json`（SHA256 `5d26bb9f…6beb`），位于`${PC_ROOT}/migration/public_refetch/GSE274229_features_20261031/`。

- 发布的`feature_types`：Antibody Capture 1,686行，去重229个（39个样本各22个，4个样本各207个，1个样本0个）；与RNA基因名无冲突。
- core中有229个特征与之完全对应（忽略大小写），`measured_in_n_gse`全部为1，只来自GSE274229。其中包括hashtag、同型对照和CD298/β2M样本标签。
- 先前的数字都不对：正则214（其中209个是真ADT，另外5个是GSE169609的`TOTALSEQC0951–0955`），宽正则222（GPT报告）。
- 影响：GSE274229这229列在合并前对象中计数全为0（156,166个细胞都没有抗体计数），对总计数和log1p_cp10k没有影响；它们只是空列，`measured_in_n_gse=1`容易被误读为"已测到"。
- **新发现**：GSE169609的5个`TOTALSEQC0951–0955`特征有真实计数，每个细胞都有，中位占总计数5.7%，最高76%。编号与TotalSeq-C细胞哈希标签一致——这是根据编号判断的**假设**，GSE169609的GEO `feature_types`尚未核。这些计数进了该来源的总计数和log1p_cp10k。该来源目前处于duplicate-source hold，当前下游不受影响；若重新纳入，必须先剔除。

### 6 目标端数据状态（已观察）

- 临床扩展：`metadata/`下36个GSE目录，各有`official_geo.json`（此前写的"约37个"应为36）；原始数据未迁移。
- 空间：Xenium入库31个counts h5ad（S18 17、S16 6、S07 2、C11 2、C05/S10/S11/S12各1）。IMC Steinbock归档已下载完成，本地`SHA256SUMS`复核OK，但该SHA是下载后自己生成的，没有对照Zenodo的md5。CosMx、MERFISH、Stereo-XCR、Visium只有原始文件，未入库。

## 修改的文件

- `scripts/wehi_migration/wehi_verify_handoff.py`（新）：按PATH_MAP核源端清单及依赖。
- `scripts/wehi_migration/core_feature_type_audit.py`（新）：按发布的feature_types分类core非RNA特征并计算计数占比。
- 本日志；`STATUS.md`新增顶部条目。

## 验证

- 两个脚本在WEHI真实数据上运行完成（exit 0），输出JSON SHA见上。
- `bash scripts/validate_repository.sh`、`git diff --check`结果见提交前记录。

## 失败尝试

- 第一次GEO下载用`awk`解析GSM列表时没有去掉`\r`，结果0个URL。修正后有25/44个文件被限速截断，删除后低速重试4轮全部恢复。

## 未验证 / 未解决

- GSE169609与哪个来源重复；其5个TotalSeq特征的GEO `feature_types`。
- GSE278694排除理由（GOV-002）。
- 跨研究DOI/GEO/CXG别名及患者重复；457个待核治疗阶段；癌种本体统一。
- `followup_probe.json`未上传。
- IMC对照Zenodo官方校验值。

## 建议下一步

1. 用户补`followup_probe.json`，并确认是否知道GSE169609是哪个来源的重复。
2. 写入特征使用规则：core的跨来源分析一律只用RNA特征；把GSE274229的229个空列和GSE169609的5个TotalSeq列列入排除表（只登记，不改冻结存储）。
3. GSE278694：由用户决定"永久排除并写ADR"或"按现有QC对象重新走资格门"。

## 建议commit message

Verify the source-server data-claim handoff on WEHI and classify core antibody features
