# 运行记录：immunotherapy_4cohort_20260803_v1

## 运行身份

| 字段 | 值 |
| --- | --- |
| Scope | `Project_v2` 四队列免疫治疗整合扩展 |
| Run type | production / multi-stage / resumable |
| Start | 2026-08-03 +08:00 |
| Snapshot | 2026-08-12 10:06 +08:00 |
| Status | `blocked_final_acceptance_policy_conflict` |
| Current stage | Phase 6 Stage 13 attempt 01 failed；19-slice cell2location 与严格 merge 已完成 |
| Source commit | `not_recorded`；Project_v2 路径未检测到 Git 元数据 |
| Code worktree status | `not_available` |
| Main environments | `scanpy`、`cell2loc_env` |
| cell2location stack | cell2location 0.1.5；scvi-tools 1.3.3；torch 2.6.0；CUDA unavailable |

逻辑根 `${IMMUNOTHERAPY_RUN_ROOT}` 指向：
`Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1`。

## 目标与模态边界

整合 GSE169246、GSE176021、GSE273952、GSE316195 共 1,518,479 个观测，同时保持 sc/snRNA 与 Visium 分模态建模。

- sc/snRNA：GSE169246、GSE176021、GSE316195，共 1,447,081 cells/nuclei。
- Visium：GSE273952，共 71,398 spots、19 张切片；spot 不解释为单细胞。
- GSE169246 为 CD45 富集，GSE176021 为淋巴细胞富集，不做可信肿瘤细胞判定。
- 可信肿瘤细胞只在 GSE316195 中定义为 Epithelial 候选、inferCNV 阳性与 CopyKAT aneuploid 的严格交集；scMalignant 只作辅助支持。
- GSE273952 的跨癌种 cell2location 参考只解释保守的免疫、基质和上皮组成，不声称 ccRCC 特异肿瘤亚型。

## 输入

| Dataset | Modality | Observations | Genes | Counts source | Role |
| --- | --- | ---: | ---: | --- | --- |
| GSE169246 | scRNA | 489,490 | 23,947 | `X`，经原始计数合同核验 | TNBC CD45+ 免疫图谱 |
| GSE176021 | scRNA | 874,608 | 26,340 | `X`，经原始计数合同核验 | NSCLC 淋巴细胞图谱 |
| GSE273952 | Visium | 71,398 | 24,146 | `layers['counts']` | ccRCC 空间图谱 |
| GSE316195 | snRNA | 82,983 | 33,438 | `layers['counts']` | PDAC 全微环境及肿瘤判定 |

输入保持只读；原 QC 对象不在本运行中覆盖。

## 分析流程、方法与参数

### 1. 元数据注册与临床字段

- 229/229 个样本或文库记录完成 GEO 匹配。
- GSE169246 的论文 Supplementary Table 1 覆盖 22 位患者、78/78 个文库的疗效字段。
- GSE273952 从 Figure S8 进行颜色注册与 NNLS 参数网格核验：27 套预定义颜色核心、空间网格和平滑参数；患者接受门槛为可评估配置不少于 9 且参数共识率不低于 0.80。
- GSE273952 最终只高置信写入 12/19 位患者的 response：R=9、NR=3；其余 7 位保留 unresolved，不强行补全。

### 2. 双细胞审计

- Scrublet 按 sc/snRNA 文库独立运行，使用原始 counts。
- 1,447,081 个 cells/nuclei 中识别 6,678 个高置信双细胞，比例 0.46%。
- 双细胞从 PCA、CNV 和注释建模中排除，但保留在原 QC 对象和审计表中。
- 后续建模输入为 1,440,403 个 singlets、22,388 个共享基因。

### 3. sc/snRNA 联合整合

处理顺序：

```text
raw counts
→ library-size normalize 10,000
→ log1p
→ 4,000 HVG
→ PCA 50
→ Harmony(dataset)
→ kNN 30
→ UMAP
→ Leiden 0.3 / 0.6 / 1.0
```

- Harmony 只修正 PCA 表征，不改写表达矩阵；观察到 3 次迭代收敛。
- 原始 counts 继续保留在 `layers['counts']`，用于需要原始计数的下游模型。

### 4. 分层注释

- 大类注释联合 marker module score、Leiden 0.6 聚类共识和 cluster DE marker。
- 随后按大类独立重做 HVG、PCA、邻居、UMAP、Leiden 与 marker 审核。
- 标签格式为 `major_subtype_marker`；ambiguous、uncertain、unresolved 语义保留，不以黑名单 marker 强制命名。
- T、NK、B、Plasma、Epithelial 等谱系在数值稳定性需要时使用 `robust_dispersion_qcut_log1p`；Myeloid 使用 raw-count batch-aware Seurat v3 HVG。
- NK 最多使用 300,000 个细胞建模；其余已执行谱系使用全量可用细胞。
- 最终覆盖 1,440,403 个细胞，形成 43 个最终标签；Fibroblast、Endothelial、Mast 因可用细胞过少保留跳过证据。

### 5. GSE316195 恶性判定

- 82,549 个 GSE316195 nuclei 一对一映射回原始 QC counts；22 个文库。
- Epithelial 候选 54,384；同文库参考细胞 2,841。
- 每文库顺序运行 inferCNV → CopyKAT → scMalignant。
- inferCNV 与 CopyKAT 均使用原始 counts；scMalignant 使用 library-size 10,000 normalization，不做 log1p。
- 22/22 文库方法运行无失败；两种 CNV 方法均可评估 22,319 个候选，严格交集得到 5,648 个 `Epithelial_Malignant_CNV`。
- scMalignant 对严格可信细胞的辅助支持率为 0.775，但不进入严格交集判定。

### 6. GSE273952 Visium 独立分析

- 保留 71,398 spots、24,146 genes、19 libraries；恢复标准 AnnData spatial contract、组织图、scalefactors 和坐标。
- raw counts 上选择 4,000 HVG；10,000 normalization + log1p 后进行 PCA 50、邻居、UMAP 和 Leiden 0.3/0.6/1.0。
- 不按 library 做 Harmony，因为 library 与患者一一对应，强制批次校正会同时去除患者生物学差异。
- 当前 spot 聚类只表示 dominant marker program；定量细胞组成由 cell2location 给出。

### 7. cell2location 参考与空间模型

参考准备：

- 从最终 sc/snRNA 图谱按 dataset 平衡抽取 28,357 个 raw-count cells，覆盖 40 个标签。
- 每标签最少 100、最多 750 个参考细胞；实际最小 219、最大 750。
- 排除 uncertain labels，并把 5,648 个严格恶性细胞更新为恶性标签。
- RegressionModel 输入 28,357 cells、16,544 signature genes、40 factors；`labels_key=cell2location_label`、`batch_key=dataset`。
- RegressionModel：250 epochs、`train_size=1`、learning rate 0.002；posterior samples 1,000、posterior batch size 2,500。

每张 Visium 切片：

- 与 signature 取基因交集；当前切片使用 15,874 个共享基因，少于 200 时硬失败。
- `N_cells_per_location=30`。
- `detection_alpha=200`。
- spatial training 1,000 epochs、`train_size=1`、`batch_size=None`。
- posterior samples 1,000；posterior batch size 为该切片全部 spots。
- torch threads 12；本次为 CPU 运行。
- 导出 `q05_cell_abundance_w_sf`；要求所有值有限、非负，且每个 spot 的 abundance sum 大于 0。
- top label 仅作显示；定量分析保留完整 40-factor abundance matrix 及按 major 汇总矩阵。

### 8. 严格合并与最终审计

- 19 张切片逐切片使用独立 attempt；完整成功切片可恢复，失败 attempt 原样保留。
- 合并必须对 19 张切片和全部 spot ID 做精确覆盖，禁止缺失、额外或重复 spot。
- 计划输出全因子 q05 abundance、major abundance、轻量 H5AD、逐 spot 注释、样本汇总和 merge audit。
- Phase 6 只在严格恶性结果与 19 张切片 merge 同时完成后运行，生成 `FINAL_ANALYSIS_SUMMARY.md`、`FINAL_DELIVERABLE_AUDIT.json`、`FINAL_DELIVERABLE_INDEX.tsv`。

## 2026-08-06 05:29 历史进度快照

### 已观察事实

- 元数据、Scrublet、sc/snRNA 整合、分层注释、GSE316195 恶性合并、GSE273952 空间整合、cell2location 参考训练均已完成。
- RegressionModel 审计记录 250 epochs、1,000 posterior samples，参考 signature 合同为有限、非负。
- 19 张切片中 6 张已有 `COMPLETED.json`，即 31.6%。
- 第 6 张切片 GSM8440299 已完成并写出；第 7 张 GSM8440300 已开始 1,000-epoch spatial training，主 Python 进程持续高 CPU 运行。
- 第 1–5 张切片的 observed elapsed time 约 10–92 分钟/张，随 spot 数量变化明显，不能据此给出可靠固定 ETA。
- Phase 6 watcher 已启动，但只是在等待完整 cell2location merge；最终交付审计尚未生成。

## 已知失败与恢复

- 初次 GEO 预检遇到远端断连，retry 成功；失败 attempt 保留。
- Phase 2 首次读取 Phase 1 H5AD 时遇到文件锁；增加 `COMPLETED.txt`、h5py 只读加锁、关键 key/shape 和首末 cell ID 就绪检查后，attempt 02 成功。
- cell2location smoke 的两个失败分别来自训练 history 维度兼容和 partial signature/label-map 精确匹配；v2/v3 兼容层在验证 signature 有限、非负且 label map 覆盖后恢复，smoke attempt 03 成功。
- 这些恢复不改变 raw-count、标签或空间丰度的科学合同。

## 输出与证据

| Artifact | Logical path | Snapshot status |
| --- | --- | --- |
| 四队列审计 | `${IMMUNOTHERAPY_RUN_ROOT}/IMMUNOTHERAPY_4COHORT_AUDIT.md` | continuously updated |
| sc/snRNA final annotation | `${IMMUNOTHERAPY_RUN_ROOT}/outputs/annotation/scsn_integrated_final_annotated.h5ad` | completed |
| malignancy audit | `${IMMUNOTHERAPY_RUN_ROOT}/outputs/malignancy/GSE316195_malignancy_audit.json` | completed |
| spatial contract audit | `${IMMUNOTHERAPY_RUN_ROOT}/outputs/spatial/GSE273952_spatial_contract_audit.json` | completed |
| reference audit | `${IMMUNOTHERAPY_RUN_ROOT}/outputs/spatial/cell2location/reference/scsn_cell2location_reference_audit.json` | completed |
| regression audit | `${IMMUNOTHERAPY_RUN_ROOT}/outputs/spatial/cell2location/reference/regression_model_audit.json` | completed |
| per-slide results | `${IMMUNOTHERAPY_RUN_ROOT}/outputs/spatial/cell2location/samples/` | 6/19 completed; slide 7 active |
| final merge | `${IMMUNOTHERAPY_RUN_ROOT}/outputs/spatial/cell2location/final/` | pending |
| final acceptance | `${IMMUNOTHERAPY_RUN_ROOT}/reports/FINAL_DELIVERABLE_AUDIT.json` | pending |

## 解释边界与下一步

- cell2location 输出是多细胞状态的 spot-level q05 abundance，不是单一 cell identity。
- 跨 TNBC/NSCLC/PDAC 参考投射到 ccRCC 只支持保守组成解释；不得据此声称 ccRCC 特异恶性亚型。
- GSE273952 仍有 7/19 位患者 response 未稳定解析，不得补写为 R 或 NR。
- 2026-08-06 当时计划由后台任务完成剩余 13 张切片、严格 merge 和 Phase 6 最终审计；该计划的执行结果见下方 2026-08-08 快照。

## 2026-08-08 06:16 审计快照

- Phase 5 已结束：19/19 张切片均有 `COMPLETED.json`，cohort manifest 为 19 completed、0 failed。
- 严格 merge 已完成 71,398 spots 和 40 factors；spot 覆盖精确，missing=0、extra=0；q05 abundance finite fraction=1、negative values=0、zero-sum spots=0。该空间分支满足定量丰度合同，不能再记为进行中。
- GSE316195 malignancy audit 仍为 22 文库、82,549 cells、54,384 candidates、5,648 trusted malignant，并记录 `samples_with_method_failure=0`。
- Stage 13 attempt 01 于 2026-08-06 21:09 失败，未生成 `reports/FINAL_DELIVERABLE_AUDIT.json`。失败门槛把 11 个 `infercnv_status=skipped_insufficient_cells` 视作 required-method failures；这些文库的 CopyKAT 均为 success，属于样本内参考不足与验收策略不一致，不是 inferCNV/CopyKAT 进程崩溃。
- 快照时未发现免疫治疗/cell2location 活跃进程。因此任务不能标为 completed，当前为最终验收阻断。
- 下一步：先书面确定“参考不足 skip”是否为允许的 terminal 状态；若允许，修正最终审计的状态映射并用新 attempt 重跑 Stage 13；若不允许，应单独补充参考策略与敏感性分析。两种选择都不得覆盖已通过的 cell2location 结果。

## 2026-08-12 10:06 审计快照

- 2026-08-08 之后没有发现新的 Phase 6/Stage 13 attempt，快照时也没有免疫治疗或 cell2location 活跃进程。
- 已通过的空间结果没有变化：19/19 张切片、71,398 spots、40 factors；严格 merge 仍满足精确 spot 覆盖、所有丰度有限且非负、无 zero-sum spot。
- Stage 13 仍停留在 attempt 01。11 个 `infercnv_status=skipped_insufficient_cells` 与最终验收 required-method 语义的冲突尚未解决，`reports/FINAL_DELIVERABLE_AUDIT.json` 仍未生成。
- 因此状态继续为 `blocked_final_acceptance_policy_conflict`；不能把空间分支成功外推为四队列最终交付完成。
- 下一步仍需先形成书面验收决定，再以新 attempt 重跑 Stage 13；不得覆盖或重算已经严格通过的 19-slice cell2location merge。
