# 运行记录：20260725_120105_unified_pancancer_covarnet

## 运行身份

| 字段 | 值 |
| --- | --- |
| Scope | `Project_v3` |
| Run type | production / multi-stage |
| Start | 2026-07-25 12:01 +08:00 |
| Last verified successful write | 2026-07-30 11:18:02 +08:00 |
| Current status | `completed_with_audit_warnings` |
| Source commit | 未在运行记录中冻结，`not_recorded` |
| Current activity | 2026-07-31 进程检查未见相关活跃任务 |

本记录是对运行目录内脚本、配置、阶段记录、调度/直接运行日志和最终审计文件的事后核验。它不改写源运行目录，也不把大型对象、完整日志或细胞级数据复制进记录仓库。

## 目标

完成泛癌单细胞 discovery 队列的逐数据集 QC、human-only 共同基因构建、全量多方法整合、证据化细胞注释、CoVarNet 模块发现，并把冻结的细胞状态 mapper 和模块权重投影到三个免疫治疗队列。

## 输入与队列变化

| 检查点 | 数据集/分片 | 细胞 | 样本 | 基因 | 说明 |
| --- | ---: | ---: | ---: | ---: | --- |
| 初始 discovery 计划 | 34 datasets | 8,887,638 | 806 | 不同 | 含 6,635,520-cell 的 `GSE162498` |
| 实际逐数据集 QC | 33 datasets | 2,252,118 before QC | 797 | 不同 | `GSE162498` 单列 planning-only |
| QC pass | 33 datasets | 2,062,230 | 797 | 不同 | 移除 189,888 |
| human-only archival plan | 32 datasets / 39 shards | 2,040,048 | 794 | 11,063 common | 排除非人 `GSE211602` |
| final concat | 31 GSE / 38 shards | 1,963,745 | 773 | 11,063 | `GSE278694` 被移除，理由未记录 |

### Counts 语义

六个输入对象在矩阵审计后从 count-like `X` 复制 `.layers["counts"]`：`GSE131907`、`GSE179994`、`GSE188737`、`GSE222315`、`GSE269826`、`GSE278694`。其余数据集使用已有 counts layer。

## 阶段状态

| 阶段 | 主要处理 | 结果 |
| --- | --- | --- |
| inventory / matrix audit | 来源、物种、矩阵语义、样本键、分支审计 | 完成 |
| Phase 06 per-dataset QC | 样本内 QC、mt 阈值、Scrublet、h5ad/zarr 验证 | 33/33 完成 |
| Phase 07 archival build | human-only 基因交集、39 shards、capped branch | 完成 |
| final concat | 38 shards 拼接为全量 counts Zarr | 完成，31 GSE |
| full integration | HVG/PCA/Harmony、BBKNN、scVI、CellTypist、scIB | 完成 |
| evidence annotation | lineage graph clustering、marker aggregation、污染细化 | 完成，audit `PASS` |
| CoVarNet | sample abundance、rank survey、K9 nsNMF、网络和模块解释 | 完成，含已记录 fallback |
| immunotherapy validation | 三队列固定 mapper / fixed-W 投影 | 三个 audit 均 `PASS` |

`manifests/stage_status.tsv` 仍把若干后续阶段保留为 `pending`，因此该 manifest 不是当前完成状态的可靠来源；本记录以实际产物和最终审计为准。

## 逐数据集 QC：处理和参数

逻辑命令：

```bash
${SCANPY_PYTHON} ${PROJECT_ROOT}/Project_v3/scripts/run_phase06_sample_qc_worker.py \
  --run-root ${RUN_ROOT} \
  --gse-id <GSE> \
  --require-doublet-success \
  --update-run-record
```

| 参数 | 值 |
| --- | --- |
| floor min genes / total counts | 200 / 500 |
| strict min genes / total counts | 500 / 1,000 |
| MAD multiplier | 3 |
| mitochondrial percent | `min(25, max(15, median + 3*MAD))` |
| expected doublet rate | 0.04 |
| max removed doublet fraction/sample | 0.08 |
| doublet minimum sample cells | 500 |
| Scrublet min counts / min cells | 2 / 3 |
| Scrublet min gene variability percentile | 85 |
| Scrublet principal components | 30 |
| random seed | 0 |

没有 mt gene signal 时跳过 mt% 过滤并记录。实现细节：代码用样本 `n_cells >= 1000` 选择 strict profile，而不是用每细胞测序深度选择。

结果：

- 2,252,118 → 2,062,230 cells。
- 189,888 cells removed。
- 695 eligible / 695 Scrublet success / 0 failure。
- 每数据集写出 QC-filtered H5AD、Zarr、cell QC table 和阶段记录。

## 全量整合：处理、参数和资源

配置文件：`${RUN_ROOT}/config/full_hvg4000_multimethod_v1.json`  
配置标识：`20260727_full_hvg4000_multimethod_v1`

### 预处理和整合参数

| 参数 | 值 |
| --- | --- |
| input | `${RUN_ROOT}/objects/02_integrate/archival_shards_concat.zarr` |
| batch key | `gse_id` |
| sample key | `sampleID` |
| normalize / transform | 10,000 counts per cell / `log1p` |
| excluded technical genes | mt、ribo、Hb、ERCC、`MALAT1` |
| HVG | 4,000 from 6,500 candidates; `seurat_v3` |
| PCA | 50 |
| seed | 20260727 |
| Harmony | max iterations 20 |
| graph neighbors | 30 |
| BBKNN | neighbors within batch 2; trim 200; 50 PCs |
| scVI | latent 30; layers 2; hidden 128; epochs 30; batch 2,048 |
| CellTypist | `Immune_All_Low.pkl`; majority voting disabled |
| scIB subset | max 100,000; deterministic GSE × CellTypist stratification |
| scIB full graph metrics | enabled |

### 调度资源

提交入口：

```bash
bash ${RUN_ROOT}/scripts/submit_full_downstream_slurm.sh
```

| Job | CPU | RAM | GPU |
| --- | ---: | ---: | --- |
| metadata panel | 16 | 180G | none |
| HVG/PCA/Harmony | 64 | 850G | none |
| Harmony UMAP | 64 | 600G | none |
| BBKNN | 96 | 800G | none |
| scVI | 24 | 550G | 1 Tesla |
| CellTypist | 32 | 300G | none |
| scIB/finalize | 32 | 850G | 1 Tesla |

线程环境为 OMP/MKL/OpenBLAS/Numba 32。观察到的主 job IDs：panel 26070，失败的初次 HVG 26071，HVG retry 26077，Harmony UMAP 26078，BBKNN 26079，scVI 26080，CellTypist 26075。

### 输出和评价

- 统一 metadata panel：1,963,745 cells、773 samples、31 GSE、38 shards。
- Harmony / BBKNN / scVI graph nnz：91,983,248 / 230,498,834 / 87,695,572。
- CellTypist：97 labels，median confidence 0.4824。
- scIB 0.5.6 100k scaled total：Harmony 0.6156、BBKNN 0.2000、scVI 0.7959。
- 全图 iLISI：Harmony 0.01557、BBKNN 0.21751、scVI 0.00205。
- 全图 cLISI：Harmony 0.99164、BBKNN 0.98741、scVI 0.99272。
- 全图 graph connectivity：Harmony 0.42150、BBKNN 0.25459、scVI 0.46284。

结论边界：scVI 的子集综合分最高，BBKNN 的全图 iLISI 最高；方法之间存在指标权衡。下游注释使用 BBKNN 图。

## 证据化注释：处理和参数

输入：

```text
${RUN_ROOT}/objects/03_preprocess/20260727_full_hvg4000_multimethod_v1/integration_bbknn_full.h5ad
```

- 一级 Leiden resolution 0.1。
- lineage resolutions：Tcell 1.0、NK 0.8、ILC 0.6、Bcell 0.8、Plasma 0.6、Monocyte 0.8、Macrophage 0.9、Dendritic 0.8、Granulocyte 0.6、Mast 0.6、Epithelial 1.0、Endothelial 0.8、Fibroblast 0.9、Erythroid 0.5、Platelet 0.4、HSPC 0.7、Neural 0.6。
- seed 20260729。
- 小 cluster 合并阈值 30 cells；以 50-PC PCA centroid 合并到稳定 cluster。
- Neural 无稳定图 cluster 时使用 MiniBatchKMeans fallback。
- raw-count expression evidence：10,000 cells/chunk，10,000 library normalization + `log1p`，marker panel 194 genes，UMAP 可视化最多 180,000 cells。
- v2 污染细化：15 clusters、404,773 cells；细化后无剩余未验证 contamination clusters。

结果：

- 17 major lineages。
- 136 lineage clusters。
- 106 final annotations。
- confidence cluster counts：high 43、medium 22、low 56、unresolved 10、ambiguous 5。
- `${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/FINAL_DELIVERY_AUDIT.json` 为 `PASS`。

## CoVarNet：处理和参数

输入列映射：`sample_uid → sampleID`、`major_celltype → major`、`final_annotation → sub`。

| 参数 | 值 |
| --- | --- |
| min cells/sample | 100 |
| abundance normalization | min-max |
| correlation method | Pearson |
| final K | 9 |
| rank survey | 2–20; options `vp` |
| final NMF | nsNMF; `nrun=30`; options `v` |
| top states | 10; interpretation report top 20 |
| network correlation threshold | 0.2 |
| network FDR threshold | 0.05 |
| observed parallelism | 2 / 192 cores |

处理后 18 个不足 100 cells 的样本被移除，保留 1,962,479 cells / 755 samples。固定矩阵为 `W: 106 × 9`、`H: 9 × 755`；网络为 101 nodes / 374 edges。

恢复链：

1. R `sc_cm_recover` 因类型比较错误失败。
2. 以 `coef(nmf_final)` 读取 CM abundance，记录为 `coef_nmf_final_fallback`。
3. R 网络图因 `RColorBrewer::Set3` 最多 12 色、而 major classes 为 17 类而失败。
4. 使用 Python/networkx 从已完成的网络表重新渲染全局和分模块图；没有重新拟合 NMF。

## 三个免疫治疗队列：固定模块投影

共同参数：

| 参数 | 值 |
| --- | --- |
| mapper / module weights | frozen 106-state mapper / frozen 106 × 9 W |
| refit in validation | no |
| expression transform | count per 10,000 + `log1p` |
| prediction chunk | 5,000 cells |
| confidence threshold | 0.50 |
| sample abundance | soft frequency primary; hard label sensitivity |
| discovery scaling | min-max |
| projection | nonnegative least squares |
| paired statistics | patient-level Wilcoxon + BH where applicable |

### GSE195832

- 30,190 cells、8 samples、4 paired pre/post nivolumab patients。
- publication-style QC：counts > 1,500、counts < 15,000、genes ≥ 200、mt < 10%。
- 9/9 soft/hard median directions 一致；median rho 0.9940；没有样本超出 discovery min-max frequency 范围。
- CM06 在 4/4 配对上升，但 `q=0.641`；无 responder label。
- 当前细胞数不等于论文报告的 22,906，因此不是论文最终对象的精确复现。

### GSE123813

- 79,040 cells、23,121 genes、30 biological samples、15 complete pairs。
- 由于 sorting confounding，主分析只使用 7 个 sort-matched pairs，并对不同 sort fraction 等权重；15 个配对仅作敏感性分析。
- median soft-hard rho 0.9535；NNLS residual 0.324。
- 没有 responder 或 survival endpoint。

### GSE169246

- 489,490 cells、23,947 genes、78 samples、22 patients。
- 主分析：20 blood pre/post pairs；敏感性分析：7 blood pre/post/progression triples、11 same-site tumor pairs。
- mapper genes 3,999/4,000；median confidence 0.6906；high-confidence fraction 0.7493；median soft-hard rho 0.9362；NNLS residual 0.2992。
- 为避免重复大型对象，没有另写 mapped H5AD；cell prediction table 为权威细胞级结果。
- treatment arm/response 来自 legacy metadata，需以原论文/官方补充材料核验；`Prog` 仅表示 progression timepoint。

三个 `FINAL_AUDIT.json` 均通过。BH 校正后没有显著模块，因此当前只支持探索性方向。

## 环境

观察到的主要版本：

- Python 3.10.19，scanpy 1.11.5，anndata 0.11.4。
- scVI direct resume 使用 Python 3.11.14。
- R 4.3.3，NMF 0.28，Seurat 5.3.1，igraph 2.2.0。
- timezone：Asia/Shanghai。

环境只以逻辑名称 `scanpy`、`scvi`、`scrna_r` 记录；实际环境路径保留在外部脚本和日志中。

## 失败与恢复记录

这些失败均有后续成功结果，但应保留以解释最终产物来源：

1. 15-dataset QC validator 因 sparse-Zarr shape probe 实现失败；不是 QC 本身失败。修正 validator 后验证通过。
2. 33-dataset gene harmonization 因非人 `GSE211602` 失败；排除后 human-only 32-dataset 分支成功。
3. capped model HVG/Harmony 初版失败；v2 成功。
4. 全量 HVG/PCA/Harmony 初次 SLURM job 26071 失败；retry 26077 成功。
5. scVI job 26080 完成 30 epochs 训练后，代码把 dict 类型的 `model.history` 当作 DataFrame 调用 `.to_csv` 而失败；首个 direct resume 仍失败，r2 resume 复用已保存模型并成功。
6. CoVarNet abundance 和绘图分别使用已记录的 R/Python fallback 完成，见上一节。
7. 日志中的 locale warning 不影响计算结果。

## 审计警告

### 已观察事实

- 最终拼接不含 `GSE278694`。
- 过滤代码和 “shard intentionally deleted” 日志存在。
- 源 `stage_status.tsv` 未更新到最终完成状态。
- 运行时 source commit 未冻结。

### 尚不能确认

- `GSE278694` 被移除的科学理由。
- 恶性/CNV 分支是否属于本次运行且已完成。
- `GSE169246` legacy clinical metadata 与原始论文补充材料是否完全一致。

### 建议

- 用 ADR 或补充 run record 冻结 `GSE278694` 的范围决定。
- 对齐源 manifest，但保留历史失败和更正记录。
- 对低置信度/争议注释执行人工复核。
- 后续外部验证继续使用固定 mapper/W，并预先冻结临床 endpoint 和统计模型。

## 权威产物

- 最终 concat：`${RUN_ROOT}/objects/02_integrate/archival_shards_concat.zarr`
- multi-method final：`${RUN_ROOT}/objects/03_preprocess/20260727_full_hvg4000_multimethod_v1/integrated_multimethod_celltypist_final.h5ad`
- annotation audit：`${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/FINAL_DELIVERY_AUDIT.json`
- refined annotation object：`${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/objects/integration_bbknn_full_evidence_annotated_refined.h5ad`
- CoVarNet weights：`${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/covarnet_r/Covarnet_R_refined_v2_20260729_072200/tables/nmf_W_subcluster_by_cm.tsv`
- CoVarNet report：`${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/covarnet_r/Covarnet_R_refined_v2_20260729_072200/module_analysis/COVARNET_CM_ANALYSIS.md`
- validation audits：`${RUN_ROOT}/09_immunotherapy_validation/*/FINAL_AUDIT.json`

## 日志索引

完整日志保留在外部运行目录，本仓库只记录路径和关键失败摘要：

- QC：`${RUN_ROOT}/logs/01_qc/`、`${RUN_ROOT}/06_full_qc_integration/per_dataset/`
- concat 和 capped branch：`${RUN_ROOT}/07_unified_discovery_build/`
- full integration stage logs：`${RUN_ROOT}/logs/03_preprocess/20260727_full_hvg4000_multimethod_v1/`
- SLURM stdout/stderr：`${RUN_ROOT}/logs/scheduler/20260727_full_hvg4000_multimethod_v1/`
- scVI/scIB direct recovery：`${RUN_ROOT}/logs/direct/`
- annotation v1：`${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v1/logs/`
- annotation v2：`${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/logs/`
- CoVarNet R/log/report：`${RUN_ROOT}/08_annotation/20260729_bbknn_full_evidence_v2/covarnet_r/Covarnet_R_refined_v2_20260729_072200/`
- validation：`${RUN_ROOT}/09_immunotherapy_validation/20260730_gse195832_fixed_cm_v1/`
- validation：`${RUN_ROOT}/09_immunotherapy_validation/20260730_gse123813_fixed_cm_v1/`
- validation：`${RUN_ROOT}/09_immunotherapy_validation/20260730_gse169246_fixed_cm_v1/`

## 复现边界

- 参数、逻辑命令、环境版本和权威产物已记录。
- source commit、完整依赖 lockfile、`GSE278694` 排除理由尚未冻结，因此不能称为完全可复现。
- CellTypist 模型是免疫参考，混合肿瘤对象中的标签必须与 marker evidence 和人工复核共同解释。
- 所有临床投影均为验证性固定投影；不得用验证队列重新选择 K、重训 mapper 或重拟合 `W`。
