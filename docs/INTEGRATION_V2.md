# 整合数据集 v2：`integration_v2.zarr`（去富集、统一环境RNA校正、标准化、去批次、降维、UMAP）

构建：2026-11-01/02（UTC）。代码：`scripts/source_ingest/`（读取、QC）、`scripts/integration_v2/`（decontX、整合、作图、报告）。组成表：[inventories/integration_v2/](../inventories/integration_v2/)。上一版的数据层见[ANALYSIS_READY_ZARR_V1](ANALYSIS_READY_ZARR_V1.md)。

## 一句话

72个来源、**4,720,931个细胞**、1,258个单位、1,014名患者、28个癌种。所有细胞先统一做decontX环境RNA校正，再按此前TME标准配方（seurat_v3 HVG → PCA → 按来源Harmony）整合，得到近邻图和UMAP。没有做细胞类型注释或聚类定名。

## 位置与结构

| 项目 | 内容 |
|---|---|
| 路径 | `${PC_ROOT}/integration_v2/run_v2/integration_v2.zarr`（WEHI scratch，107 GB，319,665个文件；没有备份，30天未访问会被清理） |
| `X` | log1p(normalize_total 1e4)，基于decontX校正后的计数，float32 |
| `layers['counts']` | 原始整数计数 |
| `layers['decontx']` | decontX校正后的计数 |
| `obsm` | `X_pca`（50）、`X_pca_harmony`（50）、`X_umap` |
| `var['selected_HVG']` | 用于PCA的4,000个HVG。HVG表（每个来源的排名、出现的批次数）在`run_v2/report/hvg_v2.tsv` |
| `obs` | 沿用v1的32个统一字段（来源、单位、患者、癌种/OncoTree、取样部位、治疗、新鲜/冻存、技术、捕获范围、回溯坐标、QC标记），另加`decontX_contamination`、`decontX_status`、`total_counts_decontx` |
| `uns` | 配方、过滤规则、各来源的过滤计数、代码commit |
| 图 | `run_v2/figures/`：按癌种/来源/批次/取样部位/治疗着色的UMAP，以及16个marker基因的叠加图（图件只放在服务器和对话中，不进公开仓库） |

## 数据来源与筛选

1. **现有分析集v2**（40个来源，当前TME分析使用的细胞，排除3个不合格来源）。
2. **补充v1**（13个通过详细核查的GEO来源）中，按用户规则再排除2个：GSE160269（CD45±分选后合并）、GSE163203（Microwell-seq，中位UMI为跨来源中位数的0.16倍），余11个。
3. **补充v2**（覆盖TCGA缺失或不足的癌种：头颈鳞癌、皮肤黑色素瘤、卵巢、子宫内膜、软组织肉瘤、DLBCL、间皮瘤）：14个来源。
4. **补充v3**（常见癌种的更多独立来源：乳腺、结直肠、胃、肾、胶质瘤、非小细胞肺癌）：7个来源。

新来源的统一规则（用户确定）：
- 不使用富集数据：谱系分选（FACS/MACS/CD45±等）排除；只去死细胞的活细胞FACS保留。
- 来源中位UMI低于跨来源中位数（4,127）的1/3或高于其3倍的排除。
- 测到的基因少于15,000的排除（作者已做过基因过滤）。
- 少于3名患者的排除。
- 只用肿瘤样本（原发、转移、癌栓）；癌旁、正常淋巴结、血液、癌前病变、类器官、空间转录组样本都去掉。
- 与现有来源属同一研究（同一PMID）的排除。

按这些规则排除的候选及原因见[v2判定](../inventories/source_search/v2/candidates_detailed_check_v2.tsv)、[v3判定](../inventories/source_search/v3/candidates_detailed_check_v3.tsv)和[任务日志](../logs/codex/2026/2026-10-31_data-gap-goals.md)。

## 处理流程

| 步骤 | 内容 |
|---|---|
| QC | 现有细胞沿用原QC与双细胞判定。新来源：检出基因≥200且UMI≥500；按样本的线粒体上限min(25%, max(15%, 中位数+3MAD))；按样本Scrublet；样本QC后不足100个细胞即剔除 |
| 环境RNA | **全部72个来源统一用decontX**，按单位运行（快速Leiden聚类作为z；celda 1.18.2的`decontX()`，参数与项目既有worker相同）；没有单位失败。去掉污染率≥0.5的细胞。保留细胞的污染率中位3.7%，P90为25.6% |
| 共同基因 | 72个来源都测到的12,563个基因 |
| HVG | seurat_v3按来源分别计算（原始计数），再按scanpy batch_key的方式合并排名：取6,500个候选，去掉技术基因（MALAT1、MT-、RPS/RPL、HB*、ERCC-）后保留4,000个 |
| 标准化 | decontX计数 → normalize_total 1e4 → log1p |
| 降维 | 中心化ARPACK PCA，50个主成分，不做缩放 |
| 去批次 | Harmony，batch = 来源，theta 2，lambda 1；5轮收敛 |
| 图与嵌入 | 基于Harmony的15近邻；UMAP（scanpy默认参数） |
| 资源 | bigmem节点，7 h，峰值内存约66 GB |

## 验证

回读检查：细胞ID全部唯一；counts层为整数；X为4,720,931 × 99,071；三个嵌入的形状正确，UMAP坐标全部为有限值；选中的HVG为4,000个（[结果](../inventories/integration_v2/readback_v2.json)）。跨来源细胞计数指纹查重（72个来源，全部细胞，用原始计数）：完整指纹跨来源相同为0；top50指纹有519个细胞、76个来源对相同，分散在不相关的癌种之间（最多为GSE305078与GSE310797的61个），按低复杂度巧合处理（[结果](../inventories/integration_v2/fingerprint_summary.json)）。

## UMAP的初步阅读（marker叠加图，不是注释）

T/NK（CD3E、NKG7）、B（MS4A1）、浆细胞（JCHAIN）、髓系（CD68、LYZ）、肥大细胞（TPSAB1）、成纤维细胞/肌成纤维细胞（COL1A1、ACTA2）、内皮（PECAM1）各自形成跨来源、跨癌种混合的大团；上皮/恶性细胞（EPCAM、KRT19）位于中央并有许多按来源分开的小岛，这符合恶性细胞的肿瘤特异性；另外可见黑色素细胞（PMEL）、胶质（GFAP）和增殖（MKI67）区域。这只是描述，不是批次校正质量的定量评估。

## 组成

### composition_by_batch

| data_batch | cells | units | patients | sources |
|---|---|---|---|---|
| existing_analysis_set_v2 | 2,676,846 | 748 | 611 | 40 |
| supplement_v1 | 678,218 | 146 | 111 | 11 |
| supplement_v2 | 898,033 | 231 | 194 | 14 |
| supplement_v3 | 467,834 | 133 | 98 | 7 |

### composition_by_cancer

| analysis_main | cells | units | patients | sources |
|---|---|---|---|---|
| NSCLC | 598,822 | 190 | 157 | 9 |
| Breast cancer | 530,092 | 99 | 90 | 3 |
| Head and neck squamous cell carcinoma | 406,249 | 119 | 114 | 5 |
| Renal cell carcinoma | 357,665 | 67 | 51 | 4 |
| Cutaneous melanoma | 317,622 | 54 | 38 | 3 |
| Colorectal cancer | 291,910 | 156 | 142 | 7 |
| Pancreatic ductal adenocarcinoma | 271,427 | 65 | 64 | 3 |
| Diffuse glioma | 239,624 | 76 | 28 | 2 |
| Bladder cancer | 193,308 | 18 | 14 | 2 |
| Prostate cancer | 162,282 | 69 | 38 | 2 |
| Gastric cancer | 159,904 | 51 | 46 | 3 |
| Osteosarcoma | 154,610 | 23 | 23 | 3 |
| Ovarian cancer | 125,274 | 40 | 29 | 4 |
| Papillary thyroid carcinoma | 118,829 | 23 | 14 | 2 |
| Hepatocellular carcinoma | 112,203 | 42 | 26 | 4 |
| Esophageal cancer | 107,898 | 18 | 15 | 2 |
| Cervical cancer | 91,437 | 26 | 18 | 2 |
| Uveal melanoma | 81,980 | 17 | 17 | 2 |
| Nasopharyngeal carcinoma | 80,516 | 10 | 10 | 1 |
| Pancreatic neuroendocrine tumor | 77,517 | 17 | 17 | 1 |
| Acral melanoma | 72,011 | 16 | 15 | 2 |
| Intrahepatic cholangiocarcinoma | 49,077 | 24 | 16 | 3 |
| Endometrial cancer | 30,990 | 5 | 5 | 1 |
| Basal cell carcinoma | 27,528 | 4 | 4 | 1 |
| Renal medullary carcinoma | 22,030 | 9 | 9 | 1 |
| Soft tissue sarcoma | 14,537 | 7 | 6 | 1 |
| Diffuse large B-cell lymphoma | 13,161 | 8 | 4 | 1 |
| Mesothelioma | 12,428 | 5 | 5 | 1 |

### composition_by_sampling_site

| sampling_site | cells | units | patients | sources |
|---|---|---|---|---|
| metastatic_site | 691,945 | 148 | 131 | 21 |
| primary_organ_site | 3,550,842 | 981 | 800 | 61 |
| tumour_thrombus | 5,248 | 2 | 2 | 1 |
| unspecified | 472,896 | 127 | 104 | 9 |

### composition_by_treatment_status

| treatment_status | cells | units | patients | sources |
|---|---|---|---|---|
| on_or_after_treatment | 555,090 | 122 | 92 | 14 |
| other | 84,865 | 16 | 12 | 2 |
| unknown | 1,903,682 | 575 | 462 | 32 |
| untreated_or_pretreatment | 2,177,294 | 545 | 480 | 39 |

### composition_by_fresh_frozen

| fresh_frozen | cells | units | patients | sources |
|---|---|---|---|---|
| cryopreserved_or_frozen | 181,705 | 69 | 50 | 3 |
| fresh | 2,433,683 | 615 | 496 | 36 |
| fresh_or_cryopreserved | 28,902 | 22 | 22 | 1 |
| unknown | 2,076,641 | 552 | 446 | 33 |

### composition_by_technology

| technology | cells | units | patients | sources |
|---|---|---|---|---|
| 10x | 1,823,667 | 486 | 383 | 30 |
| BD Rhapsody WTA | 89,369 | 15 | 15 | 1 |
| Singleron GEXSCOPE | 131,049 | 9 | 5 | 1 |
| scRNA (see source registry) | 2,676,846 | 748 | 611 | 40 |

### composition_by_qc_flag

| qc_flag | cells | units | patients | sources |
|---|---|---|---|---|
| floor_cell_calling_on_unfiltered_matrix | 131,049 | 9 | 5 | 1 |
| inherited_fit_cell | 2,676,846 | 748 | 611 | 40 |
| qc_v1_pass | 1,913,036 | 501 | 398 | 31 |

各来源的明细：[composition_by_source.tsv](../inventories/integration_v2/composition_by_source.tsv)。

## 覆盖与差距

- 去掉污染率≥0.5的细胞之前合计约5.07M，之后为4.72M。
- TCGA的33个癌种中，大约22个有对应：BLCA、BRCA、CESC、CHOL（ICC）、COAD/READ、DLBC、ESCA、GBM/LGG（弥漫性胶质瘤）、HNSC、KIRC、LIHC、LUAD/LUSC、MESO、OV、PAAD、PRAD、SARC（软组织肉瘤和骨肉瘤）、SKCM、STAD、THCA、UCEC、UVM。没有的：ACC、KICH/KIRP（只有零星样本）、LAML（骨髓，不属于组织TME）、PCPG、TGCT、THYM、UCS（只有1例卵巢癌肉瘤）。
- 仍只有1个来源的癌种：鼻咽癌、PanNET、子宫内膜癌、BCC、肾髓质癌、软组织肉瘤、DLBCL、间皮瘤。

## 使用注意

1. 跨来源比较使用`X`或HVG/嵌入；单基因比较只用共同基因，未测到的基因不能当0。
2. 现有部分与新来源的上游QC不同（`qc_flag`）；GSE277524只用QC下限做细胞识别。敏感性分析可按`data_batch`/`qc_flag`分层。
3. 患者键只在来源内唯一；同一患者在不同研究中出现无法从计数矩阵排除。
4. 治疗状态和新鲜/冻存仍有大量未知（分别为1.90M和2.08M个细胞）。
5. 现有部分旧的Leiden/scVI/注释等结果可经`origin_store` + `origin_native_row`回到原存储读取；本版UMAP是新的嵌入，不替代旧服务器冻结的结果。
