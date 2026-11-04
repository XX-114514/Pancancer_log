# 分析就绪数据集 v1：`pancancer_analysis_ready_v1.zarr`

构建时间2026-11-01（UTC）。代码：`scripts/source_ingest/`；本页的组成表：[inventories/analysis_ready_v1/](../inventories/analysis_ready_v1/)。

## 一句话

把**现有分析集v2**（40个来源，排除3个不合格来源后剩余的、当前TME分析实际使用的细胞）和**补充来源v1**（13个新GEO来源，经过统一QC）合并成一个zarr：**3,818,042个细胞、53个来源、1,020个单位、793名患者、22个癌种**。X为原始整数计数，基因以symbol为键。这一步**没有**做标准化、批次校正、降维或细胞类型注释。

## 位置与规模

| 项目 | 值 |
|---|---|
| 路径 | `${PC_ROOT}/analysis_ready_v1/pancancer_analysis_ready_v1.zarr`（labshare /vast/projects/Spatial/xin，2026-11-04起） |
| 大小 | 19 GB，4,320个文件（zarr v2，元数据已consolidate；2026-11-04 X重写为4,194,304元素的chunk以减少文件数，数值与原store逐值一致；原66,518个文件的版本仍在scratch） |
| 细胞 × 基因 | 3,818,042 × 94,200 |
| 所有53个来源都测到的基因 | 12,498 |
| 读取 | `anndata.read_zarr(path)`，或按需用`zarr`/`anndata.io.read_elem`只读obs/var |

## 数据结构

| 部分 | 内容 |
|---|---|
| `X` | CSR，int32，**原始整数计数**。现有来源取自原存储里当前分析使用的计数键（`obs.origin_counts_key`：core为`X`，immune/public为`layers/counts`），**没有**扣除环境RNA；新来源取自GEO补充文件。重复的基因symbol已求和；非RNA特征（抗体/hashtag，见[排除表](../inventories/feature_exclusions/non_RNA_features_v1.tsv)）已去除 |
| `var` | 索引为基因symbol；`n_sources_measured`、`measured_in_all_sources`、`measured_in_all_existing_sources` |
| `varm['measured_in_source']` | 94,200 × 53布尔矩阵，列顺序见`uns['source_order']`。某来源没测到的基因在X中是**结构性0**，不能当作表达为0 |
| `uns` | `schema`、`code_commit`、`X`/`unmeasured_genes`/`existing_cells`/`supplement_cells`/`not_done`等文字说明（JSON字符串） |
| `obs` | 32列，见下表 |

### obs字段

| 字段 | 含义 |
|---|---|
| `cell_id` | 全局唯一ID：`来源|原细胞名`（新来源为`来源|样本|barcode`） |
| `source`、`data_batch` | GEO/CXG编号；`existing_analysis_set_v2`或`supplement_v1` |
| `unit_id`、`sample_id` | 分析单位：现有部分为资格审查单位（qualification unit），新部分为`来源|样本` |
| `patient_id`、`patient_id_rule` | 带来源前缀的患者键。规则：`recorded`（原有记录）、`source_rule`（命名规则归并，漏洞2）、`assumed_one_unit_per_patient`、`GEO_sample_record`。**只在来源内唯一，不代表跨研究的独立患者认证** |
| `analysis_main`、`analysis_fine`、`oncotree_code`、`cancer_label_resolution` | 癌种统一v1（OncoTree），见[方法](../methods/cancer-type-harmonization-v1-20261031.md) |
| `sampling_site`、`material_detail` | `primary_organ_site` / `metastatic_site` / `tumour_thrombus` / `unspecified`；原始材料描述 |
| `treatment_status`、`treatment_detail`、`treatment_rule` | `untreated_or_pretreatment` / `on_or_after_treatment` / `other` / `unknown`；证据来源（漏洞1） |
| `fresh_frozen`、`fresh_frozen_detail` | `fresh` / `cryopreserved_or_frozen` / `fresh_or_cryopreserved` / `unknown`；原始描述 |
| `technology`、`suspension`、`capture_scope` | 平台；cell/nucleus；捕获范围（见下文"注意事项"） |
| `origin_store`、`origin_counts_key`、`origin_native_row`、`origin_cell_name` | 回溯到原存储或GEO文件的坐标（现有细胞可按行号回到原zarr，取到环境RNA校正层、原注释等） |
| `ambient_correction`、`doublet_method` | 环境RNA校正状态；双细胞处理方式 |
| `author_cell_type` | 新来源中作者提供的细胞类型（GSE149614、GSE160269、GSE242230），只作参考 |
| `qc_flag` | `inherited_fit_cell`、`qc_v1_pass`、`low_depth_platform`、`floor_cell_calling_on_unfiltered_matrix`、`sorted_fractions_recombined` |
| `total_counts`、`n_genes`、`pct_mt` | 由X统一重新计算 |

## 怎样来的

### 1 现有分析集v2（40个来源，2,871,696个细胞）
- 细胞集合 = 当前TME标准整合（`20261029_TME_standard_HVG_PCA_Harmony_v1`）使用的2,946,498个细胞（`private_cell_axes_v1.npz`中的来源编号与原生行号），去掉用户确认排除的3个来源（GSE163558、CXG SPECTRUM卵巢、GSE228501；[排除记录](../inventories/data_gaps/source_exclusions_v1.tsv)）。逐来源的细胞数与运行记录完全一致（合计2,871,696）。
- 细胞→单位：原生表达单位视图（`20261026_archive_expression_unit_views_v1`）；GSE166555按作者样本桥拆成13个单位。
- 单位级元数据：材料表（`20261030_TME_CM_technical_material_design_v1`），加上本轮漏洞1/2/5的补证结果（治疗、患者键、取样部位）和新鲜/冻存证据。
- 基因：只保留各来源实际测到的基因（core用`harmonized_per_gse`的基因列表，immune/public用`varm['measured_in_source']`），统一为symbol。
- 751个单位中有3个在拟合集中没有细胞（GSE176031、GSE178341、GSE274229各1个，均不足100个细胞），所以这里是748个单位。

### 2 补充来源v1（13个来源，946,346个细胞）
- 来源：A级候选详细核查中判定"通过"的13个（[核查表](../inventories/source_search/candidates_detailed_check_v1.tsv)），按冻结的[检索规格](../methods/tissue-source-search-spec-v1-20261101.md)筛选。从GEO下载（约5.2 GB），每个来源一个只读目录并生成`SHA256SUMS`。
- 样本选择：只保留肿瘤样本（原发、转移、癌栓）；癌旁正常、正常淋巴结、血液排除，与现有分析集口径一致。GSE215120只保留肢端黑色素瘤（皮肤型排除）。
- QC v1（项目native分支既有规则）：检出基因≥200且UMI≥500；按样本线粒体比例≤min(25%, max(15%, 中位数+3MAD))；按样本做Scrublet，去掉预测的双细胞；QC后不足100个细胞的样本剔除。没有估计环境RNA。
- 逐来源处理与异常：

| 来源 | 癌种 | 细胞（QC后） | 样本/患者 | 说明 |
|---|---|---|---|---|
| GSE163203 | PTC | 7,987 | 8/8 | **Microwell-seq，测序深度很低**（中位61–491 UMI），qc_flag=`low_depth_platform` |
| GSE191288 | PTC | 25,396 | 6/3 | 每名患者左右叶各1个样本 |
| GSE162454 | 骨肉瘤 | 43,553 | 6/6 | 未治疗，新鲜 |
| GSE299023 | 骨肉瘤 | 23,272 | 6/6 | 未治疗活检；CHCO样本为活细胞冻存；活细胞FACS |
| GSE205013 | PDAC | 239,194 | 27/27 | 20例未治疗、7例化疗后；原发17、肝转移10 |
| GSE242230 | PDAC | 30,628 | 22/22 | GEO矩阵未过滤，按作者注释的细胞列表取细胞 |
| GSE149614 | HCC | 41,933 | 13/10 | 原发、门静脉癌栓、淋巴结转移 |
| GSE125449 | HCC + ICC | 9,230 | 18/18 | 作者已过滤；P编号在Set1/Set2间重复，按样本计患者（论文为19名患者） |
| GSE189903 | HCC + ICC | 77,378 | 25/7 | 多区域；患者 = 病例号 + C/H |
| GSE277524 | 膀胱 | 150,840 | 9/5 | Singleron；GEO只有未过滤矩阵，**按QC下限做细胞识别**；P03_R2在GEO上的文件解压后全是空字节，无法读取 |
| GSE160269 | ESCC | 188,398 | 118/60 | **CD45+与CD45−分别分选建库**，按患者合并 |
| GSE203067 | ESCC | 51,379 | 7/4 | 原发与淋巴结转移；术前未治疗 |
| GSE215120 | 肢端黑色素瘤 | 57,158 | 7/6 | AM3有抗PD-1治疗前/后两个样本 |

### 3 合并与验证
- 先把每个来源的obs统一成32个字段，再在磁盘上按基因outer join拼接（`anndata.experimental.concat_on_disk`），然后写入基因×来源的测量矩阵和说明。
- 独立回读（`verify_analysis_zarr_v1.py`）：细胞ID全部唯一；X全为整数；53个来源的细胞数与总计数都和输入完全一致（[结果](../inventories/analysis_ready_v1/verify_summary.json)）。
- 跨来源查重（细胞计数指纹，53个来源全部细胞）：完整指纹跨来源相同的细胞为0，NIH同一团队的GSE125449与GSE189903之间也没有。top50指纹有210个细胞、30个来源对相同，分散在不相关的癌种之间，按低复杂度巧合处理（[结果](../inventories/analysis_ready_v1/fingerprint_summary.json)）。
- 修补：GSE125449的患者键，按样本重新赋值（只改obs，X不动；`patch_obs_v1.py`）。

## 组成

### composition_by_batch

| data_batch | cells | units | patients | sources |
|---|---|---|---|---|
| existing_analysis_set_v2 | 2,871,696 | 748 | 611 | 40 |
| supplement_v1 | 946,346 | 272 | 182 | 13 |

### composition_by_cancer

| analysis_main | cells | units | patients | sources |
|---|---|---|---|---|
| NSCLC | 548,410 | 175 | 142 | 8 |
| Breast cancer | 478,336 | 73 | 64 | 2 |
| Pancreatic ductal adenocarcinoma | 301,865 | 65 | 64 | 3 |
| Esophageal cancer | 300,699 | 136 | 75 | 3 |
| Diffuse glioma | 227,269 | 44 | 18 | 1 |
| Renal cell carcinoma | 221,502 | 48 | 32 | 3 |
| Bladder cancer | 220,384 | 18 | 14 | 2 |
| Colorectal cancer | 208,818 | 125 | 119 | 5 |
| Prostate cancer | 171,691 | 69 | 38 | 2 |
| Osteosarcoma | 168,322 | 23 | 23 | 3 |
| Gastric cancer | 143,921 | 41 | 41 | 2 |
| Papillary thyroid carcinoma | 136,049 | 31 | 22 | 3 |
| Hepatocellular carcinoma | 119,028 | 42 | 27 | 4 |
| Uveal melanoma | 97,704 | 17 | 17 | 2 |
| Cervical cancer | 94,947 | 26 | 18 | 2 |
| Acral melanoma | 83,997 | 16 | 15 | 2 |
| Pancreatic neuroendocrine tumor | 83,489 | 17 | 17 | 1 |
| Nasopharyngeal carcinoma | 82,604 | 10 | 10 | 1 |
| Intrahepatic cholangiocarcinoma | 52,317 | 24 | 17 | 3 |
| Basal cell carcinoma | 31,542 | 4 | 4 | 1 |
| Ovarian cancer | 22,696 | 7 | 7 | 1 |
| Renal medullary carcinoma | 22,452 | 9 | 9 | 1 |

### composition_by_sampling_site

| sampling_site | cells | units | patients | sources |
|---|---|---|---|---|
| metastatic_site | 351,815 | 73 | 71 | 13 |
| primary_organ_site | 3,065,433 | 867 | 670 | 48 |
| tumour_thrombus | 5,746 | 2 | 2 | 1 |
| unspecified | 395,048 | 78 | 67 | 5 |

### composition_by_treatment_status

| treatment_status | cells | units | patients | sources |
|---|---|---|---|---|
| on_or_after_treatment | 434,901 | 82 | 52 | 11 |
| other | 20,858 | 4 | 2 | 1 |
| unknown | 1,235,136 | 442 | 324 | 22 |
| untreated_or_pretreatment | 2,127,147 | 492 | 428 | 31 |

### composition_by_fresh_frozen

| fresh_frozen | cells | units | patients | sources |
|---|---|---|---|---|
| cryopreserved_or_frozen | 192,966 | 69 | 50 | 3 |
| fresh | 2,590,492 | 679 | 506 | 35 |
| fresh_or_cryopreserved | 30,628 | 22 | 22 | 1 |
| unknown | 1,003,956 | 250 | 215 | 15 |

### composition_by_technology

| technology | cells | units | patients | sources |
|---|---|---|---|---|
| 10x | 682,607 | 195 | 167 | 11 |
| 10x 5' | 104,912 | 60 | 60 | 1 |
| Microwell-seq | 7,987 | 8 | 8 | 1 |
| Singleron GEXSCOPE | 150,840 | 9 | 5 | 1 |
| scRNA (see source registry) | 2,871,696 | 748 | 611 | 40 |

### composition_by_capture_scope

| capture_scope | cells | units | patients | sources |
|---|---|---|---|---|
| CD45neg_sorted_fraction_recombined | 83,486 | 58 | 58 | 1 |
| CD45pos_sorted_fraction_recombined | 104,912 | 60 | 60 | 1 |
| source_gate_v7_qualified | 2,871,696 | 748 | 611 | 40 |
| whole_tissue | 734,676 | 148 | 116 | 11 |
| whole_tissue_viable_cell_FACS | 23,272 | 6 | 6 | 1 |

### composition_by_qc_flag

| qc_flag | cells | units | patients | sources |
|---|---|---|---|---|
| floor_cell_calling_on_unfiltered_matrix | 150,840 | 9 | 5 | 1 |
| inherited_fit_cell | 2,871,696 | 748 | 611 | 40 |
| low_depth_platform | 7,987 | 8 | 8 | 1 |
| qc_v1_pass | 599,121 | 137 | 109 | 10 |
| sorted_fractions_recombined | 188,398 | 118 | 60 | 1 |

### composition_by_source

| source | data_batch | analysis_main | cells | units | patients | sources |
|---|---|---|---|---|---|---|
| GSE125449 | supplement_v1 | Hepatocellular carcinoma | 3,525 | 8 | 8 | 1 |
| GSE125449 | supplement_v1 | Intrahepatic cholangiocarcinoma | 5,705 | 10 | 10 | 1 |
| GSE127465 | existing_analysis_set_v2 | NSCLC | 21,286 | 18 | 7 | 1 |
| GSE131907 | existing_analysis_set_v2 | NSCLC | 102,685 | 32 | 32 | 1 |
| GSE132465 | existing_analysis_set_v2 | Colorectal cancer | 33,874 | 23 | 23 | 1 |
| GSE138665 | existing_analysis_set_v2 | Uveal melanoma | 10,624 | 6 | 6 | 1 |
| GSE138709 | existing_analysis_set_v2 | Intrahepatic cholangiocarcinoma | 13,144 | 5 | 4 | 1 |
| GSE139829 | existing_analysis_set_v2 | Uveal melanoma | 87,080 | 11 | 11 | 1 |
| GSE141526 | existing_analysis_set_v2 | Basal cell carcinoma | 31,542 | 4 | 4 | 1 |
| GSE148071 | existing_analysis_set_v2 | NSCLC | 70,331 | 42 | 42 | 1 |
| GSE149614 | supplement_v1 | Hepatocellular carcinoma | 41,933 | 13 | 10 | 1 |
| GSE152048 | existing_analysis_set_v2 | Osteosarcoma | 101,497 | 11 | 11 | 1 |
| GSE155698 | existing_analysis_set_v2 | Pancreatic ductal adenocarcinoma | 32,043 | 16 | 15 | 1 |
| GSE159115 | existing_analysis_set_v2 | Renal cell carcinoma | 14,076 | 8 | 8 | 1 |
| GSE160269 | supplement_v1 | Esophageal cancer | 188,398 | 118 | 60 | 1 |
| GSE161529 | existing_analysis_set_v2 | Breast cancer | 225,597 | 41 | 32 | 1 |
| GSE162025 | existing_analysis_set_v2 | Nasopharyngeal carcinoma | 82,604 | 10 | 10 | 1 |
| GSE162454 | supplement_v1 | Osteosarcoma | 43,553 | 6 | 6 | 1 |
| GSE163203 | supplement_v1 | Papillary thyroid carcinoma | 7,987 | 8 | 8 | 1 |
| GSE166555 | existing_analysis_set_v2 | Colorectal cancer | 27,156 | 13 | 12 | 1 |
| GSE176031 | existing_analysis_set_v2 | Prostate cancer | 16,273 | 26 | 14 | 1 |
| GSE178341 | existing_analysis_set_v2 | Colorectal cancer | 96,761 | 67 | 62 | 1 |
| GSE178481 | existing_analysis_set_v2 | Renal cell carcinoma | 77,764 | 15 | 9 | 1 |
| GSE179994 | existing_analysis_set_v2 | NSCLC | 140,640 | 46 | 36 | 1 |
| GSE182109 | existing_analysis_set_v2 | Diffuse glioma | 227,269 | 44 | 18 | 1 |
| GSE183904 | existing_analysis_set_v2 | Gastric cancer | 82,969 | 29 | 29 | 1 |
| GSE184362 | existing_analysis_set_v2 | Papillary thyroid carcinoma | 102,666 | 17 | 11 | 1 |
| GSE184880 | existing_analysis_set_v2 | Ovarian cancer | 22,696 | 7 | 7 | 1 |
| GSE188711 | existing_analysis_set_v2 | Colorectal cancer | 23,019 | 6 | 6 | 1 |
| GSE189357 | existing_analysis_set_v2 | NSCLC | 84,173 | 9 | 9 | 1 |
| GSE189889 | existing_analysis_set_v2 | Acral melanoma | 26,839 | 9 | 9 | 1 |
| GSE189903 | supplement_v1 | Hepatocellular carcinoma | 43,910 | 16 | 4 | 1 |
| GSE189903 | supplement_v1 | Intrahepatic cholangiocarcinoma | 33,468 | 9 | 3 | 1 |
| GSE191288 | supplement_v1 | Papillary thyroid carcinoma | 25,396 | 6 | 3 | 1 |
| GSE197461 | existing_analysis_set_v2 | Cervical cancer | 27,224 | 8 | 8 | 1 |
| GSE200997 | existing_analysis_set_v2 | Colorectal cancer | 28,008 | 16 | 16 | 1 |
| GSE201347 | existing_analysis_set_v2 | Gastric cancer | 60,952 | 12 | 12 | 1 |
| GSE203067 | supplement_v1 | Esophageal cancer | 51,379 | 7 | 4 | 1 |
| GSE205013 | supplement_v1 | Pancreatic ductal adenocarcinoma | 239,194 | 27 | 27 | 1 |
| GSE215120 | supplement_v1 | Acral melanoma | 57,158 | 7 | 6 | 1 |
| GSE222315 | existing_analysis_set_v2 | Bladder cancer | 69,544 | 9 | 9 | 1 |
| GSE242230 | supplement_v1 | Pancreatic ductal adenocarcinoma | 30,628 | 22 | 22 | 1 |
| GSE253013 | existing_analysis_set_v2 | NSCLC | 18,614 | 11 | 3 | 1 |
| GSE256136 | existing_analysis_set_v2 | Pancreatic neuroendocrine tumor | 83,489 | 17 | 17 | 1 |
| GSE256326 | existing_analysis_set_v2 | Renal medullary carcinoma | 22,452 | 9 | 9 | 1 |
| GSE273127 | existing_analysis_set_v2 | Esophageal cancer | 60,922 | 11 | 11 | 1 |
| GSE274229 | existing_analysis_set_v2 | Prostate cancer | 155,418 | 43 | 24 | 1 |
| GSE274934 | existing_analysis_set_v2 | NSCLC | 87,641 | 9 | 9 | 1 |
| GSE277524 | supplement_v1 | Bladder cancer | 150,840 | 9 | 5 | 1 |
| GSE288479 | existing_analysis_set_v2 | NSCLC | 23,040 | 8 | 4 | 1 |
| GSE294109 | existing_analysis_set_v2 | Renal cell carcinoma | 129,662 | 25 | 15 | 1 |
| GSE297041 | existing_analysis_set_v2 | Cervical cancer | 67,723 | 18 | 10 | 1 |
| GSE299023 | supplement_v1 | Osteosarcoma | 23,272 | 6 | 6 | 1 |
| GSE299340 | existing_analysis_set_v2 | Hepatocellular carcinoma | 29,660 | 5 | 5 | 1 |
| GSE305078 | existing_analysis_set_v2 | Breast cancer | 252,739 | 32 | 32 | 1 |

## 补充之后的癌种覆盖

- 有≥2个独立来源的癌种由9个增加到**16个**：NSCLC 8，结直肠5，HCC 4，PDAC、食管、肾细胞癌、骨肉瘤、PTC、ICC各3，乳腺、膀胱、前列腺、胃、葡萄膜黑色素瘤、宫颈、肢端黑色素瘤各2。
- 仍只有1个来源：弥漫性胶质瘤、PanNET、鼻咽癌、BCC、卵巢癌、肾髓质癌。这些癌种的结论不能和来源效应分开。
- **精原细胞瘤已不在数据集中**（唯一来源GSE228501已排除，候选中没有合格来源）。

## 使用注意事项

1. **跨来源分析只用共同基因**：`var.measured_in_all_sources`（12,498个）；或者按分析需要，用`varm['measured_in_source']`选择。不要把未测到的基因当作0。
2. **现有部分与补充部分的处理不完全一样**：现有细胞沿用原存储的QC、双细胞判定和细胞筛选（`inherited_fit_cell`），部分来源在原存储中有环境RNA校正层，但这里的X**没有**扣除；补充部分用QC v1，没有做环境RNA估计。对QC或环境RNA敏感的分析，应按`data_batch`、`qc_flag`做敏感性分析。
3. **组成分析要排除或单独处理**：GSE160269（分选后合并，免疫与非免疫的比例不代表组织构成）；GSE299023（活细胞FACS）。
4. **低质量或细胞识别不确定的来源**：GSE163203（低深度），GSE277524（只用QC下限识别细胞，可能混有环境RNA液滴）。建议主分析中设为敏感性项，或者对它们收紧阈值。
5. **患者是来源内的键**：跨研究的同一患者无法排除。补充部分`patient_id`来自GEO样本记录。按患者统计时用`patient_id`，不要用`unit_id`。
6. **治疗与新鲜/冻存仍有未知**：`treatment_status=unknown` 1,235,136个细胞（22个来源），`fresh_frozen=unknown` 1,003,956个细胞。临床相关分析只能使用已知的子集。
7. **现有部分的`technology`只写了"scRNA"**：具体平台见来源登记；需要按平台分层时要补这一字段。
8. 这是数据层，不是图谱：没有标准化、HVG、PCA、整合、聚类或注释。现有部分的旧注释和嵌入可以通过`origin_store` + `origin_native_row`回到原存储读取。
