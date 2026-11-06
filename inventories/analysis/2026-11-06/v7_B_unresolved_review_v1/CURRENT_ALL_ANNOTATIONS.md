# 当前整个V7注释：完整标签与工作状态

快照生产完成UTC：2026-11-06T05:55:44.678902+00:00；独立审计PASS：2026-11-06T05:57:49.171116+00:00。冻结V7全部 **4,676,787** 行、**106** 个原始标签完整保留；本视图不是新TME图谱。

[B逐群分析与六张图](README.md) · [完整可筛选Excel](V7_current_annotation_20261106.xlsx) · [全部冻结标签原始表](global/stage_g_consolidated_annotation/all_106_frozen_labels_current_status.tsv) · [当前标签/状态完整表](global/stage_g_consolidated_annotation/current_annotation_label_counts.tsv) · [旧→当前交叉表](global/stage_g_consolidated_annotation/frozen_to_current_crosswalk.tsv)

## 覆盖范围和当前QC

| 口径 | 细胞数 | 含义 |
|---|---|---|
| 冻结原始总数 | 4,676,787 | 保留原标签、置信度、strict历史资格和CNV字段 |
| 当前QC保留 | 4,668,312 | 冻结轴上的主分析保留视图 |
| 当前QC排除 | 8,237 | 本轮及前轮累计操作性排除；原记录不删除 |
| 当前QC暂存 | 238 | 不进入主分析，保留语境复查 |
| 本周期实际重建目标 | 423,375 | Mast/T6、四稀有谱系及B模型目标；集合不重叠 |
| 大类工作／候选记录 | 322,270 | 包括改到其他父类及HSPC-like等候选，不是发表通过数 |
| 重建后仍需身份复查 | 101,105 | 父类混杂、低信息或正向竞争程序等 |
| QC保留但未重建 | 4,244,937 | 沿用旧V7，不冒充本周期重新接受 |
| 本周期最终细分新接受 | 0 | 工作身份、亚型候选和状态背景均不自动升级为最终细分 |

全轴有 **53** 个当前父类／待核类别、**142** 个标签值；它们包含保留的历史/QC记录和待核占位，不是发现了同等数量的新细胞类型。实际复核 **105** 个诊断分区，不是105种亚型。B参考8,565个只用于定位，未覆盖任何旧身份。

## 当前全部父类及待核类别

下表按当前QC保留数排序。诊断模型有自己的候选召回范围，不能把本轮模型中某类的人数当泛癌全部该类丰度。父类改变没有自动覆盖冻结良恶性结论。

| 当前父类／待核类别 | 全部记录 | QC保留 | 排除 | 暂存 | 已重建 | 大类工作 | 身份复查 | 沿用未重建 |
|---|---|---|---|---|---|---|---|---|
| Tcell | 1,343,029 | 1,342,181 | 786 | 62 | 163,741 | 163,691 | 50 | 1,178,440 |
| Epithelial | 865,585 | 865,581 | 3 | 1 | 1,031 | 1,008 | 23 | 864,550 |
| Bcell | 436,343 | 431,042 | 5,301 | 0 | 84,884 | 84,447 | 437 | 346,158 |
| Macrophage | 419,057 | 419,051 | 6 | 0 | 0 | 0 | 0 | 419,051 |
| Fibroblast | 325,496 | 325,496 | 0 | 0 | 0 | 0 | 0 | 325,496 |
| NK | 193,010 | 193,008 | 2 | 0 | 0 | 0 | 0 | 193,008 |
| Plasma | 183,888 | 183,887 | 1 | 0 | 0 | 0 | 0 | 183,887 |
| Myeloid_boundary | 126,209 | 126,208 | 1 | 0 | 0 | 0 | 0 | 126,208 |
| Endothelial | 107,108 | 107,101 | 7 | 0 | 512 | 508 | 4 | 106,589 |
| Monocyte | 98,427 | 98,427 | 0 | 0 | 0 | 0 | 0 | 98,427 |
| T_NK_ILC_boundary | 97,274 | 97,186 | 88 | 0 | 0 | 0 | 0 | 97,186 |
| Neural | 88,118 | 88,116 | 2 | 0 | 0 | 0 | 0 | 88,116 |
| Ambiguous | 75,674 | 74,890 | 783 | 1 | 0 | 0 | 0 | 74,890 |
| Unknown | 53,707 | 53,707 | 0 | 0 | 29 | 0 | 29 | 53,678 |
| Mast | 45,036 | 43,677 | 1,186 | 173 | 38,923 | 38,923 | 0 | 4,754 |
| Granulocyte | 33,827 | 33,827 | 0 | 0 | 0 | 0 | 0 | 33,827 |
| Parent_unresolved | 30,389 | 30,389 | 0 | 0 | 30,389 | 0 | 30,389 | 0 |
| IL7R_lymphoid_T_ILC_boundary | 29,406 | 29,406 | 0 | 0 | 29,406 | 0 | 29,406 | 0 |
| Melanocytic | 24,625 | 24,625 | 0 | 0 | 608 | 580 | 28 | 24,017 |
| Dendritic | 19,646 | 19,645 | 1 | 0 | 0 | 0 | 0 | 19,645 |
| Myeloid | 17,679 | 17,679 | 0 | 0 | 17,679 | 16,907 | 772 | 0 |
| Nonhematopoietic_parent_review | 14,540 | 14,540 | 0 | 0 | 14,540 | 0 | 14,540 | 0 |
| HSPC_like | 12,546 | 12,546 | 0 | 0 | 12,546 | 12,436 | 110 | 0 |
| T_candidate_low_information | 5,704 | 5,704 | 0 | 0 | 5,704 | 0 | 5,704 | 0 |
| Related_lineage_boundary | 4,970 | 4,970 | 0 | 0 | 0 | 0 | 0 | 4,970 |
| B_ASC_boundary | 4,964 | 4,964 | 0 | 0 | 4,964 | 0 | 4,964 | 0 |
| ILC_T_NK_boundary | 4,288 | 4,288 | 0 | 0 | 4,288 | 0 | 4,288 | 0 |
| Melanocytic_Mural_context | 2,432 | 2,432 | 0 | 0 | 2,432 | 0 | 2,432 | 0 |
| Mural | 2,393 | 2,393 | 0 | 0 | 2,393 | 2,380 | 13 | 0 |
| KIT_RORC_context | 1,867 | 1,867 | 0 | 0 | 1,867 | 0 | 1,867 | 0 |
| Erythroid_candidate | 1,815 | 1,815 | 0 | 0 | 1,815 | 0 | 1,815 | 0 |
| Hematopoietic_differentiation_boundary | 1,228 | 1,227 | 1 | 0 | 0 | 0 | 0 | 1,227 |
| Mixed_Ery_Epithelial | 1,045 | 1,045 | 0 | 0 | 1,045 | 0 | 1,045 | 0 |
| Stromal_boundary | 813 | 813 | 0 | 0 | 0 | 0 | 0 | 813 |
| Erythroid | 809 | 794 | 14 | 1 | 794 | 789 | 5 | 0 |
| T_with_competing_RNA_review | 662 | 662 | 0 | 0 | 662 | 0 | 662 | 0 |
| Endothelial_like | 596 | 596 | 0 | 0 | 596 | 0 | 596 | 0 |
| B_ASC_context | 371 | 371 | 0 | 0 | 371 | 0 | 371 | 0 |
| Hepatocyte_like | 357 | 357 | 0 | 0 | 357 | 338 | 19 | 0 |
| Neural_neuroendocrine_like | 351 | 351 | 0 | 0 | 351 | 0 | 351 | 0 |
| Mixed_neural_Platelet | 309 | 309 | 0 | 0 | 309 | 0 | 309 | 0 |
| Myeloid_lymphoid_context | 272 | 272 | 0 | 0 | 272 | 0 | 272 | 0 |
| ILC_like_candidate | 166 | 166 | 0 | 0 | 166 | 0 | 166 | 0 |
| Platelet_Megakaryocytic | 149 | 149 | 0 | 0 | 149 | 144 | 5 | 0 |
| Epithelial_Ig_context | 122 | 122 | 0 | 0 | 122 | 0 | 122 | 0 |
| Mixed_myeloid_Platelet | 120 | 120 | 0 | 0 | 120 | 0 | 120 | 0 |
| Immature_myeloid_like | 120 | 120 | 0 | 0 | 120 | 119 | 1 | 0 |
| Stromal_Epithelial_context | 84 | 84 | 0 | 0 | 84 | 0 | 84 | 0 |
| Mixed_T_Platelet | 58 | 58 | 0 | 0 | 58 | 0 | 58 | 0 |
| Acinar_like | 48 | 48 | 0 | 0 | 48 | 0 | 48 | 0 |
| Platelet | 3 | 0 | 3 | 0 | 0 | 0 | 0 | 0 |
| B_Plasma_boundary | 26 | 0 | 26 | 0 | 0 | 0 | 0 | 0 |
| ILC | 26 | 0 | 26 | 0 | 0 | 0 | 0 | 0 |

## 当前全部标签

每个标签的QC保留、排除、暂存及工作状态完整列出。此表合并了候选亚型/群级背景的多条聚合；这些细节在Excel“当前标签”及完整TSV中逐条保留，不能由同一标签名推断所有细胞拥有同一状态。

| 当前父类／待核类别 | 标签值 | 全部记录 | QC保留 | 排除 | 暂存 | 大类工作 | 身份复查 | 沿用未重建 |
|---|---|---|---|---|---|---|---|---|
| Acinar_like | Acinar_like | 48 | 48 | 0 | 0 | 0 | 48 | 0 |
| Ambiguous | Ambiguous | 75,674 | 74,890 | 783 | 1 | 0 | 0 | 74,890 |
| B_ASC_boundary | B_ASC_identity_review | 4,964 | 4,964 | 0 | 0 | 0 | 4,964 | 0 |
| B_ASC_context | B_ASC_context | 371 | 371 | 0 | 0 | 0 | 371 | 0 |
| B_Plasma_boundary | B_Plasma_boundary | 26 | 0 | 26 | 0 | 0 | 0 | 0 |
| Bcell | Bcell_Naive_IGHD | 106,404 | 105,825 | 579 | 0 | 0 | 0 | 105,825 |
| Bcell | Bcell_Naive_TCL1A | 129,522 | 128,798 | 724 | 0 | 0 | 0 | 128,798 |
| Bcell | Bcell_broad | 84,884 | 84,884 | 0 | 0 | 84,447 | 437 | 0 |
| Bcell | Bcell_coarse_provisional | 79,965 | 76,706 | 3,259 | 0 | 0 | 0 | 76,706 |
| Bcell | Bcell_coarse_subtype_pending | 35,223 | 34,829 | 394 | 0 | 0 | 0 | 34,829 |
| Bcell | Bcell_unresolved | 25 | 0 | 25 | 0 | 0 | 0 | 0 |
| Bcell | Bcell_unresolved_CST3 | 268 | 0 | 268 | 0 | 0 | 0 | 0 |
| Bcell | Bcell_unresolved_DNAJB1 | 52 | 0 | 52 | 0 | 0 | 0 | 0 |
| Dendritic | Dendritic_coarse_provisional | 2,723 | 2,723 | 0 | 0 | 0 | 0 | 2,723 |
| Dendritic | Dendritic_coarse_subtype_pending | 272 | 272 | 0 | 0 | 0 | 0 | 272 |
| Dendritic | Dendritic_pDC_LILRA4 | 2,072 | 2,072 | 0 | 0 | 0 | 0 | 2,072 |
| Dendritic | Dendritic_unresolved | 4,697 | 4,696 | 1 | 0 | 0 | 0 | 4,696 |
| Dendritic | Dendritic_unresolved_CD83 | 4,881 | 4,881 | 0 | 0 | 0 | 0 | 4,881 |
| Dendritic | Dendritic_unresolved_SPI1 | 5,001 | 5,001 | 0 | 0 | 0 | 0 | 5,001 |
| Endothelial | Endothelial | 261 | 261 | 0 | 0 | 257 | 4 | 0 |
| Endothelial | Endothelial_Blood_PECAM1 | 21,746 | 21,745 | 1 | 0 | 0 | 0 | 21,745 |
| Endothelial | Endothelial_Venous_ACKR1_ACKR1 | 16,238 | 16,237 | 1 | 0 | 0 | 0 | 16,237 |
| Endothelial | Endothelial_broad | 251 | 251 | 0 | 0 | 251 | 0 | 0 |
| Endothelial | Endothelial_coarse_provisional | 4,449 | 4,449 | 0 | 0 | 0 | 0 | 4,449 |
| Endothelial | Endothelial_coarse_subtype_pending | 1,509 | 1,508 | 1 | 0 | 0 | 0 | 1,508 |
| Endothelial | Endothelial_unresolved | 41,379 | 41,377 | 2 | 0 | 0 | 0 | 41,377 |
| Endothelial | Endothelial_unresolved_INSR | 6,158 | 6,157 | 1 | 0 | 0 | 0 | 6,157 |
| Endothelial | Endothelial_unresolved_PLVAP | 15,117 | 15,116 | 1 | 0 | 0 | 0 | 15,116 |
| Endothelial_like | Endothelial_like | 596 | 596 | 0 | 0 | 0 | 596 | 0 |
| Epithelial | Epithelial | 1,031 | 1,031 | 0 | 0 | 1,008 | 23 | 0 |
| Epithelial | Epithelial_coarse_CNV_not_recomputed | 53,387 | 53,386 | 1 | 0 | 0 | 0 | 53,386 |
| Epithelial | Epithelial_coarse_provisional | 33,776 | 33,775 | 1 | 0 | 0 | 0 | 33,775 |
| Epithelial | Epithelial_malignancy_discordant | 229,472 | 229,472 | 0 | 0 | 0 | 0 | 229,472 |
| Epithelial | Epithelial_malignancy_not_evaluable | 163,427 | 163,426 | 0 | 1 | 0 | 0 | 163,426 |
| Epithelial | Malignant_Epithelial_high_specificity | 185,305 | 185,305 | 0 | 0 | 0 | 0 | 185,305 |
| Epithelial | Nonmalignant_Epithelial_CNV_negative | 199,187 | 199,186 | 1 | 0 | 0 | 0 | 199,186 |
| Epithelial_Ig_context | Epithelial_Ig_context | 122 | 122 | 0 | 0 | 0 | 122 | 0 |
| Erythroid | Erythroid | 794 | 794 | 0 | 0 | 789 | 5 | 0 |
| Erythroid | Erythroid_Cycling_Erythroid_TOP2A | 6 | 0 | 6 | 0 | 0 | 0 | 0 |
| Erythroid | Erythroid_coarse_provisional | 2 | 0 | 2 | 0 | 0 | 0 | 0 |
| Erythroid | Erythroid_coarse_subtype_pending | 2 | 0 | 2 | 0 | 0 | 0 | 0 |
| Erythroid | Erythroid_unresolved_PTPRC | 5 | 0 | 4 | 1 | 0 | 0 | 0 |
| Erythroid_candidate | Erythroid_candidate | 1,815 | 1,815 | 0 | 0 | 0 | 1,815 | 0 |
| Fibroblast | Fibroblast_Cycling_Fibroblast_TOP2A | 39,720 | 39,720 | 0 | 0 | 0 | 0 | 39,720 |
| Fibroblast | Fibroblast_ECM_CAF_COL1A1 | 48,829 | 48,829 | 0 | 0 | 0 | 0 | 48,829 |
| Fibroblast | Fibroblast_coarse_provisional | 26,215 | 26,215 | 0 | 0 | 0 | 0 | 26,215 |
| Fibroblast | Fibroblast_coarse_subtype_pending | 21,626 | 21,626 | 0 | 0 | 0 | 0 | 21,626 |
| Fibroblast | Fibroblast_iCAF_CXCL12 | 35,918 | 35,918 | 0 | 0 | 0 | 0 | 35,918 |
| Fibroblast | Fibroblast_iCAF_IL6 | 26,748 | 26,748 | 0 | 0 | 0 | 0 | 26,748 |
| Fibroblast | Fibroblast_myCAF_ACTA2 | 89,122 | 89,122 | 0 | 0 | 0 | 0 | 89,122 |
| Fibroblast | Fibroblast_unresolved | 37,318 | 37,318 | 0 | 0 | 0 | 0 | 37,318 |
| Granulocyte | Granulocyte_Neutrophil_mature_FCGR3B | 13,977 | 13,977 | 0 | 0 | 0 | 0 | 13,977 |
| Granulocyte | Granulocyte_coarse_provisional | 3,531 | 3,531 | 0 | 0 | 0 | 0 | 3,531 |
| Granulocyte | Granulocyte_coarse_subtype_pending | 10,779 | 10,779 | 0 | 0 | 0 | 0 | 10,779 |
| Granulocyte | Granulocyte_unresolved_S100A10 | 5,540 | 5,540 | 0 | 0 | 0 | 0 | 5,540 |
| HSPC_like | HSPC_like | 12,546 | 12,546 | 0 | 0 | 12,436 | 110 | 0 |
| Hematopoietic_differentiation_boundary | Hematopoietic_differentiation_boundary | 1,228 | 1,227 | 1 | 0 | 0 | 0 | 1,227 |
| Hepatocyte_like | Hepatocyte_like | 357 | 357 | 0 | 0 | 338 | 19 | 0 |
| IL7R_lymphoid_T_ILC_boundary | IL7R_lymphoid_T_ILC_boundary | 29,406 | 29,406 | 0 | 0 | 0 | 29,406 | 0 |
| ILC | ILC_coarse_provisional | 18 | 0 | 18 | 0 | 0 | 0 | 0 |
| ILC | ILC_coarse_subtype_pending | 6 | 0 | 6 | 0 | 0 | 0 | 0 |
| ILC | ILC_unresolved_CD63 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| ILC | ILC_unresolved_SRGN | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| ILC_T_NK_boundary | ILC_T_NK_boundary | 4,288 | 4,288 | 0 | 0 | 0 | 4,288 | 0 |
| ILC_like_candidate | ILC_like_candidate | 166 | 166 | 0 | 0 | 0 | 166 | 0 |
| Immature_myeloid_like | Immature_myeloid_like | 120 | 120 | 0 | 0 | 119 | 1 | 0 |
| KIT_RORC_context | KIT_RORC_context | 1,867 | 1,867 | 0 | 0 | 0 | 1,867 | 0 |
| Macrophage | Macrophage_Cycling_Macrophage_MKI67 | 16,129 | 16,129 | 0 | 0 | 0 | 0 | 16,129 |
| Macrophage | Macrophage_Inflammatory_IL1B_IL1B | 140,987 | 140,985 | 2 | 0 | 0 | 0 | 140,985 |
| Macrophage | Macrophage_SPP1_TREM2_SPP1 | 70,223 | 70,223 | 0 | 0 | 0 | 0 | 70,223 |
| Macrophage | Macrophage_coarse_provisional | 11,237 | 11,237 | 0 | 0 | 0 | 0 | 11,237 |
| Macrophage | Macrophage_coarse_subtype_pending | 17,536 | 17,534 | 2 | 0 | 0 | 0 | 17,534 |
| Macrophage | Macrophage_unresolved | 61,309 | 61,309 | 0 | 0 | 0 | 0 | 61,309 |
| Macrophage | Macrophage_unresolved_CD52 | 69,588 | 69,587 | 1 | 0 | 0 | 0 | 69,587 |
| Macrophage | Macrophage_unresolved_SERPINE2 | 32,048 | 32,047 | 1 | 0 | 0 | 0 | 32,047 |
| Mast | Mast | 38,923 | 38,923 | 0 | 0 | 38,923 | 0 | 0 |
| Mast | Mast_Activated_IL1RL1_IL1RL1 | 645 | 404 | 220 | 21 | 0 | 0 | 404 |
| Mast | Mast_Activated_IL1RL1_LTC4S | 1,585 | 1,024 | 484 | 77 | 0 | 0 | 1,024 |
| Mast | Mast_coarse_provisional | 967 | 683 | 260 | 24 | 0 | 0 | 683 |
| Mast | Mast_coarse_subtype_pending | 595 | 461 | 125 | 9 | 0 | 0 | 461 |
| Mast | Mast_unresolved | 2,321 | 2,182 | 97 | 42 | 0 | 0 | 2,182 |
| Melanocytic | Melanocytic | 608 | 608 | 0 | 0 | 580 | 28 | 0 |
| Melanocytic | Melanocytic_coarse_subtype_pending | 24,017 | 24,017 | 0 | 0 | 0 | 0 | 24,017 |
| Melanocytic_Mural_context | Melanocytic_Mural_context | 2,432 | 2,432 | 0 | 0 | 0 | 2,432 | 0 |
| Mixed_Ery_Epithelial | Mixed_Ery_Epithelial | 1,045 | 1,045 | 0 | 0 | 0 | 1,045 | 0 |
| Mixed_T_Platelet | Mixed_T_Platelet | 58 | 58 | 0 | 0 | 0 | 58 | 0 |
| Mixed_myeloid_Platelet | Mixed_myeloid_Platelet | 120 | 120 | 0 | 0 | 0 | 120 | 0 |
| Mixed_neural_Platelet | Mixed_neural_Platelet | 309 | 309 | 0 | 0 | 0 | 309 | 0 |
| Monocyte | Monocyte_Classical_FCN1_FCN1 | 16,626 | 16,626 | 0 | 0 | 0 | 0 | 16,626 |
| Monocyte | Monocyte_Inflammatory_IL1B_IL1B | 40,731 | 40,731 | 0 | 0 | 0 | 0 | 40,731 |
| Monocyte | Monocyte_Nonclassical_FCGR3A_FCGR3A | 11,681 | 11,681 | 0 | 0 | 0 | 0 | 11,681 |
| Monocyte | Monocyte_coarse_provisional | 1,912 | 1,912 | 0 | 0 | 0 | 0 | 1,912 |
| Monocyte | Monocyte_coarse_subtype_pending | 11,782 | 11,782 | 0 | 0 | 0 | 0 | 11,782 |
| Monocyte | Monocyte_unresolved | 15,695 | 15,695 | 0 | 0 | 0 | 0 | 15,695 |
| Mural | Mural | 2,393 | 2,393 | 0 | 0 | 2,380 | 13 | 0 |
| Myeloid | Myeloid | 948 | 948 | 0 | 0 | 881 | 67 | 0 |
| Myeloid | Myeloid_broad | 16,731 | 16,731 | 0 | 0 | 16,026 | 705 | 0 |
| Myeloid_boundary | Myeloid_boundary | 126,209 | 126,208 | 1 | 0 | 0 | 0 | 126,208 |
| Myeloid_lymphoid_context | Myeloid_lymphoid_context | 272 | 272 | 0 | 0 | 0 | 272 | 0 |
| NK | NK_FCGR3A_cytotoxic_FCGR3A | 60,881 | 60,881 | 0 | 0 | 0 | 0 | 60,881 |
| NK | NK_XCL1_CD56bright_XCL1 | 37,526 | 37,526 | 0 | 0 | 0 | 0 | 37,526 |
| NK | NK_coarse_provisional | 26,924 | 26,924 | 0 | 0 | 0 | 0 | 26,924 |
| NK | NK_coarse_subtype_pending | 15,389 | 15,388 | 1 | 0 | 0 | 0 | 15,388 |
| NK | NK_unresolved | 35,925 | 35,925 | 0 | 0 | 0 | 0 | 35,925 |
| NK | NK_unresolved_CST3 | 4,092 | 4,091 | 1 | 0 | 0 | 0 | 4,091 |
| NK | NK_unresolved_DNAJB1 | 12,273 | 12,273 | 0 | 0 | 0 | 0 | 12,273 |
| Neural | Neural_Neural_progenitor_SOX2 | 62,127 | 62,127 | 0 | 0 | 0 | 0 | 62,127 |
| Neural | Neural_coarse_provisional | 3,006 | 3,006 | 0 | 0 | 0 | 0 | 3,006 |
| Neural | Neural_coarse_subtype_pending | 2,425 | 2,423 | 2 | 0 | 0 | 0 | 2,423 |
| Neural | Neural_unresolved | 20,560 | 20,560 | 0 | 0 | 0 | 0 | 20,560 |
| Neural_neuroendocrine_like | Neural_neuroendocrine_like | 351 | 351 | 0 | 0 | 0 | 351 | 0 |
| Nonhematopoietic_parent_review | Mixed_nonhematopoietic_review | 14,540 | 14,540 | 0 | 0 | 0 | 14,540 | 0 |
| Parent_unresolved | Cycling_parent_review | 3,595 | 3,595 | 0 | 0 | 0 | 3,595 | 0 |
| Parent_unresolved | Hepatic_RNA_parent_review | 136 | 136 | 0 | 0 | 0 | 136 | 0 |
| Parent_unresolved | Low_information_parent_review | 4,521 | 4,521 | 0 | 0 | 0 | 4,521 | 0 |
| Parent_unresolved | Mixed_parent_review | 21,577 | 21,577 | 0 | 0 | 0 | 21,577 | 0 |
| Parent_unresolved | Renal_hepatic_myeloid_RNA_review | 560 | 560 | 0 | 0 | 0 | 560 | 0 |
| Plasma | Plasma_coarse_provisional | 16,074 | 16,074 | 0 | 0 | 0 | 0 | 16,074 |
| Plasma | Plasma_coarse_subtype_pending | 4,416 | 4,416 | 0 | 0 | 0 | 0 | 4,416 |
| Plasma | Plasma_unresolved | 81,085 | 81,084 | 1 | 0 | 0 | 0 | 81,084 |
| Plasma | Plasma_unresolved_KRTCAP2 | 69,779 | 69,779 | 0 | 0 | 0 | 0 | 69,779 |
| Plasma | Plasma_unresolved_PFN1 | 12,534 | 12,534 | 0 | 0 | 0 | 0 | 12,534 |
| Platelet | Platelet_coarse | 2 | 0 | 2 | 0 | 0 | 0 | 0 |
| Platelet | Platelet_coarse_subtype_pending | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| Platelet_Megakaryocytic | Platelet_Megakaryocytic | 149 | 149 | 0 | 0 | 144 | 5 | 0 |
| Related_lineage_boundary | Related_lineage_boundary | 4,970 | 4,970 | 0 | 0 | 0 | 0 | 4,970 |
| Stromal_Epithelial_context | Stromal_Epithelial_context | 84 | 84 | 0 | 0 | 0 | 84 | 0 |
| Stromal_boundary | Stromal_boundary | 813 | 813 | 0 | 0 | 0 | 0 | 813 |
| T_NK_ILC_boundary | T_NK_ILC_boundary | 97,274 | 97,186 | 88 | 0 | 0 | 0 | 97,186 |
| T_candidate_low_information | T_apparatus_parent_review | 5,554 | 5,554 | 0 | 0 | 0 | 5,554 | 0 |
| T_candidate_low_information | T_candidate_low_information | 150 | 150 | 0 | 0 | 0 | 150 | 0 |
| T_with_competing_RNA_review | T_with_competing_RNA_review | 662 | 662 | 0 | 0 | 0 | 662 | 0 |
| Tcell | Tcell | 163,166 | 163,166 | 0 | 0 | 163,116 | 50 | 0 |
| Tcell | Tcell_CD4_naive_IL7R | 519,131 | 519,131 | 0 | 0 | 0 | 0 | 519,131 |
| Tcell | Tcell_CD8_cytotoxic_CCL5 | 496,996 | 496,995 | 1 | 0 | 0 | 0 | 496,995 |
| Tcell | Tcell_Cycling_T_MKI67 | 41,870 | 41,868 | 2 | 0 | 0 | 0 | 41,868 |
| Tcell | Tcell_broad | 575 | 575 | 0 | 0 | 575 | 0 | 0 |
| Tcell | Tcell_coarse_provisional | 19,572 | 19,572 | 0 | 0 | 0 | 0 | 19,572 |
| Tcell | Tcell_coarse_subtype_pending | 14,080 | 14,080 | 0 | 0 | 0 | 0 | 14,080 |
| Tcell | Tcell_unresolved_BATF | 845 | 0 | 783 | 62 | 0 | 0 | 0 |
| Tcell | Tcell_unresolved_NKG7 | 86,794 | 86,794 | 0 | 0 | 0 | 0 | 86,794 |
| Unknown | Unknown | 53,707 | 53,707 | 0 | 0 | 0 | 29 | 53,678 |

## 全部106个冻结V7标签的保留情况

以下保留所有原始标签，人数从原V7读取，不删除已排除细胞的历史。历史strict细分资格只反映旧合同；原置信度、大类原因和细分方法的完整计数字典在Excel“冻结V7”及原始TSV中，不能把其“fine”或“strict”直接解释成本轮发表终验。

| 冻结父类 | 冻结标签 | 原人数 | 当前保留 | 排除 | 暂存 | 本周期重建 | 未重建保留 | 旧strict细分资格 |
|---|---|---|---|---|---|---|---|---|
| Ambiguous | Ambiguous | 76,072 | 75,288 | 783 | 1 | 398 | 74,890 | 0 |
| B_Plasma_boundary | B_Plasma_boundary | 4,839 | 4,813 | 26 | 0 | 4,813 | 0 | 0 |
| Bcell | Bcell_Naive_IGHD | 106,404 | 105,825 | 579 | 0 | 0 | 105,825 | 106,404 |
| Bcell | Bcell_Naive_TCL1A | 129,522 | 128,798 | 724 | 0 | 0 | 128,798 | 129,522 |
| Bcell | Bcell_coarse_provisional | 79,965 | 76,706 | 3,259 | 0 | 0 | 76,706 | 0 |
| Bcell | Bcell_coarse_subtype_pending | 35,223 | 34,829 | 394 | 0 | 0 | 34,829 | 0 |
| Bcell | Bcell_unresolved | 42,344 | 42,319 | 25 | 0 | 42,319 | 0 | 0 |
| Bcell | Bcell_unresolved_CST3 | 55,409 | 55,141 | 268 | 0 | 55,141 | 0 | 0 |
| Bcell | Bcell_unresolved_DNAJB1 | 55,667 | 55,615 | 52 | 0 | 55,615 | 0 | 0 |
| Dendritic | Dendritic_coarse_provisional | 2,723 | 2,723 | 0 | 0 | 0 | 2,723 | 0 |
| Dendritic | Dendritic_coarse_subtype_pending | 272 | 272 | 0 | 0 | 0 | 272 | 0 |
| Dendritic | Dendritic_pDC_LILRA4 | 2,072 | 2,072 | 0 | 0 | 0 | 2,072 | 2,072 |
| Dendritic | Dendritic_unresolved | 4,697 | 4,696 | 1 | 0 | 0 | 4,696 | 0 |
| Dendritic | Dendritic_unresolved_CD83 | 4,881 | 4,881 | 0 | 0 | 0 | 4,881 | 0 |
| Dendritic | Dendritic_unresolved_SPI1 | 5,001 | 5,001 | 0 | 0 | 0 | 5,001 | 0 |
| Endothelial | Endothelial_Blood_PECAM1 | 21,746 | 21,745 | 1 | 0 | 0 | 21,745 | 21,746 |
| Endothelial | Endothelial_Venous_ACKR1_ACKR1 | 16,238 | 16,237 | 1 | 0 | 0 | 16,237 | 16,238 |
| Endothelial | Endothelial_coarse_provisional | 4,449 | 4,449 | 0 | 0 | 0 | 4,449 | 0 |
| Endothelial | Endothelial_coarse_subtype_pending | 1,509 | 1,508 | 1 | 0 | 0 | 1,508 | 0 |
| Endothelial | Endothelial_unresolved | 41,379 | 41,377 | 2 | 0 | 0 | 41,377 | 0 |
| Endothelial | Endothelial_unresolved_INSR | 6,158 | 6,157 | 1 | 0 | 0 | 6,157 | 0 |
| Endothelial | Endothelial_unresolved_PLVAP | 15,117 | 15,116 | 1 | 0 | 0 | 15,116 | 0 |
| Epithelial | Epithelial_coarse_CNV_not_recomputed | 53,387 | 53,386 | 1 | 0 | 0 | 53,386 | 0 |
| Epithelial | Epithelial_coarse_provisional | 33,776 | 33,775 | 1 | 0 | 0 | 33,775 | 0 |
| Epithelial | Epithelial_malignancy_discordant | 229,472 | 229,472 | 0 | 0 | 0 | 229,472 | 0 |
| Epithelial | Epithelial_malignancy_not_evaluable | 163,427 | 163,426 | 0 | 1 | 0 | 163,426 | 0 |
| Epithelial | Malignant_Epithelial_high_specificity | 185,305 | 185,305 | 0 | 0 | 0 | 185,305 | 185,305 |
| Epithelial | Nonmalignant_Epithelial_CNV_negative | 199,187 | 199,186 | 1 | 0 | 0 | 199,186 | 199,187 |
| Erythroid | Erythroid_Cycling_Erythroid_TOP2A | 2,354 | 2,348 | 6 | 0 | 2,348 | 0 | 2,354 |
| Erythroid | Erythroid_coarse_provisional | 1,196 | 1,194 | 2 | 0 | 1,194 | 0 | 0 |
| Erythroid | Erythroid_coarse_subtype_pending | 1,649 | 1,647 | 2 | 0 | 1,647 | 0 | 0 |
| Erythroid | Erythroid_unresolved_PTPRC | 3,839 | 3,834 | 4 | 1 | 3,834 | 0 | 0 |
| Fibroblast | Fibroblast_Cycling_Fibroblast_TOP2A | 39,720 | 39,720 | 0 | 0 | 0 | 39,720 | 39,720 |
| Fibroblast | Fibroblast_ECM_CAF_COL1A1 | 48,829 | 48,829 | 0 | 0 | 0 | 48,829 | 48,829 |
| Fibroblast | Fibroblast_coarse_provisional | 26,215 | 26,215 | 0 | 0 | 0 | 26,215 | 0 |
| Fibroblast | Fibroblast_coarse_subtype_pending | 21,626 | 21,626 | 0 | 0 | 0 | 21,626 | 0 |
| Fibroblast | Fibroblast_iCAF_CXCL12 | 35,918 | 35,918 | 0 | 0 | 0 | 35,918 | 35,918 |
| Fibroblast | Fibroblast_iCAF_IL6 | 26,748 | 26,748 | 0 | 0 | 0 | 26,748 | 26,748 |
| Fibroblast | Fibroblast_myCAF_ACTA2 | 89,122 | 89,122 | 0 | 0 | 0 | 89,122 | 89,122 |
| Fibroblast | Fibroblast_unresolved | 37,318 | 37,318 | 0 | 0 | 0 | 37,318 | 0 |
| Granulocyte | Granulocyte_Neutrophil_mature_FCGR3B | 13,977 | 13,977 | 0 | 0 | 0 | 13,977 | 13,977 |
| Granulocyte | Granulocyte_coarse_provisional | 3,531 | 3,531 | 0 | 0 | 0 | 3,531 | 0 |
| Granulocyte | Granulocyte_coarse_subtype_pending | 10,779 | 10,779 | 0 | 0 | 0 | 10,779 | 0 |
| Granulocyte | Granulocyte_unresolved_S100A10 | 5,540 | 5,540 | 0 | 0 | 0 | 5,540 | 0 |
| HSPC | HSPC_coarse_provisional | 3,693 | 3,693 | 0 | 0 | 3,693 | 0 | 0 |
| HSPC | HSPC_coarse_subtype_pending | 907 | 907 | 0 | 0 | 907 | 0 | 0 |
| HSPC | HSPC_unresolved | 7,781 | 7,781 | 0 | 0 | 7,781 | 0 | 0 |
| HSPC | HSPC_unresolved_MT2A | 1,676 | 1,676 | 0 | 0 | 1,676 | 0 | 0 |
| Hematopoietic_differentiation_boundary | Hematopoietic_differentiation_boundary | 1,232 | 1,231 | 1 | 0 | 4 | 1,227 | 0 |
| ILC | ILC_coarse_provisional | 18,164 | 18,146 | 18 | 0 | 18,146 | 0 | 0 |
| ILC | ILC_coarse_subtype_pending | 10,491 | 10,485 | 6 | 0 | 10,485 | 0 | 0 |
| ILC | ILC_unresolved | 3,874 | 3,874 | 0 | 0 | 3,874 | 0 | 0 |
| ILC | ILC_unresolved_CD63 | 993 | 992 | 1 | 0 | 992 | 0 | 0 |
| ILC | ILC_unresolved_SMCHD1 | 6,164 | 6,164 | 0 | 0 | 6,164 | 0 | 0 |
| ILC | ILC_unresolved_SRGN | 5,228 | 5,227 | 1 | 0 | 5,227 | 0 | 0 |
| Macrophage | Macrophage_Cycling_Macrophage_MKI67 | 16,129 | 16,129 | 0 | 0 | 0 | 16,129 | 16,129 |
| Macrophage | Macrophage_Inflammatory_IL1B_IL1B | 140,987 | 140,985 | 2 | 0 | 0 | 140,985 | 140,987 |
| Macrophage | Macrophage_SPP1_TREM2_SPP1 | 70,223 | 70,223 | 0 | 0 | 0 | 70,223 | 70,223 |
| Macrophage | Macrophage_coarse_provisional | 11,237 | 11,237 | 0 | 0 | 0 | 11,237 | 0 |
| Macrophage | Macrophage_coarse_subtype_pending | 17,536 | 17,534 | 2 | 0 | 0 | 17,534 | 0 |
| Macrophage | Macrophage_unresolved | 61,309 | 61,309 | 0 | 0 | 0 | 61,309 | 0 |
| Macrophage | Macrophage_unresolved_CD52 | 69,588 | 69,587 | 1 | 0 | 0 | 69,587 | 0 |
| Macrophage | Macrophage_unresolved_SERPINE2 | 32,048 | 32,047 | 1 | 0 | 0 | 32,047 | 0 |
| Mast | Mast_Activated_IL1RL1_IL1RL1 | 10,942 | 10,701 | 220 | 21 | 10,297 | 404 | 10,942 |
| Mast | Mast_Activated_IL1RL1_LTC4S | 11,302 | 10,741 | 484 | 77 | 9,717 | 1,024 | 11,302 |
| Mast | Mast_coarse_provisional | 1,408 | 1,124 | 260 | 24 | 441 | 683 | 0 |
| Mast | Mast_coarse_subtype_pending | 825 | 691 | 125 | 9 | 230 | 461 | 0 |
| Mast | Mast_unresolved | 20,559 | 20,420 | 97 | 42 | 18,238 | 2,182 | 0 |
| Melanocytic | Melanocytic_coarse_subtype_pending | 24,017 | 24,017 | 0 | 0 | 0 | 24,017 | 0 |
| Monocyte | Monocyte_Classical_FCN1_FCN1 | 16,626 | 16,626 | 0 | 0 | 0 | 16,626 | 16,626 |
| Monocyte | Monocyte_Inflammatory_IL1B_IL1B | 40,731 | 40,731 | 0 | 0 | 0 | 40,731 | 40,731 |
| Monocyte | Monocyte_Nonclassical_FCGR3A_FCGR3A | 11,681 | 11,681 | 0 | 0 | 0 | 11,681 | 11,681 |
| Monocyte | Monocyte_coarse_provisional | 1,912 | 1,912 | 0 | 0 | 0 | 1,912 | 0 |
| Monocyte | Monocyte_coarse_subtype_pending | 11,782 | 11,782 | 0 | 0 | 0 | 11,782 | 0 |
| Monocyte | Monocyte_unresolved | 15,695 | 15,695 | 0 | 0 | 0 | 15,695 | 0 |
| Myeloid_boundary | Myeloid_boundary | 126,209 | 126,208 | 1 | 0 | 0 | 126,208 | 0 |
| NK | NK_FCGR3A_cytotoxic_FCGR3A | 60,881 | 60,881 | 0 | 0 | 0 | 60,881 | 60,881 |
| NK | NK_XCL1_CD56bright_XCL1 | 37,526 | 37,526 | 0 | 0 | 0 | 37,526 | 37,526 |
| NK | NK_coarse_provisional | 26,924 | 26,924 | 0 | 0 | 0 | 26,924 | 0 |
| NK | NK_coarse_subtype_pending | 15,389 | 15,388 | 1 | 0 | 0 | 15,388 | 0 |
| NK | NK_unresolved | 35,925 | 35,925 | 0 | 0 | 0 | 35,925 | 0 |
| NK | NK_unresolved_CST3 | 4,092 | 4,091 | 1 | 0 | 0 | 4,091 | 0 |
| NK | NK_unresolved_DNAJB1 | 12,273 | 12,273 | 0 | 0 | 0 | 12,273 | 0 |
| Neural | Neural_Neural_progenitor_SOX2 | 62,127 | 62,127 | 0 | 0 | 0 | 62,127 | 62,127 |
| Neural | Neural_coarse_provisional | 3,006 | 3,006 | 0 | 0 | 0 | 3,006 | 0 |
| Neural | Neural_coarse_subtype_pending | 2,425 | 2,423 | 2 | 0 | 0 | 2,423 | 0 |
| Neural | Neural_unresolved | 20,560 | 20,560 | 0 | 0 | 0 | 20,560 | 0 |
| Plasma | Plasma_coarse_provisional | 16,074 | 16,074 | 0 | 0 | 0 | 16,074 | 0 |
| Plasma | Plasma_coarse_subtype_pending | 4,416 | 4,416 | 0 | 0 | 0 | 4,416 | 0 |
| Plasma | Plasma_unresolved | 81,085 | 81,084 | 1 | 0 | 0 | 81,084 | 0 |
| Plasma | Plasma_unresolved_KRTCAP2 | 69,779 | 69,779 | 0 | 0 | 0 | 69,779 | 0 |
| Plasma | Plasma_unresolved_PFN1 | 12,534 | 12,534 | 0 | 0 | 0 | 12,534 | 0 |
| Platelet | Platelet_coarse | 1,130 | 1,128 | 2 | 0 | 1,128 | 0 | 1,130 |
| Platelet | Platelet_coarse_provisional | 46 | 46 | 0 | 0 | 46 | 0 | 0 |
| Platelet | Platelet_coarse_subtype_pending | 1,067 | 1,066 | 1 | 0 | 1,066 | 0 | 0 |
| Related_lineage_boundary | Related_lineage_boundary | 4,970 | 4,970 | 0 | 0 | 0 | 4,970 | 0 |
| Stromal_boundary | Stromal_boundary | 813 | 813 | 0 | 0 | 0 | 813 | 0 |
| T_NK_ILC_boundary | T_NK_ILC_boundary | 97,349 | 97,261 | 88 | 0 | 75 | 97,186 | 0 |
| Tcell | Tcell_CD4_naive_IL7R | 519,131 | 519,131 | 0 | 0 | 0 | 519,131 | 519,131 |
| Tcell | Tcell_CD8_cytotoxic_CCL5 | 496,996 | 496,995 | 1 | 0 | 0 | 496,995 | 496,996 |
| Tcell | Tcell_Cycling_T_MKI67 | 41,870 | 41,868 | 2 | 0 | 0 | 41,868 | 41,870 |
| Tcell | Tcell_coarse_provisional | 19,572 | 19,572 | 0 | 0 | 0 | 19,572 | 0 |
| Tcell | Tcell_coarse_subtype_pending | 14,080 | 14,080 | 0 | 0 | 0 | 14,080 | 0 |
| Tcell | Tcell_unresolved_BATF | 156,631 | 155,786 | 783 | 62 | 155,786 | 0 | 0 |
| Tcell | Tcell_unresolved_NKG7 | 86,794 | 86,794 | 0 | 0 | 0 | 86,794 | 0 |
| Unknown | Unknown | 53,771 | 53,771 | 0 | 0 | 93 | 53,678 | 0 |

## 全部已复核分区及候选程序

[全部105个诊断群完整清单](global/stage_g_consolidated_annotation/all_reviewed_cluster_inventory.tsv)和Excel“已复核群”包含模型内群号、命名marker、工作父类、候选细分、群级背景、目标/参考数、仍需复查人数及理由；此前Mast/T6的簇级背景和四稀有谱系的理由没有丢失。

[此前Mast/T6逐群图表](../v7_qc_clean_reannotation_v1/clean_refit_review_v2/README.md) · [此前Ery/Platelet/HSPC/ILC全部59群](../v7_other_unresolved_review_v1/ALL_CLUSTER_REVIEW.md)。166个ILC-like候选的67个type2/99个type3正向程序证据继续单独保留在[既有候选程序表](../v7_other_unresolved_review_v1/global/stage_g_ILC_candidate_followup/candidate_programme_summary.tsv)，这是程序候选，不是最终ILC2/3注释，不能因缺TCR转录本当Lin阴性。

## 冻结CNV与良恶性轴

七列冻结CNV/良恶性字段逐细胞原样保留，独立比较通过。本周期有144,308个工作父类改变的记录标为语境需重核；没有对其生成新的恶性/正常结论。Excel“冻结CNV轴”与[完整聚合](global/stage_g_consolidated_annotation/frozen_malignancy_axis_current_context_counts.tsv)展示原判定/层级、当前父类及QC。HSPC-like不能据此称正常HSC，新改的上皮/黑色素候选也不自动成为恶性细胞。

## 私有逐细胞文件与验证

完整逐细胞文件留在服务器，不公开细胞/患者行：

```text
${PROJECT_ROOT}/runs/20261106_v7_B_unresolved_review_v1/
  stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz
  stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.tsv.gz
```

两种格式均覆盖全部4,676,787行。pickle写出/读回比较，私有TSV全部分块核对细胞轴、原标签、当前标签与QC，所有文件SHA、工作目标不重叠、对照不覆盖和冻结V7 SHA均独立通过。[完成清单及私有文件大小/SHA](global/stage_g_consolidated_annotation/complete.json) · [独立合并审计](global/stage_g_consolidated_annotation/independent_consolidation_audit.json)。

## 下一步的具体优先级

先对B02/B04/B09/B10等父类混合分区按来源和正向身份组件纯化，再分别从counts建图；不能为了消除unresolved让一次整合图替细胞定父类。B05核对ASC正向组件、B程序与污染估计，B08逐细胞确认T装置及低信息。对已支持B父类的分区，再做naive/memory/GC候选、来源/分辨率/整合敏感性和逐细胞证据验收；IFN、应激、cycling和同种型维持独立状态。未重新整理的旧输入仍有早期来源偏倚，本轮不能作为干净TME重算结果或新CM证据。
