# 原有已命名标签保留审计

不再要求原有亚型重新通过v3的大类组件门槛后才准保留。30个不含占位词的原标签逐一恢复；这30个包含亚型、状态和恶性语境分类，不等于30种新验收亚型。

| 原标签 | QC保留 | 原名保留 | 此前重聚类例外 | 其中原v3未纳入、此次恢复原名 |
|---|---:|---:|---:|---:|
| Bcell_Naive_IGHD | 105,825 | 105,825 | 0 | 783 |
| Bcell_Naive_TCL1A | 128,798 | 128,798 | 0 | 2,715 |
| Dendritic_pDC_LILRA4 | 2,072 | 2,072 | 0 | 464 |
| Endothelial_Blood_PECAM1 | 21,745 | 21,745 | 0 | 2,036 |
| Endothelial_Venous_ACKR1_ACKR1 | 16,237 | 16,237 | 0 | 550 |
| Epithelial_malignancy_discordant | 229,472 | 229,472 | 0 | 12,382 |
| Epithelial_malignancy_not_evaluable | 163,426 | 163,426 | 0 | 2,779 |
| Erythroid_Cycling_Erythroid_TOP2A | 2,348 | 326 | 2,022 | 3 |
| Fibroblast_Cycling_Fibroblast_TOP2A | 39,720 | 39,720 | 0 | 6,214 |
| Fibroblast_ECM_CAF_COL1A1 | 48,829 | 48,829 | 0 | 6,189 |
| Fibroblast_iCAF_CXCL12 | 35,918 | 35,918 | 0 | 2,245 |
| Fibroblast_iCAF_IL6 | 26,748 | 26,748 | 0 | 1,023 |
| Fibroblast_myCAF_ACTA2 | 89,122 | 89,122 | 0 | 21,393 |
| Granulocyte_Neutrophil_mature_FCGR3B | 13,977 | 13,977 | 0 | 447 |
| Macrophage_Cycling_Macrophage_MKI67 | 16,129 | 16,129 | 0 | 2,120 |
| Macrophage_Inflammatory_IL1B_IL1B | 140,985 | 140,985 | 0 | 6,723 |
| Macrophage_SPP1_TREM2_SPP1 | 70,223 | 70,223 | 0 | 3,399 |
| Malignant_Epithelial_high_specificity | 185,305 | 185,305 | 0 | 23,368 |
| Mast_Activated_IL1RL1_IL1RL1 | 10,701 | 10,701 | 0 | 192 |
| Mast_Activated_IL1RL1_LTC4S | 10,741 | 10,741 | 0 | 559 |
| Monocyte_Classical_FCN1_FCN1 | 16,626 | 16,626 | 0 | 576 |
| Monocyte_Inflammatory_IL1B_IL1B | 40,731 | 40,731 | 0 | 1,683 |
| Monocyte_Nonclassical_FCGR3A_FCGR3A | 11,681 | 11,681 | 0 | 7,650 |
| NK_FCGR3A_cytotoxic_FCGR3A | 60,881 | 60,881 | 0 | 13,424 |
| NK_XCL1_CD56bright_XCL1 | 37,526 | 37,526 | 0 | 8,278 |
| Neural_Neural_progenitor_SOX2 | 62,127 | 62,127 | 0 | 3,930 |
| Nonmalignant_Epithelial_CNV_negative | 199,186 | 199,186 | 0 | 4,510 |
| Tcell_CD4_naive_IL7R | 519,131 | 519,131 | 0 | 1,466 |
| Tcell_CD8_cytotoxic_CCL5 | 496,995 | 496,995 | 0 | 2,003 |
| Tcell_Cycling_T_MKI67 | 41,868 | 41,868 | 0 | 4,120 |

唯一发生上述例外的原已命名标签及依据：
- `Erythroid_Cycling_Erythroid_TOP2A`：2,022个细胞已在此前rare-four counts清理重聚类中被分流为其他父类或父身份未定。保留其真实重聚类结果；仍支持原父身份的326个细胞保留原名。此前证据见[四类清理结果](../v7_other_unresolved_review_v1/README.md)。这不是用v3粗粒度门槛删掉旧亚型。

[原亚型与v3粗粒度父类差异表](stage_h_selective_subtype_restore_v4/retained_annotations_parent_disagreement_for_review.tsv)及[所有继承标签的大类诊断语境](stage_h_selective_subtype_restore_v4/retained_annotations_broad_screen_context.tsv)保留审查依据。v3诊断父类与原标签存在差异的其他记录，保留原名并在独立字段`broad_parent_difference_needs_review`标示，不把粗粒度组件命中直接当作已经推翻原亚型。候选细分不自动升级成确定亚型；父身份仍混合者继续保留在全轴审计对象中。
