# 59个诊断群的完整工作注释与原因

以下为清理后res0.4计算分区；不是59个被接受的亚型。左侧名称是簇级工作描述；同群仍有逐细胞竞争证据者计入“需复查”，少数直接正向T/非HB红系路由可与簇描述不同。最终细分统一未接受。

[返回结果](README.md) · [完整数值TSV](all_59_cluster_review_with_evidence.tsv)

## Erythroid

| 群 | n | 工作描述 | 命名marker | 仍需身份复查 | 来源数 / 最大来源占比 | 原因 |
|---|---:|---|---|---:|---|---|
| E00 | 1,674 | T ／ globin-rich RNA | CD3D、CD3E、CD3G、TRBC2、LCK | 50 | 42 / 16.2% | CD3/TCR与DE一致，HB不证明红系。 |
| E01 | 1,417 | Ery candidate [review] | AHSP、SLC4A1、GYPA、ALAS2 | 1,098 | 31 / 27.7% | 非HB红系只在部分细胞支持；低深度且有T背景。 |
| E02 | 1,257 | Ery / epithelial [review] | AHSP、GYPA、EPCAM、KRT8、KRT18 | 1,045 | 31 / 25.4% | 红系及上皮组件混在同一分区，须逐细胞拆解。 |
| E03 | 888 | Ery candidate [review] | AHSP、GYPA、RHAG、ALAS2 | 717 | 32 / 29.8% | 有非HB红系子集，群内异质；当前无正向DE。 |
| E04 | 709 | Myeloid ／ C1Q/CD14 | CSF1R、C1QA、MS4A7、CD14 | 51 | 31 / 32.0% | 髓系组件与C1Q/CD14差异一致。 |
| E05 | 535 | Epithelial ／ KRT/EPCAM | EPCAM、KRT8、KRT18、KRT19 | 7 | 28 / 35.7% | KRT/EPCAM及ELF3/CLDN4差异支持上皮。 |
| E06 | 423 | B / ASC context [review] | CD74、MZB1、JCHAIN、IGKC | 371 | 30 / 21.0% | Ig高不是B/ASC身份；B、ASC与T成分须另核。 |
| E07 | 388 | Mural ／ RGS5/MYH11 | RGS5、PDGFRB、ACTA2、MYH11 | 9 | 26 / 23.2% | RGS5/PDGFRB/ACTA2/MYH11支持血管壁父类。 |
| E08 | 357 | Hepatocyte-like ／ ALB/HPX | ALB、APOA1、HPX、ASGR1 | 19 | 4 / 93.0% | ALB/APOA1/HPX/ASGR1支持肝样背景；不判良恶性。 |
| E09 | 354 | Epithelial ／ KRT19/MUC1 | KRT19、KRT18、EPCAM、MUC1 | 13 | 24 / 41.0% | KRT/MUC1及NAPSA/SFTA2差异支持上皮。 |
| E10 | 261 | Endothelial ／ CLDN5/VWF | CLDN5、VWF、PECAM1、KDR | 4 | 22 / 34.5% | CLDN5/VWF/PECAM1/KDR支持内皮。 |
| E11 | 201 | Myeloid ／ CSF3R/FPR1 | CSF3R、FPR1、TYROBP、S100A8 | 14 | 18 / 49.8% | CSF3R/FPR1及髓系程序；暂不强定粒细胞亚型。 |
| E12 | 189 | Erythroid ／ cycling | GYPA、SLC4A1、RHAG、AHSP、ALAS2 | 5 | 11 / 71.4% | 非HB结构/合成支持完整红系；cycling另轴。 |
| E13 | 142 | Epithelial ／ cilia context | EPCAM、KRT8、KRT18、FOXJ1、RSPH4A | 3 | 10 / 87.3% | KRT/EPCAM及纤毛差异；不接受组织亚型。 |
| E14 | 114 | Neural / endocrine [review] | PCSK1N、CPE、CHGB、SNAP25 | 114 | 12 / 40.4% | 神经/内分泌共享程序，组织父类仍待核。 |
| E15 | 112 | Melanocytic ／ PMEL/MLANA | PMEL、MLANA、MITF、TYR | 6 | 10 / 75.9% | PMEL/MLANA/MITF/TYR及差异支持黑色素谱系。 |
| E16 | 48 | Acinar-like [review] | PRSS1、PRSS2、CPA1、CTRC、CEL | 48 | 5 / 81.2% | 腺泡RNA一致，但小群且来源集中，完整细胞身份待核。 |
| E17 | 29 | Small source cluster [review] | TNFRSF11B、SYT1、SOSTDC1 | 29 | 1 / 100.0% | 仅29细胞单一来源；差异基因不构成稳定父类。 |

[marker完整统计](Erythroid/markers/measured_marker_expression.tsv) · [RNA DEG](Erythroid/markers/top15_RNA_DEG.tsv) · [原CellTypist/scATOMIC](Erythroid/markers/cluster_reference_top5.tsv) · [疾病背景聚合](Erythroid/cluster_disease_context.tsv)

## Platelet

| 群 | n | 工作描述 | 命名marker | 仍需身份复查 | 来源数 / 最大来源占比 | 原因 |
|---|---:|---|---|---:|---|---|
| P00 | 538 | Mural ／ RGS5/ACTA2 | RGS5、PDGFRB、ACTA2、MYH11、TAGLN | 1 | 27 / 30.5% | 血管壁程序一致；NRGN不特异。 |
| P01 | 448 | Mural ／ RGS5/ACTA2 | RGS5、PDGFRB、ACTA2、MYH11、TAGLN | 1 | 22 / 62.7% | 低深度仍有RGS5/ACTA2/MYH11组合。 |
| P02 | 411 | Mural ／ RGS5/ACTA2 | RGS5、PDGFRB、ACTA2、MYH11、TAGLN | 0 | 27 / 44.5% | ACTA2/TAGLN/MYH11差异与血管壁组合一致。 |
| P03 | 309 | Neural / platelet [review] | NCAM1、PCSK1N、PPBP、PF4、ITGA2B | 309 | 29 / 40.1% | 神经背景与血小板组件并存；不整群定血小板。 |
| P04 | 149 | Platelet / MK lineage | PF4、PPBP、ITGA2B、GP9、GP1BA、TUBB1 | 5 | 17 / 35.6% | 分泌＋结构组件强；先停血小板/巨核谱系。 |
| P05 | 123 | Epithelial / Ig [review] | EPCAM、KRT18、KRT8、IGKC | 122 | 8 / 85.4% | PPBP/Ig高而结构不足，伴上皮程序。 |
| P06 | 120 | Myeloid / platelet [review] | CD14、LST1、TYROBP、PF4、GP9 | 120 | 23 / 15.0% | 髓系与血小板组件并存，可能附着/环境RNA/混合。 |
| P07 | 77 | Neural / endocrine [review] | PCSK1N、CHGB、CPE、SNAP25、NEFL | 77 | 7 / 45.5% | 神经/内分泌组合，血小板结构不足。 |
| P08 | 61 | T / platelet [review] | CD3D、CD3E、TRBC2、PF4、GP9 | 58 | 17 / 57.4% | CD3/TCR与血小板结构共存，保留混合复查。 |
| P09 | 59 | Mural ／ RGS5/ACTA2 | RGS5、PDGFRB、ACTA2、MYH11、TAGLN | 0 | 15 / 62.7% | RGS5/PDGFRB/ACTA2/MYH11支持血管壁。 |
| P10 | 23 | Mural ／ RGS5/ACTA2 | RGS5、PDGFRB、ACTA2、MYH11、TAGLN | 0 | 2 / 73.9% | 血管壁程序强，伴腺泡差异，小群来源集中。 |

[marker完整统计](Platelet/markers/measured_marker_expression.tsv) · [RNA DEG](Platelet/markers/top15_RNA_DEG.tsv) · [原CellTypist/scATOMIC](Platelet/markers/cluster_reference_top5.tsv) · [疾病背景聚合](Platelet/cluster_disease_context.tsv)

## HSPC

| 群 | n | 工作描述 | 命名marker | 仍需身份复查 | 来源数 / 最大来源占比 | 原因 |
|---|---:|---|---|---:|---|---|
| H00 | 3,171 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 7 | 35 / 71.4% | primitive与髓系背景；低深度且AML集中，不能定正常HSC。 |
| H01 | 2,881 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 8 | 41 / 66.9% | SPINK2/CD34/FLT3背景；当前无正向DE，不强分亚型。 |
| H02 | 2,510 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 39 | 25 / 71.6% | primitive＋髓系及S期RNA；分化与增殖分开。 |
| H03 | 1,288 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 5 | 18 / 49.8% | MEIS1/GATA2与primitive背景；共享marker且当前无正向DE。 |
| H04 | 1,101 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 13 | 15 / 76.1% | 髓系priming与MKI67/TOP2A增殖，不能据此判doublet。 |
| H05 | 748 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 32 | 8 / 93.7% | primitive背景，单来源占94%；不命名正常MPP。 |
| H06 | 598 | Endothelial-like [review] | CLDN5、CD34、FLT1、EMCN、RAMP2 | 596 | 38 / 31.6% | 内皮样组件与共享CD34/GATA2；仍有造血候选。 |
| H07 | 526 | Mural ／ RGS5/MYH11 | RGS5、PDGFRB、ACTA2、MYH11 | 2 | 27 / 30.2% | 血管壁程序一致，非HSPC子型。 |
| H08 | 277 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 2 | 6 / 89.2% | primitive背景，来源集中；不作功能性HSC判定。 |
| H09 | 271 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 3 | 8 / 97.0% | primitive与血液病语境，不作离散祖细胞亚型。 |
| H10 | 173 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 0 | 2 / 80.3% | primitive/髓系组合，主要两来源；不强定GMP。 |
| H11 | 160 | Neural / endocrine [review] | PCSK1N、CHGB、CPE、SNAP25 | 160 | 15 / 55.6% | PCSK1N/CPE/CHGB等神经/内分泌程序。 |
| H12 | 126 | HSPC-like ／ myeloid axis | SPINK2、CD34、FLT3、GATA2、MPO | 1 | 3 / 95.2% | 髓系/primitive与S期背景，来源集中。 |
| H13 | 120 | Immature myeloid-like | LYZ、MPO、ELANE、CTSG、S100A8 | 1 | 7 / 78.3% | MPO/ELANE/CTSG/LYZ未成熟髓系背景，非正常HSPC亚型。 |
| H14 | 84 | Stromal / KRT [review] | COL1A1、COL3A1、KRT5、KRT14、MEIS1 | 84 | 4 / 75.0% | 基质/KRT/MEIS1共享并有IFN，转化/混合未定。 |
| H15 | 38 | Myeloid ／ CTSK/ACP5 | TYROBP、SPI1、CSF1R、CTSK、ACP5 | 2 | 4 / 84.2% | 髓系＋CTSK/ACP5程序，未证明破骨功能。 |

[marker完整统计](HSPC/markers/measured_marker_expression.tsv) · [RNA DEG](HSPC/markers/top15_RNA_DEG.tsv) · [原CellTypist/scATOMIC](HSPC/markers/cluster_reference_top5.tsv) · [疾病背景聚合](HSPC/cluster_disease_context.tsv)

## ILC

| 群 | n | 工作描述 | 命名marker | 仍需身份复查 | 来源数 / 最大来源占比 | 原因 |
|---|---:|---|---|---:|---|---|
| I00 | 12,867 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 8,898 | 51 / 40.5% | IL7R/KLRB1共享且T apparatus丰富；有正向T者转出。 |
| I01 | 9,283 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 8,955 | 46 / 63.0% | 低深度IL7R/RORA，共享程序，未检出TCR不证明ILC。 |
| I02 | 8,353 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 7,623 | 49 / 21.9% | IL7R/KLRB1与T背景，当前无正向DE，暂留边界。 |
| I03 | 4,758 | KIT/XCL lymphoid [review] | IL7R、KIT、KLRB1、XCL1、XCL2、TRBC1 | 4,375 | 50 / 22.7% | KIT/XCL、ILC候选与T/NK证据并存，不能整群分ILC3。 |
| I04 | 3,352 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 2,794 | 42 / 55.7% | 淋巴/T背景，IL7R/GATA3不能单独定ILC2。 |
| I05 | 2,509 | Melanocytic / mural [review] | PMEL、MLANA、MITF、ACTA2、TAGLN | 2,439 | 46 / 31.1% | 黑色素/血管壁程序；ILC组件阳性也须先有造血父类。 |
| I06 | 1,881 | KIT/RORC low-depth [review] | KIT、RORC、MLANA、PMEL | 1,867 | 46 / 22.4% | 低深度KIT/RORC与黑色素候选，不能仅凭KIT/RORC定ILC3。 |
| I07 | 565 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 456 | 45 / 18.1% | IFN另轴；父身份仍需T/ILC边界复核。 |
| I08 | 453 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 363 | 19 / 68.7% | 有T/淋巴身份并伴肝RNA；不据ALB直接转肝细胞。 |
| I09 | 414 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 360 | 43 / 18.8% | cycling另轴；不把增殖群叫新ILC亚型。 |
| I10 | 341 | Melanocytic ／ PMEL/MLANA | PMEL、MLANA、MITF、TYR | 22 | 13 / 87.1% | PMEL/MLANA/MITF/TYR支持黑色素谱系。 |
| I11 | 292 | Myeloid / lymphoid [review] | LYZ、TYROBP、FCER1G、MS4A7、IL7R | 272 | 40 / 13.7% | 髓系与IL7R淋巴程序混合，按细胞复核。 |
| I12 | 155 | Melanocytic ／ MITF/KIT | MLANA、MITF、TYR、PMEL | 0 | 4 / 98.1% | MITF/MLANA/KIT及黑色素背景，来源高度集中。 |
| I13 | 67 | IL7R lymphoid [review] | IL7R、KLRB1、RORA、CD3E、TRBC2 | 29 | 1 / 100.0% | 单来源、原始库量高且技术不可一概称UMI；T/NK边界复核。 |

[marker完整统计](ILC/markers/measured_marker_expression.tsv) · [RNA DEG](ILC/markers/top15_RNA_DEG.tsv) · [原CellTypist/scATOMIC](ILC/markers/cluster_reference_top5.tsv) · [疾病背景聚合](ILC/cluster_disease_context.tsv)
