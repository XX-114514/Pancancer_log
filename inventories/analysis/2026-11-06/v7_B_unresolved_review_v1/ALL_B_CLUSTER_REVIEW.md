# B未定群全部16个诊断分区

这是本轮157,888个注释目标的全表；8,565个参考仅帮助建模，单列人数。群级父类/候选背景与逐细胞竞争复查后的工作记录分开。大类工作不等于最终细分通过，状态不是逐细胞功能证明。

[报告](README.md) · [机器可读完整决策](global/review_decisions/B_cluster_decisions.tsv) · [逐细胞结果的聚合](global/stage_f_B_working_annotations/B_old_to_working_annotation_counts.tsv)

| 群 | 目标 | 参考 | 大类工作目标 | 身份待复查 | 目标来源 | 工作父类／候选 | 命名marker | 保留未定的原因或接受边界 |
|---|---:|---:|---:|---:|---:|---|---|---|
| B00 | 41,402 | 2,377 | 41,166 | 236 | 51 | Bcell / memory_like_candidate | MS4A1;CD79A;BANK1;CD27;TNFRSF13B;AIM2 | B身份组合明确；CD27/TNFRSF13B/AIM2支持memory-like候选；FOS/JUN和应激单列，非全群细分验收。 |
| B01 | 23,127 | 1,259 | 23,029 | 98 | 49 | Bcell / naive_like_candidate | MS4A1;CD79A;BANK1;TCL1A;IGHD;IGHM;FCER2 | TCL1A/IGHD/IGHM/FCER2组合支持naive-like；以B大类工作记录，未给全群细胞强配最终naive标签。 |
| B02 | 21,812 | 15 | 235 | 21,577 | 50 | Parent_unresolved / not_assigned | CST3;MGST1;EPCAM;KRT8;CLDN5;PMEL | 父身份混合且信息量偏低；不能由CST3或HLA留作B或判DC，部分正向B/T可单独路由。 |
| B03 | 15,066 | 257 | 15,018 | 48 | 51 | Bcell / not_assigned | MS4A1;CD79A;BANK1 | 全RNA仍有一致B身份；当前DE无过阈值正向基因并非无身份，较低信息量支持先停B大类。 |
| B04 | 14,761 | 31 | 221 | 14,540 | 51 | Nonhematopoietic_parent_review / not_assigned | AKAP12;IGFBP7;CLU;EPCAM;CLDN5;MLANA | 上皮/内皮/血管壁/黑色素等不同父类混合；UMI中位4971，不能全归因低深度，需按父类和来源拆开重建。 |
| B05 | 5,529 | 3,623 | 565 | 4,964 | 48 | B_ASC_boundary / ASC_like_component_candidate | MZB1;DERL3;JCHAIN;PRDM1;XBP1;MS4A1 | 有明确ASC分泌成分但仍有B身份；仅目标marker绘图，不让3623个对照抬高证据，不强称成熟度/过渡或doublet。 |
| B06 | 8,486 | 420 | 7,978 | 508 | 51 | Myeloid / cDC2_macrophage_boundary_candidate | LYZ;TYROBP;LST1;CSF1R;CD163;CLEC10A;CD1C | 髓系大类成立；cDC2与巨噬程序并存，不把整个分区强称一种DC或巨噬亚型。 |
| B07 | 8,245 | 19 | 8,048 | 197 | 51 | Myeloid / LAMP3_DC_like_candidate | LAMP3;FSCN1;CCR7;TYROBP;IL4I1 | LAMP3/FSCN1/CCR7支持DC样迁移状态候选；保留髓系大类与DC候选分列，尚非全群最终DC细分。 |
| B08 | 6,027 | 26 | 473 | 5,554 | 51 | T_candidate_low_information / not_assigned | CD3G;CD3D;CD3E;TRAC;TRBC2;BCL11B | T方向富集但仅393个目标满足保守T组件规则，不能给6027个全部强改T；正向细胞分别路由，其余待复查。 |
| B09 | 4,630 | 7 | 109 | 4,521 | 50 | Parent_unresolved / not_assigned | MS4A1;KRT8;CLDN5;TYROBP | 无过阈值正向DE且多父类弱信号，仍属父身份未定，不能当未命名B亚型。 |
| B10 | 3,804 | 207 | 209 | 3,595 | 48 | Parent_unresolved / not_assigned | MKI67;TOP2A;UBE2C;CENPF | 增殖程序主导但父类混合；cycling不是B身份，不能命名为增殖B并据此接受细分。 |
| B11 | 3,115 | 128 | 3,068 | 47 | 49 | Bcell / memory_like_context_not_final | MS4A1;CD79A;BANK1;IFIT1;IFIT3;MX1;ISG15 | B身份明确、IFN程序突出；按B大类加IFN群级背景记录，不造新的功能性B亚型。 |
| B12 | 896 | 160 | 888 | 8 | 42 | Bcell / GC_like_candidate | MS4A1;CD79A;RGS13;BCL6;AICDA;MEF2B | RGS13及BCL6/AICDA/MEF2B组合支持GC-like候选，仍有naive/memory背景，保留候选而非全群终验。 |
| B13 | 580 | 2 | 20 | 560 | 33 | Parent_unresolved / not_assigned | PDZK1;APOM;KHK;APOE;SPP1 | 组织/髓系RNA混合，B支持弱；不能用ALB或APOE单独定肝细胞/巨噬，需来源内拆分。 |
| B14 | 251 | 0 | 251 | 0 | 5 | Endothelial / not_assigned | CLDN5;PECAM1;FLT1;TIE1;PROCR;MMRN1 | 内皮组合明确，可纠正B父标签；92.8%来自同一来源，血管细分及泛癌推广未验收。 |
| B15 | 157 | 34 | 21 | 136 | 4 | Parent_unresolved / not_assigned | ALB;APOA1;ASGR1;APOE | 肝样RNA明显但小群来源偏斜、污染估计偏高且有B/髓系成分；不强给全群肝细胞或恶性标签。 |

全部群级状态、英文理由、最终细分FALSE及工作标签在TSV完整保留。开放分区中的正向B/T个体可单独路由，因此某群“工作父类”不必等于其全部目标的逐细胞最终工作父类。无阳性组件不等于阴性，剩余记录保留复查。
