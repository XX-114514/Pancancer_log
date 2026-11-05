# 修订后的全部基因 panel

共 112 个证据模块，覆盖旧 162 项。此资源是文献引导的候选注册表，不是已经校准的自动细胞分类器。

“核心组件”中的分组用于并列检查不同证据，不表示每个细胞必须检出所有基因；不会用一个共享基因替代整个身份。“降低基因”只在可测且有适当对照时解释，掉失/未测不算阴性。

参考等级：author_exact_signature=作者完整功能签名；curated_literature_informed=依据论文/代码整理的项目组合；curated_reference_gap=尚缺逐项核验的适用原始参考，禁止自动终判。

## parent

| 模块 | 父类/适用 | 全部正向基因 | 期望降低（需对照） | 命名与边界 | 依据 |
|---|---|---|---|---|---|
| T | 独立入口 | CD3D, CD3E, CD3G, TRAC, TRBC1, TRBC2, TRDC, TRGC1, TRGC2, LCK, CD247, CD2 | 无自动负向 | CD3与TCR组件分开；αβ/γδ均可；低覆盖不作阴性。 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal); [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| NK | 独立入口 | NCR1, KLRD1, KLRF1, NKG7, GNLY, PRF1, FCGR3A, NCAM1, XCL1, XCL2 | 无自动负向 | NK受体与效应联合；TCR/CD3另核，绝不由细胞毒性直接定NK。 | [YangNK2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6718415/); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| ILC | 独立入口 | IL7R, IL1R1, GATA3, RORC, RORA, TBX21, KIT, KLRB1, PTGDR2, IL1RL1 | 无自动负向 | 仅ILC候选入口；必须检查Lin/TCR与NK/Mast边界，不能单靠IL7R。 | [Bjorklund2016](https://pubmed.ncbi.nlm.nih.gov/26878113/) |
| B | 独立入口 | CD79A, CD79B, MS4A1, CD19, CD22, CD37, BANK1, CD74, CD40 | 无自动负向 | 成熟B组合；CD79A/CD37与ASC共享。 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal); [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| ASC | 独立入口 | MZB1, DERL3, PRDM1, XBP1, SDC1, TNFRSF17, JCHAIN, CD38 | 无自动负向 | 抗体分泌细胞入口；不由单个XBP1或JCHAIN判浆细胞。 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal); [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| Myeloid | 独立入口 | LYZ, LST1, SPI1, TYROBP, FCER1G, CSF1R, CTSD | 无自动负向 | 广义髓系入口而非单核亚型；部分DC/粒细胞须并行召回。 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal) |
| DC | 独立入口 | FLT3, HLA-DRA, CD74, CLEC9A, XCR1, CD1C, FCER1A, LILRA4, CLEC4C | 无自动负向 | DC父身份由下列分支并行核查后汇总；该粗模块不是所有DC必需的硬AND。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Mast | 独立入口 | TPSAB1, TPSB2, CPA3, KIT, MS4A2, HDC, HPGDS, LTC4S | 无自动负向 | 先确认Mast；无稳定亚型时停大类，HDC/LTC4S不证明激活。 | [Tauber2023](https://pubmed.ncbi.nlm.nih.gov/37462672/); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Basophil | 独立入口 | IL3RA, FCER1A, CLC, MS4A2, HDC, ENPP3 | 无自动负向 | basophil邻近谱系候选，原Mast/Basophil_like移出；须完整对照。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Neutrophil | 独立入口 | CSF3R, FCGR3B, S100A8, S100A9, FPR1, CXCR2, CEACAM8 | 无自动负向 | 成熟与未成熟分支并查；LYZ不是髓系污染。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Eosinophil | 独立入口 | EPX, PRG2, CLC, SIGLEC8, IL5RA | 无自动负向 | 嗜酸粒细胞与basophil/Mast组合对照；CLC非独占。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Epithelial | 独立入口 | KRT8, KRT18, KRT19, KRT7, EPCAM, KRT5, KRT14, MUC1 | 无自动负向 | 允许组织/EMT例外；恶性另轴，不含melanocyte。 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal); [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Endothelial | 独立入口 | PECAM1, CDH5, VWF, KDR, FLT1, CLDN5, RAMP2, EMCN, ESAM | 无自动负向 | 血管/淋巴分支并查，不要求所有亚型表达VWF。 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Fibroblast | 独立入口 | DCN, LUM, DPT, COL1A1, COL1A2, COL3A1, PDGFRA | 无自动负向 | 先证实fibroblast；肿瘤来源才能考虑CAF名称。 | [Buechler2021](https://www.nature.com/articles/s41586-021-03549-5); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Pericyte | 独立入口 | RGS5, CSPG4, PDGFRB, MCAM, NOTCH3, COX4I2 | 无自动负向 | mural谱系，不是fibroblast必然子型。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf); [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Smooth_muscle | 独立入口 | MYH11, CNN1, DES, TAGLN, ACTA2 | 无自动负向 | 收缩程序不能独立证明SMC；MYH11及完整背景联合。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf); [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Mesothelial | 独立入口 | MSLN, UPK3B, CALB2, KRT8, KRT18, KRT19 | 无自动负向 | 独立邻近谱系，来源及腔面背景需要核查。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Erythroid | 独立入口 | GYPA, AHSP, SLC4A1, GATA1, KLF1, HBB, HBA1, HBA2, TFRC, ALAS2 | 无自动负向 | 不能只用HB；分化阶段可选，不强制亚型。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| Megakaryocytic | 独立入口 | GP9, ITGA2B, GP1BA, PF4, PPBP, TUBB1, RGS18, NRGN | 无自动负向 | 巨核/血小板共享入口；附着、环境RNA、核/复杂度另查。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| HSPC | 独立入口 | SPINK2, HLF, MEIS1, CD34, PROM1, GATA2, HOPX | 无自动负向 | 造血祖细胞候选；AML不能自动叫正常HSC，成熟谱系排查。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| Neural_glial | 独立入口 | SOX10, PLP1, MPZ, SLC1A3, MBP, MOG, RBFOX3, SNAP25, S100B | 无自动负向 | 仅广义神经/胶质入口，以Schwann/少突/神经元等分支OR复核；恶性独立。 | [Neftel2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/) |
| Melanocytic | 独立入口 | PMEL, MLANA, TYR, DCT, SOX10, MITF | 无自动负向 | 从Neural/Epithelial移出；MITF-low/AXL-high肿瘤可弱表达，低检出不排除来源；恶性另判。 | [Tirosh2016](https://shaleklab.com/wp-content/uploads/2024/03/Tirosh-et-al-Science-2016.pdf) |

## branch

| 模块 | 父类/适用 | 全部正向基因 | 期望降低（需对照） | 命名与边界 | 依据 |
|---|---|---|---|---|---|
| Monocyte | Myeloid | FCN1, CD14, VCAN, LYZ, LST1, LILRB1, S100A8, S100A9, CTSD | 无自动负向 | 单核/巨噬连续谱；必须与粒细胞/DC比较。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf); [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal) |
| Macrophage | Myeloid | C1QA, C1QB, C1QC, APOE, CSF1R, CD68, MSR1, MARCO, APOC1, FOLR2 | 无自动负向 | 允许炎症/单核来源群低C1Q；不得以此筛掉整个亚群。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf); [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal) |
| T_CD4 | T | CD4, IL7R, LTB | 无自动负向 | CD4 RNA掉失时保持T分支待核；IL7R不代替CD4。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_CD8 | T | CD8B, CD8A | 无自动负向 | CD8A/B分开检查；CD8A也见NK；不要求每个细胞同时检出。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_Treg | T | FOXP3, IL2RA, IKZF2, IKZF4, CTLA4, TIGIT, TNFRSF18, TNFRSF4 | 无自动负向 | FOXP3与配套调节组件；先T身份，不以TIGIT/CTLA4或模型多数票独定Treg。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_Tfh_like | T_CD4 | BCL6, CXCR5, CXCL13, ICOS, PDCD1, TOX2, SH2D1A | 无自动负向 | 保留Tfh-like；与Treg/Tph/激活程序比较，不凭CXCL13命名Tfh。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_Tph_like | T_CD4 | CXCL13, PDCD1, ICOS, MAF, IL21 | CXCR5 | 需CXCR5相对低及Tfh对照；RA功能结论不直接移植TME，仅Tph-like候选。 | [Rao2017](https://www.nature.com/articles/nature20810); [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_Th1_like | T_CD4 | TBX21, CXCR3, IFNG, IL12RB2, STAT4 | 无自动负向 | 方向程序，不由IFNG倒推CD4。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_Th17_like | T_CD4 | RORC, CCR6, IL17A, IL17F, KLRB1 | 无自动负向 | IL17低检出保留不确定；KLRB1不等于Th17。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_gd | T | TRDC, TRGC1, TRGC2, KLRB1, KLRG1 | 无自动负向 | γδ TCR组合；VDJ可用时核查。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| T_MAIT | T | SLC4A10, TRAV1-2, KLRB1, ZBTB16, CXCR6 | 无自动负向 | MAIT-like候选；TCR/蛋白证据增强，不能凭KLRB1判定。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| ILC1 | ILC | TBX21, IFNG, IL12RB2, CXCR3 | 无自动负向 | ILC1与NK必须边界复核；无受体/谱系证据时停ILC-like。 | [Bjorklund2016](https://pubmed.ncbi.nlm.nih.gov/26878113/) |
| ILC2 | ILC | GATA3, PTGDR2, IL1RL1, RORA, IL5, IL13 | 无自动负向 | 类型2程序配套受体；IL7R不能作为唯一类型证据。 | [Bjorklund2016](https://pubmed.ncbi.nlm.nih.gov/26878113/) |
| ILC3 | ILC | RORC, IL23R, CCR6, IL22, KIT, KLRB1 | 无自动负向 | RORC/IL23R等联合；KIT/KLRB1仅辅助。 | [Bjorklund2016](https://pubmed.ncbi.nlm.nih.gov/26878113/) |
| B_naive_like | B | TCL1A, IGHD, IGHM, IL4R | 无自动负向 | naïve-like候选，MS4A1不再是亚型anchor。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038); [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal) |
| B_memory_like | B | CD27, TNFRSF13B, AIM2, BANK1, IGHG1, IGHA1 | 无自动负向 | CD27阴性记忆也可能存在；类别转换另轴，需整体参考。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| B_GC | B | BCL6, AICDA, RGS13, MEF2B, LMO2 | 无自动负向 | GC组合；不把cycling当GC身份。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| B_atypical_like | B | TBX21, ITGAX, FCRL5, ZEB2 | 无自动负向 | atypical-like/TBX21-FCRL5-associated；不推断年龄或TAAB功能。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| Mono_classical | Monocyte | FCN1, CD14, S100A8, S100A9, VCAN, LYZ | 无自动负向 | 组合支持，检查髓系连续谱。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Mono_FCGR3A | Monocyte | FCGR3A, MS4A7, CX3CR1, LST1 | 无自动负向 | FCGR3A-associated单核；非经典终点须比较CD14程序。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Mono_intermediate | Monocyte | CD14, FCGR3A, HLA-DRA, MS4A7 | 无自动负向 | 要求双程序的群内/同细胞证据，不由HLA-DRA独定。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| DC_cDC1 | DC | CLEC9A, XCR1, CADM1, BATF3, IRF8 | 无自动负向 | 作者正文组合，分组件为本项目实现。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) |
| DC_cDC2 | DC | CD1C, FCER1A, CLEC10A, CD1E | 无自动负向 | CD1C需配套证据；DC3另查。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) |
| DC_DC3 | DC | CD1C, CD14, CD163, LILRB1, S100A8, S100A9 | 无自动负向 | DC3-like需CD1C+单核相关组合和鉴别；不能CD1C单项。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) |
| DC_pDC | DC | LILRA4, CLEC4C, TCF4, GZMB, JCHAIN, IRF7, MZB1 | 无自动负向 | LILRA4/CLEC4C与TCF4优先；GZMB/JCHAIN是共享辅助。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) |
| DC_AS | DC | AXL, SIGLEC6, CD2, CD5 | 无自动负向 | 稀有AS-DC候选，需真实细胞与批次稳定性；不强制拆群。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) |
| DC_LC_like | DC | CD1A, CD207, EPCAM | 无自动负向 | LC-like候选；EPCAM在此不是自动污染排除。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/); [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Endo_lymphatic | Endothelial | PROX1, FLT4, LYVE1, CCL21, PDPN, MMRN1 | 无自动负向 | LYVE1/PDPN不独占；需内皮与淋巴组合。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Endo_arterial | Endothelial | EFNB2, GJA5, SOX17, HEY1, DKK2 | 无自动负向 | 来源/器官层级验证。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Endo_venous | Endothelial | ACKR1, NR2F2, PLVAP | 无自动负向 | 仅ACKR1-associated venous候选；PLVAP不能替代命名marker。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Endo_capillary | Endothelial | CA4, RGCC, CD36, CLDN5, EMCN | 无自动负向 | 器官特异，不能覆盖全部毛细血管类型。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Schwann | Neural_glial | MPZ, SOX10, PLP1, S100B, S100A1 | 无自动负向 | 神经鞘候选；SOX10单项不够。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Oligodendrocyte | Neural_glial | MBP, MOG, PLP1, OLIG1, OLIG2 | 无自动负向 | 正常/恶性另判；脑参考需细化。 【参考缺口】 | [Neftel2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/) |
| Neuron | Neural_glial | RBFOX3, SNAP25, SYT1, TUBB3 | 无自动负向 | 神经元候选，不由TUBB3单项判定；需组织适配脑/外周神经参考。 【参考缺口】 | [Neftel2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/) |
| Epi_basal | Epithelial | KRT5, KRT14, TP63, KRT17, SFN | 无自动负向 | 器官/来源与良恶性单独核查；RNA不等于正常细胞。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Epi_AT2 | Epithelial | SFTPC, SFTPB, ABCA3, SFTPA1, SFTPA2 | 无自动负向 | 器官/来源与良恶性单独核查；RNA不等于正常细胞。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Epi_AT1 | Epithelial | AGER, CAV1, EMP2, CAV2, PDPN | 无自动负向 | 器官/来源与良恶性单独核查；RNA不等于正常细胞。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Epi_club | Epithelial | SCGB1A1, SCGB3A1, KRT8, KRT19 | 无自动负向 | 器官/来源与良恶性单独核查；RNA不等于正常细胞。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Epi_ciliated | Epithelial | FOXJ1, TPPP3, PIFO, CAPS, TUBB4B | 无自动负向 | 器官/来源与良恶性单独核查；RNA不等于正常细胞。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Epi_hepatocyte_like | Epithelial | ALB, ASGR1, HNF4A, HP, APOA1, TTR, APOA2, CYP3A4 | 无自动负向 | ALB与ASGR1/HNF4A肝细胞方向；环境ALB与去分化另核，不能由EPCAM缺失排除。 | [Aizarani2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6687507/) |
| Epi_cholangiocyte_like | Epithelial | KRT19, CFTR, KRT7, EPCAM, KRT8, KRT18, MMP7 | 无自动负向 | KRT19与CFTR成熟胆管方向；KRT7不再独立定义，肿瘤和前体细胞可偏离。 | [Aizarani2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6687507/) |
| Epi_renal_tubular_like | Epithelial | SLC34A1, LRP2, CUBN, HAVCR1, EPCAM | 无自动负向 | 旧组织panel保守改写；缺已逐条核验的器官原始参考，默认reference_gap，不作自动终判。 【参考缺口】 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal) |
| Epi_prostate_luminal_like | Epithelial | KLK3, KLK2, ACP3, MSMB, NKX3-1, EPCAM, KRT8, KRT18 | 无自动负向 | 前列腺组织内KLK3/KLK2与ACP3/MSMB方向；NKX3-1为项目辅助项，不能由此判良恶性或覆盖所有肿瘤状态。 | [Henry2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6411034/) |

## state

| 模块 | 父类/适用 | 全部正向基因 | 期望降低（需对照） | 命名与边界 | 依据 |
|---|---|---|---|---|---|
| Cytotoxic | T, NK, ILC | PRF1, GZMB, NKG7, GNLY, CCL5, FGFBP2, GZMA | 无自动负向 | 细胞毒程序，不能确定CD8/NK。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| GZMK_program | T, NK | GZMK, GZMA, CCL5, IL7R, CXCR3 | 无自动负向 | GZMK-associated，不命名CD8。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| Naive_like | T | TCF7, LEF1, CCR7, SELL, MAL, LTB | 无自动负向 | naïve/central-memory边界需另判，CD4/CD8独立。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| Exhaustion_like | T | PDCD1, HAVCR2, TOX, LAG3, TIGIT, CXCL13, ENTPD1 | 无自动负向 | exhausted-like程序；Treg/激活对照，功能耗竭未证明。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| Resident_like | T, NK | ITGAE, ZNF683, CXCR6, CD69, RGS1 | 无自动负向 | 驻留样程序；CD69单项不证明驻留。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| IFN_response | * | ISG15, IFIT1, IFIT3, MX1, OAS1, IFIT2 | 无自动负向 | 跨谱系状态，不作身份否决。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| Cycling | * | MKI67, TOP2A, UBE2C, TYMS, STMN1 | 无自动负向 | 增殖状态；不把platelet升级为megakaryocyte。 | [Lodi2025](https://github.com/lambrechtslab/Pan-cancer_Lodietal) |
| Inflammatory | Myeloid, Monocyte, Macrophage | IL1B, CXCL8, TNF, CCL3, CCL4, CXCL2 | 无自动负向 | 炎症程序，可与C1Q/SPP1共存。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| MHCII | * | HLA-DRA, HLA-DRB1, CD74, CIITA, HLA-DPA1, HLA-DPB1 | 无自动负向 | 呈递相关表达，不证明呈递功能；apCAF须先fibroblast。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) |
| DC_maturation | DC | LAMP3, CCR7, FSCN1, CCL19, CD274 | 无自动负向 | 成熟迁移样状态；来源cDC身份可未知，不强行溯源。 | [Ma2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/) |
| Mac_C1Q | Macrophage | C1QA, C1QB, C1QC, APOE, FOLR2 | 无自动负向 | C1Q/APOE-associated，不等于M2。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Mac_SPP1 | Macrophage | SPP1, TREM2, GPNMB, APOC1, LGALS3, FABP5 | 无自动负向 | SPP1程序，不证明抑制功能。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Mac_MRC1 | Macrophage | MRC1, CD163, FOLR2, SELENOP | 无自动负向 | marker-associated，不直接M2。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| NK_FCGR3A | NK | FCGR3A, FGFBP2, PRF1, GNLY, GZMB | 无自动负向 | FCGR3A相关，不作为NK父身份唯一证据。 | [YangNK2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6718415/) |
| NK_XCL | NK | XCL1, XCL2, NCAM1, KLRC1, GZMK, IL7R | 无自动负向 | RNA称XCL-associated；不声称测得CD56蛋白bright。 | [YangNK2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6718415/) |
| Mast_IL1RL1 | Mast | IL1RL1, AREG, LTC4S, HDC | 无自动负向 | 只描述IL1RL1相关表达，取消activated断言。 | [Tauber2023](https://pubmed.ncbi.nlm.nih.gov/37462672/) |
| Mast_CMA1 | Mast | CMA1, CPA3, TPSAB1, TPSB2 | 无自动负向 | CMA1与CPA3不同；不按RNA直接等同蛋白MCTC。 | [Tauber2023](https://pubmed.ncbi.nlm.nih.gov/37462672/) |
| Mast_CPA3 | Mast | CPA3, HPGDS, HDC | 无自动负向 | 描述marker，不强制Mast成熟亚型。 | [Tauber2023](https://pubmed.ncbi.nlm.nih.gov/37462672/) |
| Fib_contractile | Fibroblast | ACTA2, TAGLN, MYL9, POSTN, TPM2, COL11A1 | 无自动负向 | 先排mural/SMC；肿瘤背景才可称myCAF-like。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf); [Buechler2021](https://www.nature.com/articles/s41586-021-03549-5) |
| Fib_inflammatory | Fibroblast | IL6, CXCL1, CXCL2, CCL2, CXCL12, CXCL14, HAS1 | 无自动负向 | 多趋化子程序需保留差异，非所有CXCL12群均iCAF。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Fib_ECM | Fibroblast | COL1A1, COL1A2, FN1, MMP11, THBS2, COL3A1, MMP2 | 无自动负向 | ECM状态与基础身份相关，CAF需来源。 | [Buechler2021](https://www.nature.com/articles/s41586-021-03549-5) |
| Endo_tip_like | Endothelial | ESM1, APLN, KDR, DLL4, ANGPT2 | 无自动负向 | 血管生成/tip-like，不和血管分支互斥。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Epi_KRT19 | Epithelial | KRT19, KRT8, KRT18, EPCAM, MUC1 | 无自动负向 | 取消泛癌luminal身份断言。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| Epi_MSLN | Epithelial | MSLN, WFDC2, MUC1 | 无自动负向 | MSLN-associated，mesothelial需对照；非通用secretory身份。 | [Travaglini2020](https://www.nature.com/articles/s41586-020-2922-4) |
| EMT_like | Epithelial | VIM, FN1, ITGA5, LAMC2, COL17A1 | 无自动负向 | 上皮身份/同细胞证据及stromal混合另查。 | [Neftel2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/) |
| HSPC_myeloid_primed | HSPC | MPO, ELANE, CTSG, AZU1, PRTN3 | 无自动负向 | myeloid-primed，不证明离散GMP或正常细胞。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| HSPC_ery_primed | HSPC | GATA1, KLF1, GYPA, TFRC, PPBP | 无自动负向 | 分开红系/巨核方向，不证明双潜能MEP。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| HSPC_lymphoid_primed | HSPC | DNTT, IL7R, TCF7, LTB, CD7 | 无自动负向 | 需祖细胞身份，成熟淋巴对照。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| Ery_early_like | Erythroid | GATA1, KLF1, TFRC, CA1 | 无自动负向 | 可选分化阶段；TFRC非特异。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| Ery_late_like | Erythroid | AHSP, SLC4A1, HBB, HBA1, HBA2 | 无自动负向 | 环境HB另核，不强分早晚。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| Neut_immature_like | Neutrophil | MPO, ELANE, CTSG, AZU1, PRTN3, CEACAM8 | 无自动负向 | 成熟程度程序，HSPC/AML对照。 | [Velten2017](https://doi.org/10.1038/ncb3493) |
| Neut_OLR1_ARG1 | Neutrophil | OLR1, ARG1, S100A8, S100A9 | 无自动负向 | 取消PMN-MDSC终判；抑制功能需实验。 | [MDSC2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4935811/) |
| Neural_progenitor_like | Neural_glial | SOX2, NES, FABP7, VIM | 无自动负向 | 发育样/恶性状态另判，不能默认正常祖细胞。 | [Neftel2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/) |
| ASC_SDC1 | ASC | SDC1, MZB1, DERL3, PRDM1 | 无自动负向 | SDC1-associated，不单独证明成熟。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| ASC_secretory | ASC | MZB1, DERL3, XBP1, PRDM1, JCHAIN, SEC11C | 无自动负向 | 分泌程序不能独立区分plasmablast与plasma。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| Stress_response | * | HSPA1A, HSPA1B, DNAJB1, HSPH1, HSP90AA1 | 无自动负向 | 应激另轴，不在身份panel全局blacklist后丢失该读出。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| NK_adaptive_like | NK | KLRC2, LAG3 | FCER1G, SYK, SH2D1B, ZBTB16 | 仅候选；降低组件需同NK对照，不要求所有adaptiveNK都降低；移除TRDC正向。 | [Schlums2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4612277/); [YangNK2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6718415/) |

## isotype

| 模块 | 父类/适用 | 全部正向基因 | 期望降低（需对照） | 命名与边界 | 依据 |
|---|---|---|---|---|---|
| IgA | B, ASC | IGHA1, IGHA2 | 无自动负向 | 同种型直接转录证据；环境Ig、BCR和双细胞另核，不决定成熟/记忆身份。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| IgG | B, ASC | IGHG1, IGHG2, IGHG3, IGHG4 | 无自动负向 | 同种型直接转录证据；环境Ig、BCR和双细胞另核，不决定成熟/记忆身份。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| IgM | B, ASC | IGHM | 无自动负向 | 同种型直接转录证据；环境Ig、BCR和双细胞另核，不决定成熟/记忆身份。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| IgD | B, ASC | IGHD | 无自动负向 | 同种型直接转录证据；环境Ig、BCR和双细胞另核，不决定成熟/记忆身份。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| IgE | B, ASC | IGHE | 无自动负向 | 同种型直接转录证据；环境Ig、BCR和双细胞另核，不决定成熟/记忆身份。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |

## boundary

| 模块 | 父类/适用 | 全部正向基因 | 期望降低（需对照） | 命名与边界 | 依据 |
|---|---|---|---|---|---|
| T_NK_boundary | T, NK | CD3D, CD3E, TRAC, NCR1, KLRD1, NKG7 | 无自动负向 | 同细胞共享效应/真混合/错入父类分开；不叫NKT或污染。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf); [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
| B_ASC_transition | B, ASC | MS4A1, CD19, MZB1, PRDM1 | 无自动负向 | 连续分化或混合复核，绝不由CD79A/CD37判B污染。 | [YangB2024](https://doi.org/10.1016/j.cell.2024.06.038) |
| Mono_Mac_transition | Monocyte, Macrophage | FCN1, CD14, C1QA, APOE | 无自动负向 | 分化连续谱与混合两种解释分开。 | [Qian2020](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf) |
| Megakaryocyte_review | Megakaryocytic | ITGA2B, GP1BA, TUBB1, PF4, PPBP | 无自动负向 | 有核巨核/血小板需要核状态、复杂度/组织证据，不由cycling自动升级。 | [Velten2017](https://doi.org/10.1038/ncb3493) |

## reference_state

| 模块 | 父类/适用 | 全部正向基因 | 期望降低（需对照） | 命名与边界 | 依据 |
|---|---|---|---|---|---|
| Chu_Treg_functional_14 | T | FOXP3, IKZF2, IKZF4, IL2RA, ENTPD1, CCR4, ICOS, IL10RA, TGFB1, TIGIT, CTLA4, LAG3, HAVCR2, PDCD1 | 无自动负向 | Table S6原14基因完整保留；功能签名分数不替代Treg身份组合。 | [Chu2023](https://www.nature.com/articles/s41591-023-02371-y) |
