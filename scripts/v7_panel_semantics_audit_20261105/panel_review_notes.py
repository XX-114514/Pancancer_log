"""Manual semantic review, not a replacement marker classifier."""
NOTES = {
'Tcell': {
'CD4_naive':'panel 不含 CD4；IL7R/CCR7/TCF7 不能区分 CD4 与 CD8，也不能独立分清 naïve 与 central-memory。先确认 CD4 分支，再用 naïve-like 描述并检查记忆/激活证据。',
'CD4_Tfh':'CXCL13/PDCD1/ICOS 可与激活或耗竭程序共享；必须确认 CD4，并比较 BCL6/CXCR5 等完整 Tfh 程序及组织背景。仅 CXCL13 不足以认定 Tfh。',
'CD4_Th1':'先确认 CD4；TBX21/IFNG/CXCR3 也见于 CD8/NK。保留 Th1-like 程序候选，不能由炎症程序倒推 CD4。',
'CD4_Th17':'先确认 CD4；RORC/CCR6/KLRB1 与 MAIT/ILC 共享。需成套程序及非常规 T 对照，不以 KLRB1 单基因命名。',
'CD4_Treg':'保留 Treg 候选，但 FOXP3/IL2RA/IKZF2 应合看，CTLA4/TIGIT 是共享激活/调节证据；确认 T/CD4 分支与群内一致性。Tcell_6 为实际漏评案例。',
'CD8_naive':'虽含 CD8A/B，canonical 却不要求 CD8；CCR7/SELL/TCF7 与 CD4-naïve 共享。CD8 身份单独确认，naïve-like 程序单列。',
'CD8_GZMK':'panel 不含 CD8A/B；GZMA 可替代 GZMK 成为唯一 canonical。先确认 CD8，实际有 GZMK 特征才用 GZMK 命名；不能与 NK 细胞毒程序互斥竞分。',
'CD8_cytotoxic':'panel 不含 CD8A/B；CCL5/NKG7/GNLY/PRF1 不能确定 CD8，CD4 细胞毒 T 也可能表达。分开 CD8 身份和细胞毒状态。',
'CD8_exhausted':'panel 不含 CD8A/B；TIGIT/CTLA4 相关调节程序不能替代 CD8/耗竭证据。确认分支后联合 PDCD1/HAVCR2/TOX 等程序，使用 exhausted-like 并避免功能实证外推。',
'Tissue_resident':'CD69 可反映激活，CXCR6/ZNF683 等支持驻留样程序但不直接证明组织驻留；标记 resident-like 状态并保留 CD4/CD8/非常规身份。',
'MAIT':'KLRB1/ZBTB16/SLC4A10 可提名 MAIT-like；不能仅靠共享基因确认 MAIT。结合 TCR（可用时）、TRAV1-2 与完整程序及 γδ/常规 T 对照。',
'gamma_delta':'TRDC/TRGC1/TRGC2 比细胞毒程序更直接，保留候选；检查同细胞 TCR 表达、αβ T 背景和可用 VDJ，低覆盖不能直接视为阴性。',
'NK_contamination':'GNLY/NKG7/KLRD1 与细胞毒 T 共享；移出污染/亚型竞赛，作为 cytotoxic/NK-associated 程序。T 身份强时不得据此否定 T；NK 转归须重新评估 TCR/CD3 与 NK 证据。'},
'NK':{
'FCGR3A_cytotoxic':'确认 NK 父身份后描述 FCGR3A/FGFBP2 相关细胞毒程序；FCGR3A/GNLY/PRF1 不独占 NK，不能压过强 TCR 证据。',
'XCL1_CD56bright':'XCL1/XCL2/KLRC1 不等价于 CD56 蛋白 bright。RNA-only 可标 XCL1/XCL2-enriched NK 或 CD56bright-like 候选，若直接称 bright 需相应蛋白证据。',
'Tissue_resident':'CD69/ITGAE/CXCR6/ZNF683 是驻留样候选程序；需确认 NK，排除 T/ILC 分支混用，不能从解离 RNA 推断已证明驻留。',
'Adaptive_KLRC2':'严重方向/身份问题：FCER1G 被列为正向 canonical，但经典 adaptive NK 子集常降低 FCER1G；TRDC 是 γδ T 线索。暂停 adaptive 命名，改为 KLRC2-associated 候选并核查方向性签名和 TCR。',
'NKT_like':'CD3D/CD3E/TRAC 加 NKG7/KLRD1 不能在 NK 父类内证明 NKT，更不能证明 iNKT。作为 T/NK 边界复核，必要时回到 T 分支和 TCR/CD1d 证据。'},
'ILC':{
'ILC1':'TBX21/IFNG/KLRD1 与 T/NK 共享，需先确认非 T/B、非典型 NK 的 ILC 身份；不能仅看 Th1 程序。',
'ILC2':'IL7R/KLRB1 非特异；联合 GATA3/PTGDR2/IL1RL1 等并核查 T、Mast、basophil 对照。候选 ILC2 不等于已确认。',
'ILC3':'KIT/KLRB1 不足以确认 ILC3；需 RORC/IL23R/CCR6 等一致程序及 T/Mast 对照。ILC_3 当前只有共享基因通过。',
'NK_like':'ILC 和 NK 有生物关系，但 NK-like 不是可与 ILC1 互斥判定的终点；保留边界状态，核查父身份。'},
'Bcell':{
'Naive':'MS4A1/CD79A/B 是 B 身份基因；MS4A1 不应独立充当 naïve anchor。需 IGHD/TCL1A/IGHM 等与 memory/GC 对照；IGHD/IGHM 也非绝对特异。',
'Memory':'CD27 不是所有记忆 B 必备，IGHG1 是同种型证据；避免混合记忆、类别转换和激活三条轴，比较 naïve/GC/ASC 程序。',
'Germinal_center':'BCL6/AICDA/RGS13 组合可保留 GC 候选；验证 B 父身份、样本一致性和非循环状态导致的差异。',
'Age_associated':'TBX21/ITGAX/FCRL5 可支持 atypical/ABC-like 程序，不能从表达反推年龄来源；用 TBX21/FCRL5-associated B 更中性。',
'Plasma_contamination':'ASC 分化和 B→plasma 连续谱不等于污染；XBP1 是分泌/UPR 线索，需 MZB1/JCHAIN/PRDM1/SDC1 等及成熟 B 对照。'},
'Plasma':{
'IgA':'同种型轴；JCHAIN 不是 IgA 专属，必须直接支持 IGHA1/2 才能标 IgA。与成熟/增殖状态可同时成立。',
'IgG':'同种型轴；需直接 IGHG 证据，且旧 panel 只含 IGHG1/3，遗漏其他 IgG 亚类不等于阴性。',
'Mature_SDC1':'MZB1/DERL3 是分泌细胞程序，不能独立证明成熟或 SDC1 阳性；需 SDC1 实际支持并结合分化背景。',
'Plasmablast':'MZB1/JCHAIN/PRDM1 与成熟 plasma 共享；仅分泌程序不能区分 plasmablast/plasma。结合分化/增殖、蛋白和 BCR 证据（可用时）。',
'Bcell_contamination':'CD79A/CD37 为共享 B 谱系线索，不能据此否定 plasma；检查 MS4A1/CD19 等成熟 B 程序、ASC 程序与同细胞共表达。'},
'Monocyte':{
'Classical_FCN1':'联合 FCN1/S100A8/A9/CD14，区分 neutrophil、DC3 和 macrophage 连续谱；名称含 FCN1 时需其实际支持。',
'Nonclassical_FCGR3A':'FCGR3A/MS4A7 也见于 macrophage/NK；须有 monocyte 身份，并比较 CD14 与非经典整体程序。',
'Intermediate':'不能由 HLA-DRA 或 CD14/FCGR3A 任一单项确认；需要同群/同细胞 CD14 与 FCGR3A 及连续分化背景。',
'Inflammatory_IL1B':'炎症状态，不独立定义 monocyte 亚型；IL1B/CXCL8 需实际支持，与 IFN/分化类型可共存。',
'Macrophage_contamination':'单核到巨噬细胞分化是连续过程；重命名为 macrophage-associated evidence，判断父类转归/过渡/混合，不能一律污染。',
'DC_contamination':'多个 DC 分支混在同一 panel，不能要求各支全部共表达；在髓系子集内分支比较并核查 CD1C、CLEC9A、LILRA4 等不同模块。'},
'Macrophage':{
'C1QC_APOE':'可保留 C1Q/APOE-associated 描述；不能由 APOE 单项证明 C1QC 型或免疫抑制功能。',
'SPP1_TREM2':'SPP1/TREM2/GPNMB 程序需与 macrophage 父身份分开；名中双 marker 不代表两者都被证明。',
'MRC1_CD163':'可描述 MRC1/CD163-associated，不应直接等同 M2 或抑制功能；检查实际 marker 和组织背景。',
'Inflammatory_IL1B':'炎症状态轴，不能与 C1Q/SPP1 身份候选硬互斥。',
'MHCII':'抗原呈递相关状态，不能作为独立 macrophage 身份，也不证明实际呈递功能。',
'Monocyte_contamination':'FCN1/S100A8/A9 可为分化/炎症信号；优先检验连续谱与群内混合，不直接定污染。',
'DC_contamination':'需髓系内 DC 分支与巨噬细胞证据联合比较；共享呈递/分化程序不是 doublet 证明。'},
'Dendritic':{
'cDC1':'CLEC9A/XCR1/CADM1 组合可保留；BATF3/IRF8 的共享表达不可独立确定 cDC1。',
'cDC2':'CD1C/FCER1A/CLEC10A 可提名；与 DC3/monocyte/basophil 等比较，不能只凭 CD1C。',
'DC3':'CD1C 与 cDC2 共享；S100A8/CD14 与 monocyte 共享。需要完整组合和 cDC2/monocyte 对照，不能凭单个 canonical。',
'pDC':'联合 LILRA4/CLEC4C/TCF4 等；GZMB/JCHAIN 与 NK/plasma 等共享，GZMB 单项不能定义 pDC。',
'Migratory_LAMP3':'成熟迁移相关 DC 状态可保留，但不能等同独立发育来源；同时保留可确定的 cDC 身份。',
'MHCII':'DC 多支普遍表达的功能程序；不能因 MHCII 高就定义新的 DC 亚型。',
'Myeloid_contamination':'DC 属于髓系框架，LYZ/LST1/S100A8 不构成通用反证；转为 monocyte-associated/炎症边界检查。'},
'Granulocyte':{
'Neutrophil_mature':'FCGR3B/CSF3R/CXCR2 与成熟谱系程序联合，避免 S100A8/A9 共享信号单独命名。',
'Neutrophil_immature':'MPO/ELANE/CTSG 是髓系颗粒程序；需与 GMP/髓系祖细胞、恶性髓系细胞比较后判成熟阶段。',
'PMN_MDSC':'PMN-MDSC 含免疫抑制功能主张，不能仅用 OLR1/ARG1 转录证明；RNA-only 标 OLR1/ARG1-associated neutrophil 或 PMN-MDSC-like 候选。',
'Eosinophil':'CLC/EPX/PRG2 组合可保留，需 basophil/mast 对照和源数据覆盖核查。',
'Myeloid_contamination':'粒细胞本属髓系；LYZ/S100A8/A9 是共享信号。该名称逻辑不当，改为 monocyte/macrophage boundary review。'},
'Mast':{
'TPSAB1':'TPSAB1/TPSB2/KIT/HDC 多为 Mast 大类身份证据；不能因没有谱系内富集就否认 Mast，也不自动构成亚型。',
'CPA3':'CPA3/HPGDS/HDC 与其他 Mast 程序共享；CPA3 不等同 CMA1，不得将此 panel 当 MCTC 蛋白分型；需实际差异和共表达。',
'Activated_IL1RL1':'名称含 activation 和 IL1RL1，但 HDC/LTC4S 可通过而 IL1RL1 不通过。去除未经支持的 activated 主张，优先 marker-enriched 描述并单独验证激活程序。',
'Basophil_like':'basophil 应作为邻近谱系边界复核，不是当然的 Mast 子型；CLC/MS4A2/HDC 共享，需 KIT/TPSAB1/B2 与 IL3RA 等组合。'},
'Epithelial':{
'Basal_squamous':'组织/分化表型，需来源背景及 basal 程序；不能推断恶性。',
'Luminal_KRT19':'KRT8/18/19/EPCAM 是广泛上皮程序，不足泛癌通用 luminal 分类；结合组织特异谱系。',
'Secretory_MSLN':'MSLN/WFDC2/MUC1 不是跨癌通用 secretory 身份；还需 mesothelial 等对照。',
'Alveolar_AT2':'肺泡表型组合可保留，需肺组织和同群证据；肿瘤中是 AT2-like 分化，恶性单独判断。',
'Alveolar_AT1':'CAV1/EMP2 等非特异；联合 AGER/PDPN 与肺背景，排除血管/基质。',
'Club':'SCGB1A1/SCGB3A1 组合可保留为 club-like，需气道背景并检查分泌物环境 RNA。',
'Ciliated':'FOXJ1/TPPP3/PIFO 等组合可保留，结合纤毛程序与组织来源。',
'Hepatocyte_like':'ALB/APOA1/TTR 易受环境 RNA/肝来源影响；需完整肝细胞程序及同细胞表达。',
'Cholangiocyte':'KRT7/KRT19/EPCAM 不是胆管专属；必须有肝胆背景和其他谱系对照。',
'Renal_tubular':'EPCAM 不能充当肾小管特异 anchor；需 SLC34A1 等分段/来源证据；KIM1→HAVCR1 的历史映射保留审计。',
'Prostate_luminal':'KLK3/ACPP/NKX3-1 需前列腺背景和同细胞证据；ACPP/ACP3 等版本别名按 stable ID 审核。',
'Melanocytic':'melanocyte 是独立谱系，不应塞入 epithelial 子型；PMEL/MLANA/TYR 与恶性判断分开。',
'EMT':'状态程序；VIM/FN1 也可能来自 stromal 混合，需保留上皮/肿瘤身份和同细胞证据，不能单独定义谱系。'},
'Endothelial':{
'Blood':'层级比 arterial/venous/capillary 更粗；广泛 endothelial marker 不能与这些子型同级竞分。',
'Lymphatic':'PROX1/LYVE1/CCL21 组合结合 endothelial 身份，避免 LYVE1 macrophage/PDPN fibroblast 单项误导。',
'Tip_angiogenic':'血管生成/tip-like 状态，可能横跨血管分支；KDR 单项不足。',
'Venous_ACKR1':'PLVAP 不特异静脉也不替代 ACKR1；标 ACKR1 时必须有其证据，SELE/VCAM1 激活程序另列。',
'Arterial':'EFNB2/SOX17/GJA5 等联合可保留，需同谱系对照与跨来源一致性。',
'Capillary':'RGCC/CA4/CD36 有组织差异，IFI27 可能状态依赖；需来源分层，不用一个通用 panel 覆盖所有器官毛细血管。'},
'Fibroblast':{
'myCAF':'ACTA2/TAGLN/MYL9 与 pericyte/SMC 共享；先排 mural/SMC，确认肿瘤背景后称 myofibroblastic CAF-like。',
'iCAF':'炎症/趋化程序，可与 ECM/myofibroblastic 程序共存；不能把所有 CXCL12 正常 fibroblast 自动叫 CAF。',
'apCAF':'MHCII/CD74 表达支持呈递相关程序；确认 fibroblast 同细胞身份、排除 APC 混合，功能呈递仍需额外证据。',
'ECM_CAF':'COL1A1/FN1 属广泛 fibroblast/ECM 程序；需区别基础身份与肿瘤特异状态，不能自动称 CAF。',
'Pericyte':'pericyte 是 mural 邻近谱系，不是必然 fibroblast 子型；重建父类树，检查 RGS5/PDGFRB/CSPG4 与 fibroblast 证据。',
'Smooth_muscle':'SMC 是邻近谱系，应单独父类/分支审查；MYH11/CNN1/DES 联合，区别 myCAF。',
'Mesothelial':'mesothelial 是邻近谱系而非通用 fibroblast 子型；MSLN/CALB2/KRT 与来源结合，不自动判 epithelial contamination。'},
'Erythroid':{
'Early':'分化阶段轴；TFRC 不专属红系，需 GATA1/KLF1/GYPA 等大类身份后判阶段。',
'Late':'HB 易环境污染；需 AHSP/SLC4A1/GYPA 组合与细胞质量，早晚分期不强制。',
'Megakaryocyte_contamination':'红系/巨核祖细胞可相关，PPBP/PF4 也易环境污染；先区分巨核分化、血小板粘附、环境 RNA 与 doublet。'},
'Platelet':{
'Platelet_PPBP':'主要是大类身份，不必强行亚型；PPBP/PF4 环境污染和 platelet-cell aggregate 需同细胞/QC 证据。',
'Megakaryocyte':'与 platelet 共享基因，需结合核/转录复杂度、组织和完整巨核程序；不能仅凭 PF4/PPBP 区分。',
'Cycling_Megakaryocyte':'增殖程序不能把 platelet 自动升级为 megakaryocyte；先证实有核巨核/祖细胞身份，再单列 cycling。'},
'HSPC':{
'HSC_MPP':'CD34 与 endothelium 共享，需 SPINK2/HLF/MEIS1 等及造血背景；白血病来源不能自动当正常 HSC/MPP。',
'GMP':'MPO/ELANE/CTSG 也见于 immature neutrophil/恶性髓系；需祖细胞身份后再判 myeloid-primed。',
'MEP':'GATA1/KLF1/GYPA/PPBP 偏红系/巨核程序，不能独立证明双潜能 MEP；先用 erythroid/megakaryocytic-primed 候选。',
'Lymphoid_progenitor':'IL7R/TCF7/LTB 见于成熟淋巴；需祖细胞程序/DNTT 与成熟 T/B 对照。',
'Myeloid_contamination':'HSPC 髓系分化或白血病程序不等于污染；联合原始病种、祖细胞证据及 myeloid 程序解释。'},
'Neural':{
'Oligodendrocyte':'MBP/PLP1/MOG 支持少突谱系，但旧 Neural 大类偏 Schwann；需补齐父类覆盖并独立恶性判断。',
'Schwann':'S100B/SOX10 与 melanocyte 共享；联合 MPZ/PLP1，并以 PMEL/MLANA/TYR 程序审查边界。',
'Neural_progenitor':'SOX2/NES/VIM/FABP7 不证明正常神经祖细胞；胶质肿瘤/其他肿瘤也可能表达，需 CNV/来源与谱系证据。',
'Neuron':'SNAP25/RBFOX3/SYT1/TUBB3 组合可保留，但旧 Neural 大类不覆盖典型 neuron；父类召回需审查。',
'Melanocytic':'与旧 Neural 排除 PMEL/MLANA/TYR/DCT 的父规则冲突；应设独立 melanocytic 谱系，不能继续放在 Neural 亚型。'}
}

MAJOR_NOTES={
'Tcell':'保留 CD3D/CD3E/TRAC 身份模块；LCK/CD247 与 NK 共享，TYROBP 不宜一律作 T 反证。旧代码只要求来源含≥2基因，并未要求每个细胞实际检测≥2身份基因。',
'NK':'旧正向面板主要细胞毒程序，需 NK 受体/谱系组合及 TCR/CD3 反证；排除面板缺 TRAC，不应以 RNA dropout 自动认定 TCR 阴性。',
'Bcell':'CD74 为 APC 共享，CD37/CD79A 与 plasma 共享；需成熟 B 与 ASC 分支而非互相污染。',
'Plasma':'ASC 程序较合理；MS4A1 需结合分化/群内混合，不能忽略 plasmablast 连续谱。',
'Monocyte':'S100A8/A9 与 neutrophil 共享，LILRB1/CTSD 广泛髓系；应联合 monocyte 特异组合。',
'Macrophage':'C1Q/APOC1/MARCO 偏成熟/部分状态，可能漏检炎症/单核来源群；补充髓系树与不同状态覆盖。',
'Dendritic':'将 cDC1/cDC2/pDC 标志混合平均会稀释某分支；旧 pDC 仅 GZMB 较直接，需分支 OR 结构后合并 DC 父身份。',
'Mast':'KIT/TPSAB1/B2/CPA3/MS4A2 组合可保留，核查 basophil 与 KIT+祖细胞；反证作为边界/QC，不等同 doublet。',
'Granulocyte':'旧面板偏 neutrophil，不能覆盖 eosinophil 等所有粒细胞；父类名应与实际覆盖一致。',
'Epithelial':'泛角蛋白/EPCAM 覆盖偏 epithelial；不覆盖 melanocyte，EMT 或低角蛋白肿瘤需独立判断。',
'Endothelial':'KDR/VWF/PECAM1 等有分支差异，补充 lymphatic/organ-specific coverage；正负均需来源可测性。',
'Fibroblast':'ECM 正向面板偏常规 fibroblast，不涵盖 mural/mesothelial；PTPRC/EPCAM 等仅作边界证据。',
'Erythroid':'三个 HB 高表达可主导分数，需 AHSP/GYPA 等非 HB 证据；PTPRC 单项排除要区分有核前体/混合/环境 RNA。',
'Platelet':'PPBP/PF4 易环境 RNA/附着；需组合、QC 与复杂度；可将大类作为合理终点。',
'HSPC':'CD34/GATA2 不独占 HSPC，需祖细胞组合和血液疾病背景；旧 vascular guard 只降低假阳性，不证明剩余正常。',
'Neural':'旧大类偏 Schwann/glia，不能作为全部 neuron/oligodendrocyte 的共同正向 panel；与 Melanocytic 子型定义冲突。',
'ILC':'IL7R/KIT/KLRB1/RORC/IL1R1 偏部分 ILC，并与 T/Mast 共享；须分支 OR、TCR/NK/Mast 对照，父类本身也需复核。'}

def review(major,name):
    specific=NOTES.get(major,{}).get(name)
    if name.startswith('Cycling_'):
        return 'state_program','分离状态',specific or 'MKI67/TOP2A/UBE2C/STMN1 是增殖程序，不提供谱系身份。先确认父类，再单列 cycling；不与身份亚型争夺唯一标签。'
    if name.startswith('IFN_'):
        return 'state_program','分离状态','IFN 响应是可跨谱系共存的状态，不能提供父身份，也不应阻止已有亚型；作为独立程序输出。'
    if 'contamination' in name:
        return 'boundary_review','停止自动污染结论',specific or '此 panel 只能提出邻近谱系/异源表达候选。需父身份、同细胞共表达、子集结构、环境 RNA、doublet 与样本来源证据，不能根据字符串或共享 marker 直接定污染。'
    assert specific is not None,(major,name)
    if name in ['NKT_like','NK_like','Melanocytic','Pericyte','Smooth_muscle','Mesothelial','Basophil_like']:
        role='lineage_boundary'
    elif name in ['Tissue_resident','EMT','MHCII','Inflammatory_IL1B','Tip_angiogenic','Migratory_LAMP3','Activated_IL1RL1','PMN_MDSC','Adaptive_KLRC2','myCAF','iCAF','apCAF','ECM_CAF']:
        role='state_or_context_program'
    elif major=='Plasma' and name in ['IgA','IgG']:role='isotype_axis'
    elif name in ['Blood','TPSAB1','CPA3','Platelet_PPBP']:role='identity_or_descriptive_program'
    else:role='conditional_identity_or_differentiation'
    return role,'条件保留或重构',specific
