# ccRCC现有身份：CM邻域与空间通讯

用户要求身份疑点先记录，按现有身份继续。因此主表固定为原`spatial_pilot_v1/identity_candidates_PRIVATE.tsv.gz`，SHA见[保留回执](../inventories/analysis/2026-11-09/ccrcc_followup_v1/current_identity_forward_v1/identity_preservation.json)。不复用新增聚类/一致性掩码，不做新亚型冻结，状态分析不作前置。

输入：已验收的同一477 RNA/27蛋白源及原358,841 QC细胞；空间原生质心微米不缩放。RNA和蛋白按feature_type独立，通讯只用RNA。6类原身份151,833细胞进入CellChat，其余原QC细胞仍参与最近邻几何背景。

## 空间组成完整12槽

每个B/T/NK接收细胞取最近6个原QC细胞，排除自身，截50/100μm。统计EC—B边、EC—NK边及同时有B和NK接收邻居的EC数；后者是原CM两边以外的新组织结构假设。250/500μm网格×floor(log2 RNA UMI)内置换B/T/NK身份999次；EC/位置不动，局部各免疫身份数量固定。3关系×2半径×2网格=12，单侧(1+null≥obs)/1000，统一BH/BY。

## 空间表达完整108槽

发送类Blood_EC/Mural/Lymphatic，接收B/T/NK，轴CCL19–CCR7、CXCL12–CXCR4、CCL21–CCR7；2半径×2网格共108。观测为6NN边上配体/受体均≥1原生RNA UMI的比例。以4个可测基因的整行二元profile，在身份×网格×RNA深度层内置换，保留同细胞多基因结构与局部阳性量。999次；整数边数≥观测含并列，加一P，缺测保留NA/P1参与全108 BH/BY。它检验局部空间表达耦合，不是功能信号/因果激活。

## 空间CellChat

复用既有2.2.0环境和数据库。原身份六类、原477RNA CP10K/log1p，逐barcode对齐原生微米坐标，一份组织几何；不重新训练分类器、不改QC/标签、不做PPI补值。数据库配体、受体复合体及非空调节亚基需全部实测；仅蛋白信号。truncatedMean trim0.1、do.DE FALSE/min.cells10、raw.use TRUE、population.size FALSE、k.min10、nboot100、固定seed，interaction.range100/250×contact.range10/20共4变体。geometry-only scale.distance=1.5/最小正组间距离；meta同时有labels与group，沿用已验证接口。[CellChat官方空间教程](https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.Rmd)

保留原生prob/P及其原生候选表。另从保存的100次超越计数给出(1+b)/101的有限MC敏感性，并对所有模型槽/四变体共同BH/BY；原工具的并列规则无法由保存P重建，因此此列不冒称新的精确空间零模型，原生P=0也不作确定证据。真正的局部空间背景检验是上述108家族。保存数组独立读回、候选计数及全模型FDR核对后才记录模型验收。

## 记录与范围

身份疑点登记为非阻断事项；此前独立EC门槛、RNA细分和一致性视图留档，不作为主身份或主结果。T细簇状态分裂不能当作大类身份不一致。当前不重新注释，不把这些问题当启动通讯的条件。

同一供者、同一切片；细胞/EC/半径/网格/四个变体不是患者重复。有限面板缺IGHD/XCL1/CCL21，不能声称精确IGHD-B或XCL1-NK模块及CCL21信号复现。现有身份的潜在混合会影响解释，按用户要求后续独立复核。

实际结果：[ccRCC主报告](../inventories/analysis/2026-11-09/ccrcc_followup_v1/README.md)。代码在`scripts/ccrcc_followup_20261109/`；主入口`ccrcc_current_identity_forward_v1.py`，未来CellChat接续用`run_current_identity_CellChat_v2.py`（等待不创建生产目录，避免启动竞态）；v1控制器历史保留，现有模型不重启。矩阵/细胞/边/坐标/模型仅在${RUN_ROOT}。
