# 多初始化状态程序：第一轮真实诊断

## 实际结果

2026-10-14 01:45:28—01:45:37 UTC（09:45 UTC+8）真实诊断exit0；01:47:41独立审计完成。4来源均衡12000细胞、3seed、173因子/6512特征，所有输出SHA及完整键已核查。519匹配cosine/Jaccard、957评分相关的1914统计、638来源/身份R2和全部173稳定规则独立复算通过；cosine最大误差1.95e-7，R2最大1.23e-7。数值验收不等于生物状态验收。

[安全描述候选表](all173_candidates_interpretation_safe.tsv)、[完整匹配](pairwise_factor_matching.tsv)、[原身份/复杂度诊断](source_identity_complexity_diagnostics.tsv)、[审计](independent_audit.json)。预先规则[计划](../Spectra_multiseed_diagnostics_plan_v1.json)，全部结果保留，未选择最有利seed。

|固定门槛|通过数 /173|
|---|---:|
|3对权重cosine均≥.9|45|
|3对top50Jaccard均≥.5|32|
|适用major去来源评分相关均≥.7|33|
|3对匹配cycle一致|156|
|上述条件同时满足|8|

8个探索稳定候选：1shared、5T、2Mac。**这不是8个最终状态，也不是165个没有生物意义的程序。** 初始prior多、factor多、minibatch、所选细胞和有限优化均可能影响结果；目前仅同数据不同初始化稳定性，没有完成来源留出泛化或cNMF交叉方法验证。

## 为什么还不能正式命名

- SP_T_cell_150的top genes包括KLRB1/CCR6/RORC/IL17A等，可考虑TH17相关程序，但现有broad T注释不能据先验直接变成正式TH17细胞，仍需检测覆盖、原亚型内变化和迁移核查。
- SP_T_cell_151的参考重叠为IL2–STAT5集合，LTB/KLF6/RORA/CTLA4/ICOS等不构成唯一IL2刺激证明；其source-nodeR2最高.283。
- SP_T_cell_161和SP_Macrophage_168虽稳定，其top genes含基础转录/蛋白合成或伴侣蛋白；参考集合名称IL12/LPS不能直接作为特异刺激标签。
- SP_T_cell_166和SP_Macrophage_173两个稳定自由/候选因子存在复杂度标记，source-nodeR2分别.537/.472；可能包含身份及生物活动，不能因为相关直接删除或称技术伪影。
- SP_global_021偏糖胺聚糖合成；sourceR2最高.182，仍不意味着跨来源机制验证成立。

v1描述字段发现3条top50与所有prior零交集时选了字典序名称。原v1结果保留，新安全表将其改为NO_PRIOR_OVERLAP_IN_TOP50并明确descriptor_is_confirmed_state=false；稳定规则/拟合/统计完全不改，不把Treg或LPS误赋给零交集因子。

## 对此前通讯结果的解释

此前不是全部通讯筛查都“失败”：原CM完整20000次置换有178个主q记录，147保留区间，主要CRC；只是尚未达到同明确端点两来源共同q/区间和独立机制验证。[原结果](../../2026-10-12/CM_LR_permutation_precision_v2/RESULTS.md)。

证据链断点包括：细分端点共现/可比患者不足；组成共变不能直接推LR机制；丰度和同群细胞数影响信号；共享刺激/接收自身RNA使原8主BH线索变为0；参考mean-log程序67.8%记录与复杂度高度相关；全局删LR又损失真实程序基因；IFNG读出/控制敏感且来源信息集中；空间panel和接收程序不足。控制可能同时是混杂或中介，衰减不等于不存在生物关系。[程序准备诊断](../../2026-10-13/state_programme_readiness_v1/RESULTS.md)、[细分有限两轮](../../2026-10-13/CM_fine_IFNG_finite_round_v1/RESULTS.md)。

来源/癌种限定生态和多癌种同向效应有研究价值。保留原9CM/172身份与数字，不以换检验、阈值或新方法抹去历史；新状态问题优先研究同一身份内活动及多细胞协调，再检查候选接口是否解释接收响应。

## 下一步

cNMF完整600run+全部K consensus继续原句柄，T部分完成、Mac进行中；比较稳定性/重建及身份/深度，不按LR或临床p挑K。之后固定多基因模块，用每来源完整实测gene轴UCell评分，标记未测/低覆盖，连续/多状态sidecar独立于身份。原9CM其他谱系必须补充，随后状态相关细胞对/LR和接收响应，状态完成后DIALOGUE及新空间核查。全部尚未称完成。

主要方法依据：[Spectra原研究](https://doi.org/10.1038/s41587-023-01940-3)结合类型与过程先验学习程序；[cNMF原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC6639075/)以多次分解共识发现身份和活动。以上论文支持方法选择，不验证本数据的任何状态或通讯轴。
