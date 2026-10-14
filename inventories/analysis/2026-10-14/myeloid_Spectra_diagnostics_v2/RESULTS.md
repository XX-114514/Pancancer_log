# 髓系5大类完整三seed诊断

26723cells，Mac/Mono/Myeloid/DC/Neutrophil五原major、4575features、177factors；原已完成第一seed逐SHA复用、后两seed同参数data1续算。06:22:41—06:24:13 UTC固定诊断完成，06:24:13—06:24:34独立531factor匹配、2283major-score关系、761回归/原生支持、177规则全部通过。

{"factors": 177, "cycle": 155, "mincos9": 11, "minJ50": 9, "weightpass": 1, "score_anymajor": 8}

最终0shared/lineage/context合格候选。唯一权重稳定AS_myeloid_148，scopeMac，mincos.908/top50Jaccard.754，但最低source-residual score rho.4056<.7。其prior标签Mac_LPS_response不代表已确认LPS刺激或state。raw数据/身份/score路径与旧原模型保留，不改稳定性门槛。

这是当前分解未达到预设重复稳定性要求，不能把整个髓系数据判坏或把生物programme判不存在；继续无先验cNMF和原身份内/原生rank测量。Source/cancer限定效应仍有价值，source-high不自动artifact。数值audit PASS不升级状态或通讯机制。

独立source分布quantiles本步骤未再次计算；无cell级p或LR/临床p筛选。其余4fit blocks/原9CM/完整stateUCell/DIALOGUE/空间仍待。
