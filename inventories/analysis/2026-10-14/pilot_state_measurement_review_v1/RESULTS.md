# 冻结试点程序的测量与身份内变化审阅

全部98programme、495source-major记录保留，不选最有利K/基因，不新增LR/临床检验，不正式命名states。由已审source-major/source-node量化汇总派生，495支持比例和原node内部spread独立标准库复算通过；没有在本步骤再次逐cell重算UCell。

|K15示例程序（其他K仍全保留）|五来源中检测>=5签名基因的细胞比例范围|加权原node内q90-q10 / 来源大类q90-q10的中位数|abs RNA rho最大值|
|---|---|---|---|
|T Treg-related|97.39—99.99%|0.710|0.597|
|T cytotoxic-related|61.16—92.99%|0.826|0.393|
|T G2M-related|15.10—55.75%|1.026|0.528|
|Mac IFN-like|97.59—100%|0.919|0.589|
|Mac T-core/cytotoxic-related|4.12—8.80%|1.009|0.095|
|Mac G2M-related|54.30—80.25%|1.180|0.430|

**检测>=5仅代表测量支持，不是状态阳性比例。** spread比值不是方差分解，可以超过1；不据此证明独立状态或机制。RNA相关描述技术/生物关联，不能判因果。

已有程序在原subcluster内有变异，说明只看细分身份/数量可能遗漏表达变化；是否可正式叫状态还需签名、RNA深度、谱系身份、来源和跨方法审阅。T细胞的cytotoxic也可能含CD8/NK/γδ亚型差异；Treg-related不能将所有高分T重新标为Treg。

新增“other_immune_lineage_signature_review”只在审阅表标记：Mac top50含>=3T_core、或T含>=3Myeloid_core，基于已冻结panels。原字典/评分/稳定性/foreign_RNA字段完全不变。这是混合RNA/身份/功能歧义提醒，不能自动删除细胞；已有Mac T-core签名的CD3/TCR相关基因需核查而非直接命名cytotoxic macrophage。

所有缺TPM alternate（Spectra没有该对照）及恒定分数相关继续NA，汇总nanmedian告警仅为空对照，不填0、不重新运行数据。没有通过此表晋升任何状态、调整阈值或重新检验LR。
