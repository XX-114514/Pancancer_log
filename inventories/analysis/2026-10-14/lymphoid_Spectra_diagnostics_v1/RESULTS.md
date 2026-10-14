# 首T/淋巴块完整三种子诊断：固定条件未通过

15,644cells/3majors/207factors/4984features；三seed真实完成05:01:07 UTC。固定诊断05:03:22—05:04:20完成；05:08:25—05:08:40独立621factor匹配、1503score关系、501major-factor回归/原生支持和全部207资格规则通过。

207中cycle187、mincos>=.9为25、minJaccard>=.5为18、cycle+cos+Jaccard4；49个factor至少1major的score/基因支持达到条件，但与上述4无交集。最终0shared/lineage主候选，0context候选。**这是当前模型划分未达到预设重复稳定性要求，生物programme是否存在仍待其他方法和固定签名测量。** 不改旧五来源8/173诊断，也不通过此结果调整门槛。

4个权重稳定因子：AS_147 global（各major最低source-residualrho.158—.305）；AS_152 T（prior IL2/STAT5描述、rho.662）；AS_161 T（TCR描述、rho.670）；AS_162 T（IL12描述、rho.528）。prior描述不是确认的刺激/state；具体genes全保留。来源效应/身份/RNA深度和基因支持已分列。3major中支持>=5genes的中位fraction接近1，但该fraction不是state阳性比例。

独立复算采用dummy-lstsq、scalar Spearman和sparse incidence；source分布quantiles本步骤未再次独立复算。没有cell级p或LR/临床p筛选。其余谱系、跨method和全原9CM状态目标仍待。
