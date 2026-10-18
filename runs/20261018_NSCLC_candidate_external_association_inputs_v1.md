# 原候选完整外部来源模型输入

2026-10-18 09:27:11.789056 UTC /17:27:11.789056 UTC+8在server私有冻结d4141c94fe11d2fbccd764fbc2bff503d8a068e7e32132cb6f5c790bb158d752。已绑定真实RNA c705d400和原评分c04e2169，尚未读取新外部关联/P或临床结局。[完整合同](../methods/NSCLC-candidate-external-source-association-v1.md)。

新源码prepare_NSCLC_candidate_external_associations_v1.py重建全部同文库机会，保留全部原计划。旧程序—程序名单要求双方整套programme的source gate，本问题没有用发送者程序，故不继承该额外门槛；接收原评分两个gate、RNA/质量/父比例/全分类块保持。逐文库保存所有失败flag，逐单位保存同病例聚合、原paired library键、RNA分母和完整控制，逐模型保存资格与可识别性及冻结家族大小；不计算效应/P。

代码检查前发现旧assay单位表没有annotated_scope_cells列，因此从原已接受完整library metadata按严格相同轴读取全注释细胞分母，并核metadata SHA与刚接受RNA一致。此为部署前I/O对齐，未读真实效果或改变单位/门槛。源码语法及路径检查待最终复核，实际执行及资格审计尚待；不能把冻结计划当来源复现。

## 实际生产与完整审核依赖

源码b7944d3f6eb03057977f5debc289028c5e19309733e7a3423f1b9972d891a66c于2026-10-18 09:30:35.988573 UTC /17:30:35.988573 UTC+8真实启动，09:30:41.455245 UTC /17:30:41.455245 UTC+8完成，终态9c7d900f8abc230b185536859381646e8f8285914504d2fe1b37f72d3b7899e3。全部候选×原情境×两家族资格、每文库失败标记和完整case聚合已生成，真实统计效果尚未计算。

完整输入审核源码1bcd5cc014d2c5268bb8afbda02c62801a68209b3b432d174d8d32fa82aa6790于09:32:48.587184 UTC /17:32:48.587184 UTC+8启动，独立重建文库flag、每case数据/同库关系及全部资格和完整分母；design/zero guards复用原已接受kernel，不称新随机流或P验收。其正式接受尚待后续核验。

## 全部新单位与资格接受

完整输入于2026-10-18 09:32:57.837581 UTC /17:32:57.837581 UTC+8通过，证明1719ebcf16d4634bb803620c048e028605e1ab076356a255d694bc56b22dd46a。全部1040计划、2408 case行、6760文库flag和mixed-block排除逐一重建，最大聚合误差7.11e-15；原rank kernel用于设计/zero guard，不声称新推断随机流验证。两个冻结家族的资格分母在读取新外部效果前接受，后续正式全家族见[推断run](20261018_NSCLC_candidate_external_associations_v1.md)。
