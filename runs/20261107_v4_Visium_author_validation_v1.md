# V4 CM的现有Visium原文方法检验

[方法](../methods/V4-Visium-CoVarNet-author-validation-20261107.md) · [当前检查点](../inventories/analysis/2026-11-07/v4_Visium_author_validation_v1/README.md) · [科学脚本](../scripts/v4_visium_author_validation_20261107/)。对象是2026-11-06 V4的全部14CM，不是旧CT75或论文中同编号CM。

初次取证2026-11-06 22:29:32 UTC/2026-11-07 06:29:32 UTC+8；用户继续后00:54 UTC接续。实际逻辑run`${PROJECT_ROOT}/runs/20261107_v4_Visium_author_validation_v2`，v1准备中因numpy进度计数序列化失败保留。v2准备2026-11-07 01:04:54—01:10:31 UTC完成，01:10:50 UTC独立输入审核通过：50标签、45,385参考＋9,000留出、13,548共同实测基因/41来源，随机500个原始CSR行逐值复核，全部参考计数有限/非负/整数。

准入53切片/3来源/165,552现有spots：GSE175540 24、GSE203612 10、GSE273952 19。PDX2个排除，两个canonical QC文件未取得者单列；不代表所有目录内档案均已审计。细胞参考与空间source accession不同，患者跨研究重叠或技术重复未重新认证。输出公开只含聚合计数、方法、通用代码及脱敏收据。

GPU环境真实导入和CUDA矩阵乘法已通过，旧CPU环境未使用。首次参考拟合01:10:53—01:11:17 UTC在训练前因installed scvi wrapper重复`devices`关键字失败；移除多余参数的新`cell2location_V4_v2.py`复用已接受参考，旧失败脚本/退出码保留。01:13:45 UTC取证有1个真实参考拟合进程、0个空间模型完成；不能把队列注册当作拟合结束。较晚检查点优先。

有限模型controller接续250epoch参考→53个full-batch、30,000epoch/N5/alpha20/1000posterior空间模型（2并行）；有限评分controller接续留出同标签重放→固定作者函数→每切片独立审核；最终controller按全请求家族汇总并独立复核，失败/不可估计不删分母。每段记录时间、退出码与输入版本。不是定时自动GitHub推送；模型和原始日志私有，公开同步仍经审阅。

作者60近邻有向inverse-square Moran与独立矩阵公式在合成fixture验证通过（22个非空CM分数，误差≤4.44e-16）。合成测试初次缺Path导入，仅测试器在R执行前失败，新版本恢复，未影响真实模型。

复现：Python准备/输入审核/资格/汇总脚本通过`runpy.run_path`注入`PROJECT_ROOT: Path`；`cell2location_V4_v2.py RUN_ROOT reference|section [SECTION_CODE]`；`score_author_Visium.R RUN_ROOT BASE_CM_RUN SECTION_CODE`；`audit_score_Visium.py`同样参数。`BASE_CM_RUN`为原主v2。作者R从固定commit获取并放入`${RUN_ROOT}/author/CoVarNet.R`。不要把失败的原devices版本当已接受生产版本；实际阶段是否完成以较晚真实终态为准。

当前为生产中，尚无真实全范围空间结论，亦未完成CellCharter、独立身份、病例层、接触或机制验证。
