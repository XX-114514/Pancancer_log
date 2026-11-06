# 2026-11-06 V4注释CoVarNet原文CM发现

用户明确要求在当前修正注释上先按CoVarNet原文识别CM。[方法及边界](../methods/CoVarNet-V7-selective-author-CM-20261106.md)。输入与历史CM、新TME注释分开。

主运行逻辑位置`${PROJECT_ROOT}/runs/20261106_v7_selective_CoVarNet_paper_v2`；保留v1作来源范围敏感性。v1输入阶段实际2026-11-06 12:31:55 UTC/20:31:55 UTC+8开始；v2实际12:36:17/20:36:17开始。来源范围暂存依据GEO公开矩阵范围，在查看CM效应前冻结，不以CM方向选择排除。

v2输入560,871非上皮细胞、164样本、13来源、45节点，12:37:32 UTC输入独立审计通过。K2–20各30次全部成功；12:40:49 UTC网络生产终止，12:40:57 UTC最终独立R/Python核验通过。按严格CCC峰规则选K14；990对相关，18条全局网络边，10个跨父类网络。结果不等于14个独立可靠的多细胞生态系统。

公开[完整结果](../inventories/analysis/2026-11-06/v7_selective_CoVarNet_author_v1/README.md)保留全部因子、rank调查、top15及网络过滤后成员、边、拓扑/置换P和审计，并完成15张图。11个非空网络中1个只含髓系，CM01/02/05为空；CM10/12/14含近零W，CM06/12/14最终成员相同。诊断阈值只标记不筛除。v1诊断12:37:48完成生产，12:40:42全部独立核验通过，选K19。两个范围不进行编号直接互换，未过来源稳健性。

调用合同：研究环境Python通过`runpy.run_path`向输入脚本注入`PROJECT_ROOT=Path(os.environ['PROJECT_ROOT'])`；向Python网络审计注入`RUN_ROOT=Path(os.environ['RUN_ROOT'])`。将对应配置保存为`${RUN_ROOT}/config/paper_config.json`，固定作者R源放`${RUN_ROOT}/author/CoVarNet.R`。输入及审计成功后逐K执行：

```bash
Rscript --vanilla run_author_rank.R "${RUN_ROOT}" 2
# 对预先登记的K=2..20逐一执行；每K内部30次。不要重复启动已完成模型。
Rscript --vanilla assemble_author_modules.R "${RUN_ROOT}"
Rscript --vanilla audit_author_models.R "${RUN_ROOT}"
# 随后执行注入RUN_ROOT的audit_author_network.py。
```

每K单线程、8个K并行；有限controller独立记录起止及退出码，未创建长期自动任务。NMF0.28与作者测试版本差异保留。服务器locale警告未阻止模型/审计，原日志留私有。代码/配置SHA随结果保存。公开仓库不包含私有输入，Cloud无法假定具备真实执行条件。
