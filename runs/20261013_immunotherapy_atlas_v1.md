# 免疫治疗专用单细胞对象与下游 v1

- 授权日期：2026-10-13（UTC+8）。
- 当前阶段（2026-10-14 11:49 UTC+8）：RNA/15CM完成；logspace_v2运行，16复用+11新success+1skip，真实每库20—47秒；完整400库与最终Zarr待验收。
- 输入：20261012_public_preprocessed_zarr_v2完整技术合集，Project_v5及历史免疫队列的既有QC/counts与作者临床资料。
- 来源角色：独立验证/内部复现/机制参照/产品分支/隔离逐原研究登记。不得将Spatial、动物、bulk计为单细胞；CAR-T产品和免疫富集单独标记。
- 输出：原项目卷新独立run下计数权威对象、降维/图/聚类、major/fine/state/单marker编号描述、逐样本CNV证据与精确回填、CM新发现及旧9CM固定投影；持久化前需阶段实测容量与完整验收。
- 注释：复用当前P3证据与命名方法，身份/恶性/状态分别保存，未决保留；不以单marker硬给成熟亚群。
- CNV：复用malignant_inference核心，原文库原生counts和同文库参考；skips不当normal；可用时既有CopyKAT证据只单独记录。
- CM：R CoVarNet新版本，rank survey与稳定性不强制K9；患者/材料/来源内共变，未捕获谱系为缺失。旧W固定投影与新W发现分开。
- 新增公共矩阵下载：本任务不恢复，优先现有Project_v5合格数据。新授权允许下游计算及独立版本。
- 详见[持续任务日志](../logs/codex/2026/2026-10-13_immunotherapy-atlas-downstream.md)。启动代码/参数/实际范围和终态按证据追加。

## 已执行节点

- 来源盘点完成：21条来源、5,274,578细胞登记，角色分层详见[汇总](../inventories/runtime/2026-10-13/immune_atlas_v1/source_preflight_public.json)。不等于临床验收或独立验证。
- GSE305078：原生count检查PASS，32文库DD/DecontX后台运行；独立复核pending。
- 下游pilot：2026-10-13 13:47:03 UTC / 21:47:03 UTC+8开始，PID 341235，12个免疫主队列/待审来源，每来源600细胞；冻结代码aef9684。范围之外的机制参照/CAR产品/基线重叠没有混入该模型。

- 生产链实际启动：2026-10-13 14:28:19 UTC / 22:28:19 UTC+8，PID139233，冻结7cc78f9；待额外QC独立验收后依次进入冻结、范围pilot、全量生物学与最终对象。
- CoVarNet R接口合成测试PASS：[报告](../inventories/runtime/2026-10-13/immune_atlas_v1/synthetic_CM_API_public.json)；真实细胞0，未声称真实CM完成。

- 首次驱动在临床原库别名门槛退出；GSE316195用原生GSM精确crosswalk修复，保存物理GM库名。恢复于15:12:51 UTC / 23:12:51 UTC+8启动，PID369104，冻结1eef352，只读复核/复用已审计准备数据。
- GSE305078全量独立QC报告：[验收汇总](../inventories/runtime/2026-10-13/immune_atlas_v1/GSE305078_preprocessing_acceptance_public.json)。预计主对象9来源2,670,057细胞/400原库；正式下游和最终权威仍pending。

- 15:18:14 UTC / 23:18:14 UTC+8队列冻结PASS：[实际范围](../inventories/runtime/2026-10-13/immune_atlas_v1/cohort_freeze_v2_public.json)，9来源2,670,057细胞/400物理库。新范围pilot已启动；完整生物学与最终Zarr尚未完成。

- 正式9来源pilot PASS：5,400细胞/22,351共享features，[实际指标](../inventories/runtime/2026-10-13/immune_atlas_v1/integration_pilot_v3_public.json)。15:20:53 UTC / 23:20:53 UTC+8全量2,670,057细胞生产整合已开始；没有声称全量下游完成。

- 大样本后端线程上限修复已实测：[模拟12,000坐标PyNNDescent检查](../inventories/runtime/2026-10-13/immune_atlas_v1/large_neighbor_backend_API_public.json)。v4生产恢复15:38:26 UTC / 23:38:26 UTC+8，PID74376；冻结QC/pilot复用，无新的全量完成声明。方法详见[合同](../methods/immune-specific-atlas-v1.md)。

- 15:52 UTC / 23:52 UTC+8健康取证：驱动74376/子进程74379仍运行，GSE169246329,728/489,360原生块已处理；生产全量整合未完成。后续阶段由已部署的验收执行链继续，正式权威发布pending。

- 全量真实模型已完成：[报告](../inventories/runtime/2026-10-14/immune_atlas_v1/integration_integrate_public.json)，2,670,057×22,351、12粗聚类。注释loess失败保留；00:41:42 UTC / 08:41:42 UTC+8恢复PID382735，真实[HVG检查](../inventories/runtime/2026-10-14/immune_atlas_v1/HVG_recovery_API_v5_public.json)通过。

- 00:49 UTC / 08:49 UTC+8生产谱系检查点：B16、Endothelial15、Epithelial33、Erythroid13群，397,121细胞；Fibroblast在处理。仅部分注释，正式CNV/CM/最终权威pending，目标active。

- 01:32:46 UTC / 09:32:46 UTC+8全量注释完成：[脱敏汇总](../inventories/runtime/2026-10-14/immune_atlas_v1/annotation_public.json)，2,670,057细胞/270子群；ID唯一实际核验。完整原生marker精炼运行中。免疫富集和T谱系宽范围含多组织，仅探索，追加材料分层。

- 全原生marker：[完成报告](../inventories/runtime/2026-10-14/immune_atlas_v1/marker_refinement_public.json)、[完整270群目录](../inventories/runtime/2026-10-14/immune_atlas_v1/subcluster_catalog_aggregate_public.tsv)、[目录审计](../inventories/runtime/2026-10-14/immune_atlas_v1/annotation_catalog_audit_public.json)。
- [材料实测](../inventories/runtime/2026-10-14/immune_atlas_v1/CM_material_preflight_public.json)通过；宽范围/材料CM主模型后LOSO因不可观测节点失败，新discovery_v2恢复合同保存旧拟合，输入先全细胞重计。逐库CNV保持独立运行。
- [旧P3固定W门槛](../inventories/runtime/2026-10-14/immune_atlas_v1/P3_fixed_CM_projection_gate_public.json)要求多基因映射；本轮新发现与固定9CM验证分开，未宣称后者完成。

- 02:16:29 UTC / 10:16:29 UTC+8：[实际检查点](../inventories/runtime/2026-10-14/immune_atlas_v1/downstream_checkpoint_20261014T0217_public.json)，6范围CM完成、4/400库CNV success。转移分支NMF成功但无支持网络，修复原CoVarNet空表行名后新版本继续，不降低科学门槛。

- [全15范围CM完成](../inventories/runtime/2026-10-14/immune_atlas_v1/CM_recovery_v3_public.json)及[独立接口审计](../inventories/runtime/2026-10-14/immune_atlas_v1/CM_artifact_audit_v1_public_recovered.json)PASS：10范围有支持网络成员、5无支持成员，6准备不足；均非临床独立验证。待原CNV完整覆盖400库通过后自动打包及全行/复制/读回审计。

- 03:07:44 UTC原CNV与待打包驱动控制停止；[输入纠正门槛](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_input_contract_correction_v1_public.json)。[真实读性能](../inventories/runtime/2026-10-14/immune_atlas_v1/native_read_benchmark_v2_public.json)排除raw read作为主要耗时；新版本保证counts保留/独立log层/2进程×4原库worker、科学AST与完整映射门槛。旧图/旧结果保留不发布，主RNA/CM不重复。

- 03:26:08 UTC正确logspace版本PID153444启动；[AST最小纠正实际验收](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_logspace_core_patch_public.json)、[首4真实文库表达/计数/profile](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_logspace_first_real_checks_public.json)PASS。新计算继续全400库，旧raw整数直接推断不复用；最终pack使用明确版本CNV结果。

- [真实坐标匹配性能与全gene等价性](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_gene_lookup_benchmark_v1_public.json)：209.178→0.0139秒，基因集合/计数/推断不改；16正确logspace_v1库保留，v2按精确ID/方法/表达门槛复用。

- 03:47:25 UTC v2恢复PID251775（4worker），03:49:46 UTC[真实加速/复用检查点](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_logspace_v2_first_recovery_checks_public.json)：28/400库有状态，规范表达全条目误差门槛继续PASS，目标active。

- 最终接口审计入口`audit_immune_final_interfaces.py`已准备：从实际data4成品重读节点/CNV/图与15个CM频率，依赖逐行原生计数和逐文件复制校验先PASS。当前尚未正式发布。

- 04:09:28 UTC[全部400库工作进程完成](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_logspace_v2_workers_completed_public.json)：149success/251skip、4×exit0，无method fail；最终精确合并和authority仍pending。

- [全量CNV精确映射](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_logspace_v2_public.json)PASS；[归档类型故障](../inventories/runtime/2026-10-14/immune_atlas_v1/authority_pack_metadata_type_guard_v1_public.json)实际确认，旧[打包终态](../inventories/runtime/2026-10-14/immune_atlas_v1/CNV_logspace_v2_pipeline_terminal_public.json)保留。`run_immune_authority_recovery.py`不重拟合科学模型，只恢复candidate_v2归档/验收。

- 04:16:58 UTC[元数据全行转换](../inventories/runtime/2026-10-14/immune_atlas_v1/authority_metadata_conversion_v1_public.json)PASS；[恢复快照](../inventories/runtime/2026-10-14/immune_atlas_v1/authority_metadata_recovery_v1_start_public.json)实际启动，archive candidate_v2继续。参见[下游交接协议](../docs/IMMUNE_ATLAS_V1_HANDOFF.md)。

- [全基因写入检查点](../inventories/runtime/2026-10-14/immune_atlas_v1/authority_fullgene_pack_progress_v1_public.json)：前5来源完成，超过百万行；这是处理进度，最终发布仍pending。

- [全部9来源native数组写入](../inventories/runtime/2026-10-14/immune_atlas_v1/authority_fullgene_pack_progress_v2_public.json)完成，实际候选大小及最终500GB预留预核fit；完整全行和成品读回审计仍pending。

- 05:03:15 UTC[完整原生计数全行审计](../inventories/runtime/2026-10-14/immune_atlas_v1/independent_immune_fullgene_audit_public.json)PASS，X最大误差0；[正式容量门槛](../inventories/runtime/2026-10-14/immune_atlas_v1/immune_final_storage_budget_public.json)fits。data4 copy与成品接口审计继续。
