# Codex任务：最新任务状态与进展重新梳理

- 用户要求：再次阅读最新任务状态，梳理相对上一轮解释的新进展。
- 取证：2026-10-12T13:27:02 UTC / 2026-10-12 21:27:02 UTC+8；最新科学发布节点为21:17:00 UTC+8。
- 开始HEAD：e4dc7ec517f084eb932b474b3366cd5861e06a5e。唯一既有未跟踪文件为旧研究审阅日志，保留。暂存区为空。
- 范围：只读状态、既有run、公开汇总与终态；不启动统计、下载、重注释或数据处理。核对“CM识别后”和“免疫治疗单细胞数据”最新聊天消息，与文件证据交叉检查，不向其他聊天发消息。

## 方法和检查入口

读取AGENTS、README、STATUS、ROADMAP、TODO、methods索引、RUNTIME_SYNC、10月12日日志及持续维护的Codex日志。以最终RESULTS/终态/audit优先于producer初始pending和旧TODO状态。

重点证据：

- [新公共数据最终Zarr验收](../../../inventories/runtime/2026-10-12/FINAL_ZARR_ACCEPTANCE_v1.md)
- [原CM全资源提取](../../../inventories/analysis/2026-10-12/CM_edge_LR_extraction_v2/RESULTS.md)
- [CM统一20000次置换](../../../inventories/analysis/2026-10-12/CM_LR_permutation_precision_v2/RESULTS.md)
- [固定候选靶表达](../../../inventories/analysis/2026-10-12/target_association_v1/RESULTS.md)及[ADAM12端点](../../../inventories/analysis/2026-10-12/ADAM12_native_endpoint_probe_v1/RESULTS.md)
- [大类表达](../../../inventories/analysis/2026-10-12/major_pair_native_assembly_v1/RESULTS.md)、[TF输入](../../../inventories/analysis/2026-10-12/major_receiver_TF_preflight_v3/RESULTS.md)、[全量筛查](../../../inventories/analysis/2026-10-12/major_LR_receiver_TF_screen_v1/RESULTS.md)、[联合wild复核](../../../inventories/analysis/2026-10-12/major_LR_receiver_TF_wild_v1/RESULTS.md)
- [巨噬—T初轮](../../../inventories/analysis/2026-10-12/macrophage_T_multisource_response_v1/RESULTS.md)、[共同程序条件结果](../../../inventories/analysis/2026-10-12/macrophage_T_state_sensitivity_v1/RESULTS.md)、对应全部27744汇总记录/终态/审计
- [来源与治疗适用性](../../../inventories/analysis/2026-10-12/source_cancer_context_review_v1/RESULTS.md)
- [空间邻近run](../../../runs/20261011_spatial_conditional_neighborhood_v1.md)、[空间接收程序run](../../../runs/20261011_spatial_receiver_programme_exposure_v1.md)
- [完整基因v17](../../../inventories/analysis/2026-10-12/fullgene_readiness_v17/public_checkpoint_v17.json)、[患者证据覆盖审计](../../../inventories/analysis/2026-10-12/patient_block_provenance_overlay_v1/independent_audit.json)

使用git status/log、rg定位、Python标准库JSON和gzip/TSV选择读取与计数、SHA256，未重算科学模型。聊天读取工具即使includeOutputs=false仍含很长历史，后续在工具内只选最新助手消息，避免继续输出原始历史。一次不存在的public_summary文件和相对helper路径已改为实际audit/工作区路径，无源文件影响。

## 已观察的新增进展

1. 新数据技术合集：9来源、792文库、3180973细胞×49985特征，43.84GB；全部逐行审计及标准读回通过。951701采用污染校正，108635预测doublet额外排除。精确指纹0命中不能消除研究/文库别名重复。来源/作者计数槽/临床正式验收待，source_final_AC及downstream_release_eligible false；原权威库未替换，新注释/整合/CoVarNet/LIANA未执行。
2. 原336候选的完整基因v17准备通过：3来源、164文库、285360选细胞、1121原节点单元、139558612聚合基因值检查；不是新增通讯发现。1290样本患者证据覆盖层新增审阅汇总402文库/175来源内患者块，未改原块，不宣称跨GSE独立患者已核清。
3. 靶表达全1860内部模型已审；ADAM12–ITGB1来源内M1为唯一同时满足诊断BH/稳定区间案例，完整共享程序M2不稳定；原生发送ADAM12稀疏、ITGA9/ITGB1无同细胞检出，机制优先级降低。
4. 原CM68方向/332文库/8GSE全资源提取完成；26994关联记录统一延长到20000置换，主178 q通过、147区间支持（34正/113负）、119方向LR接口。重新计数来源为CRC145、乳腺28、胃3、肺1、GBM1。无同一明确方向/LR/CM在两个GSE同时q及区间通过；是来源内模块关联，不是驱动或独立验证。
5. 空间45接口29可估计/16不可测、主BH0；局部接收程序70模型22可估计/48不可估计，名义p<.05及主BH均0。相对上次“尚未做空间检验”已实质推进，但未得到验证支持。
6. 大类输入/TF覆盖扩至11源、1201793细胞、318有向材料情境；59131 LR单位×3模型=177393全量筛查BH/BY0。5477实际联合wild复核仍三模型BH/BY0，不能把2023名义记录算阳性轴。
7. 巨噬—T pooled新探索覆盖5GSE/5癌种、212库/140块（120已审、20样本替代），27744模型。初轮主BH8/BY3、亚型2/1、全基因11/4；加入发送HSF1及接收同配体RNA后主0/0、亚型0/0、全基因3/0。最新生产13:09:42 UTC和audit13:09:52 UTC均exit0；13:10:29终态原进程退出。24均衡fixture8580组件及选定阳性全TF补充审计通过，不称所有系数完全独立重演。
8. IFNG→IFNGR1_IFNGR2→IRF1主/全基因β.575/.605；五来源IRF1方向正、逐来源删除方向保持。接口联合q主.239260、全基因.049718、亚型1，BY全基因.482737。亚型变体lead TF为HOXA9，不能把lead变化误说同IRF1方向反转。全基因另两条为ICAM1–IL2RG、TGM2–SDC4，三条仅敏感性边界BH不算稳健验证。
9. 新ccRCC GSE29410925库/15来源链接患者标识，response/timepoint缺失，当前对象作者标签缺失/来源待审，不能用于现成耐药预测。

## 解释、下一步与未验证

解释：主要进步为底座可用、分母与患者证据更完整、原CM来源内支持及对共同程序的区分。初轮8结果衰减不能证明解离伪影或排除真实胞外作用；两新增控制尚未拆解贡献。IFNG值得固定独立响应读出追踪，主/亚型校正与真正外部/蛋白/空间验证仍缺。旧DEFB1及336记录保留，不被新候选自动替换。

当前计划：固定独立IFNγ/抗原呈递程序和其他刺激对照、检查IRF1共享基因、拆解HSF1/自身RNA控制、回映原9CM亚型贡献、IFNGR1/2同细胞与空间支持、完成新队列来源/时点/标签/临床验收后再冻结留出验证。当前已核两聊天均idle，最新已记录生产审计终止；未直接枚举服务器所有进程，不能声明全服务器无任务。

未验证：独立新来源通讯轴、受体特异因果、正式状态、CM驱动、免疫治疗预测；新318万未直接并入上述既有缓存的科学结果。数值审计和独立生物验证严格分开。

重要源文件SHA：共同程序RESULTS e71e2e49e38db009103187fdde6f63ff7e9978a4e00f3ec36be0f7ba9e201d28；CM精度RESULTS 91514f40ca8ce4d5d389fa11e26bd0dad6c5d573afec353ea4e1bf20366fba65。

修改仅本任务日志；不重写STATUS或科学结果，不暂存旧草稿。验证所有新增相对链接/脱敏/计数并运行仓库检查，使用git commit --only明确路径避免共享index竞争，按授权普通push及独立SHA核验。实际结果后补，不预称发布成功。

建议commit：`records: reconcile latest CoVarNet and data preparation progress`。

提交前实际检查：仓库验证exit0、0 errors/1既有路径warning；本日志相对链接及公开路径检查通过。IFNG三变体联合q和178来源计数来自现有全量汇总重新读取，不是新的模型估计。
