# 完整来源谱系及功能基因资格核查

## 要求、起点和方法

持续推进有生物学或临床价值的细胞协同/互作证据。肺癌同源有限功能背景诊断已经结束；本轮先核全部既有来源的可比较性，不追加旧临床/LR阈值。[方法合同](../methods/complete-external-lineage-and-gene-qualification-v1.md)。原52来源/1290文库/172节点和免疫9来源/400文库/270节点分臂处理；不称61独立队列。

## 实际事件与失败修正

- 2026-10-16 15:50:53 UTC /23:50:53 UTC+8：完整谱系机会v1启动；读取旧TNBC未映射sampling sidecar触发断言，真实退出。15:52:14 UTC依赖审计启动，随后因生产已退出而停止。两者保留日志，未写新科学结果，未修改数据。
- 15:54:37 UTC /23:54:37 UTC+8：v2改为已验收双官方token sampling sidecar v2，并核其终态/文件SHA；范围、细胞门槛、节点定义不变。生产代码SHA 8940b97e57f2a3406a0a9cd0c6950c8cbc0f883836c329cfeb0771ca0e46ec1d，helper 7a207d7453a0493795780f8c040def43638fbd29abc2b8637c9c1e8fa96bced6。
- 15:55:40.016596 UTC /23:55:40.016596 UTC+8：全部来源/文库、大类、细候选、分层和明确配对捕获生产完成。
- 15:56:08.189961 UTC /23:56:08.189961 UTC+8：独立全图谱细胞编码重计数、全部池/来源/分层集合交集及患者配对核对通过。审计代码SHA a3e03ab39ad31c8f45377b07cb0de07bcc92109e78a5713bc904b85a199b6cdc，生产终态SHA 663a0f4d0e481dec196f91dadbf98dd4d170d24878816c40c1470a514bbfcb52。
- 16:02:05.168725 UTC /2026-10-17 00:02:05.168725 UTC+8：全部来源功能/程序物理基因覆盖检查真实启动，代码SHA ef3b8323f984e9b3ea00e06f7c78e2829e79f81230749c8ef83c7be52cb5005e，输入helper dfe121fe1b09aa3785c23fa6e0a9abeb990807ef1d64e4ce2965e0b342f72766。最终结果及独立覆盖审计待。

16:01—16:03 UTC再次取证原收据和真实/proc句柄，确认谱系生产与审计均已终止且PASS文件存在。之前工具输出截断只代表观察丢失，未因此重启。上述早期事件在本轮补录，不能当时已提交的日志。

## 证据、执行和边界

服务器scanpy Python，单线程OMP/BLAS，使用PANCANCER_PROJECT_ROOT/PANCANCER_PREFLIGHT_ROOT。执行preflight_complete_external_lineage_opportunities_v2.py、audit_complete_external_lineage_opportunities_v2.py及preflight_all_external_functional_gene_coverage_v1.py。启动收据绑定源代码SHA和PID/startticks，终态以私有文件而非等待超时判断。

谱系完整证明：${RUN_ROOT}/complete_external_lineage_opportunities_20261016_v2/private_independent_opportunity_audit.json；基因覆盖目录：${RUN_ROOT}/complete_external_functional_gene_coverage_20261016_v1。患者、完整表/科学汇总保留server及获准当前对话；公共Git无患者行、基因权重或新效果。

本轮尚无新表达评分、P/BH、临床模型或直接互作确认；技术资格不能提升为科学验证。下一独立核全部基因覆盖，再连接捕获/来源角色以冻结实际可执行外部问题。无下载、删除或原身份/2026年10月CM（v20261010）变更；完整goal active。

建议commit：analysis: qualify complete external lineage and native gene coverage without effect selection。

## 16:06 UTC /2026-10-17 00:06 UTC+8 完整基因覆盖验收补录

生产16:03:26.573209 UTC结束；独立审计16:05:26.786789 UTC /00:05:26.786789 UTC+8通过全部48800功能、3416程序、800原肺癌校准及61轴。核全部文件SHA、原/剩余参考门槛、原53系数槽和3无系数槽、别名多匹配/物理碰撞及所有失败，未重用生产覆盖函数。审计代码SHA 67e586020fac9c99a3c313d18c98c6325c4557c71cb560c02b0f5a9c77a678fb，证明SHA 0e0407ab6e6992c9be6056abb33944f184a78dcab522c6b904a494a3d78853c1，生产终态SHA 2b89d8172b5f9b9616448fc158933927dc9d73f388775951568c624105bd4a4a。此前pending是当时状态。

下一连接完整原context的准确端点、宽谱系与多基因候选机会，登记每个context真实发现来源。头颈癌优先进入原生功能的治疗前后测量准备，其余来源完整保留；尚未计算新外部效应或重开原头颈癌临床家族。
