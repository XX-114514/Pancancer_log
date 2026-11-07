# 当前V4 CM的Visium原文方法检验启动

用户要求用服务器现有Visium按CoVarNet原文检验全部当前CM。初次2026-11-06 22:29:32 UTC开始取证，用户继续后2026-11-07 00:54 UTC接续；逻辑run为`20261107_v4_Visium_author_validation_v2`，不改原V4标签或CM。

41个身份核心基因完整实测来源、13,548共同非MT基因、50标签、45,385训练＋9,000留出，01:10:31 UTC准备完成、01:10:50独立输入核验通过。原53切片来自GSE203612 10、GSE175540 24、GSE273952 19，共165,552现有spot；PDX及两处未取得的canonical QC对象单列。GPU参考回归01:16:35 UTC仍有真实进程，空间模型0完成，不能提前说CM得到空间支持。

每张后续30,000epoch/N5/alpha20/1000posterior，原作者60近邻有向inverse-square Moran与Spearman主分数、same-label测量/谱系与参考谱诊断、独立审核及全14家族汇总均注册为有期限接续。合成100spot公式测试通过，不冒充真实样本。CellCharter、独立身份/病例/机制尚未完成。

[检查点](../../../inventories/analysis/2026-11-07/v4_Visium_author_validation_v1/README.md) · [run](../../../runs/20261107_v4_Visium_author_validation_v1.md) · [任务日志](../../codex/2026/2026-11-07_v4-visium-author-validation.md)。失败及兼容修正保留；矩阵和个体行私有，公开汇总与通用代码。
公开检查点内容提交`7cb1cab2e2f6a21bc62b76ea9c514cbd6199bb4c`已普通push，独立ls-remote精确核验，见[收据](../../../inventories/analysis/2026-11-07/v4_Visium_author_validation_v1/publication_checkpoint_v1.json)。仓库0错误/1项既有历史警告，新30文件扫描与13份部署代码SHA核对通过。01:28 UTC参考仍实际训练，真实空间结论未产生。

后续01:29:45 UTC／09:29:45 UTC+8参考250轮/1000posterior完成，01:33:10独立读回通过；完整50标签×3面板留出诊断交付，CM03 NK F1约0.64–0.67且与其他NK/T混淆，不能当独立身份真值。01:42:08 UTC主模型两张仍拟合、0完成/0失败；CellCharter独立CPU分支三张完成且独立审核通过，另保留一张HVG零信息导致的训练前失败及限定恢复。26份部署代码SHA一致；主模型/原CM/标签不改，全范围空间和机制结论仍待。

01:46:11 UTC较晚固定快照：主模型两张继续，niche已8张独立通过；低信息恢复仅3个spot对niche不可评估，主原轴不改。参考完成/niche里程碑内容提交`3981a3fb94d09cab24a21f88ca5c797b6e5c3b04`已普通push并独立核验，有限生产/评分/报告队列继续；不是全范围CM空间支持或机制结论。
