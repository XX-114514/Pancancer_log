> **2026-11-09 完整终态已核实**：[53/53模型/评分/独立审核、53/53 niche及全14CM结果](../../2026-11-09/v4_CM_validation_review_v1/README.md)均已完成；下文运行中数字为历史检查点，不再代表当前状态。

# V4 CM的Visium原文方法检验：生产中

**同V4参考已完成250轮回归/1,000后验并通过独立审核；两张真实空间模型正在GPU拟合，53张全范围及最终CM空间结论尚未完成。** 本页是带时间的检查点，不是成功验证报告。对象是2026-11-06 V4 CM01–CM14，优先解释CM03，其余完整对照，原注释/CM定义不改。

最新固定快照为[01:46:11 UTC／Melbourne 12:46:11](runtime_checkpoint_v6.json)：两张主模型仍拟合，0完成/0失败；**8张CellCharter niche已通过独立审核**，来源级技术汇总见[新表](niches/early_independent_audit_by_source_v2.tsv)。这仍不是任何CM的最终空间支持结论。[本里程碑公开同步收据](publication_checkpoint_v2.json)绑定已核验内容提交；服务器有限队列继续运行，本页不是实时监控。

[实际运行状态](runtime_checkpoint_v1.json) · [参考独立核验](audit/reference_input_audit.json) · [原生输入范围](inputs/source_scope_summary.tsv) · [50标签参考覆盖](reference/node_reference_coverage.tsv) · [来源实测资格](reference/source_measurement_qualification.tsv) · [方法](../../../../methods/V4-Visium-CoVarNet-author-validation-20261107.md) · [run](../../../../runs/20261107_v4_Visium_author_validation_v1.md) · [代码](../../../../scripts/v4_visium_author_validation_20261107/)。

最新取证为[01:24:56 UTC快照](runtime_checkpoint_v2.json)：参考回归仍在真实GPU生产，空间模型0完成/0失败。[13份部署代码独立SHA读回](audit/server_code_binding.json)与公开代码逐字节一致，包含尚未执行、仅在最终真实审核后生成报告的有限接续器；原检查点保留。

[公开同步收据](publication_checkpoint_v1.json)绑定已独立核验的内容提交；同步成功不表示空间检验已完成。

## 实际输入范围

较晚里程碑：[01:31参考完成/两张GPU拟合](runtime_checkpoint_v3.json) · [01:37双分支状态](runtime_checkpoint_v4.json) · [拟合参考独立审核](audit/fitted_reference_audit.json) · [全部50标签×3面板留出诊断](reference/all_panel_label_heldout_separability.tsv) · [完整聚合混淆表](reference/all_panel_label_confusion.tsv) · [来源—面板对应](reference/native_gene_panel_source_scope.tsv)。

最新[01:43:33 UTC检查点](runtime_checkpoint_v5.json)：主模型两张仍真实拟合、0完成/0失败；niche原分支4张＋合格恢复1张均通过独立审核。恢复切片的1,863个原spot中仅3个在HVG内零信息，对niche不评估，1,860个完成；原轴/主模型完整保留。见[来源级早期审核汇总](niches/early_independent_audit_by_source.tsv)和[26份部署代码SHA](audit/server_code_binding_v2.json)。这些是技术/表达—空间建模结果，尚无CM最终空间支持结论。

| CM03成员 | 三种面板留出F1范围 | 主要错误去向 |
|---|---:|---|
| Bcell_Naive_IGHD | 0.811–0.843 | TCL1A naive B、MS4A1/CD79A B |
| Endothelial_Blood_PECAM1 | 0.787–0.818 | ACKR1静脉内皮等 |
| NK_XCL1_CD56bright_XCL1 | 0.637–0.670 | FCGR3A NK、KLRD1/GNLY NK、细胞毒T等 |

每标签有200个留出细胞；这是继承标签的同标签重放，未按来源/患者隔离，不是独立亚型真值或spot去卷积准确率。命名基因在三面板中可测，仍不排除邻近类型混淆。最终须完整保留这些诊断。

CellCharter0.3.7专用环境导入/API与依赖实查通过，独立CPU队列已启动，按切片原生表达scVI→6近邻/3层→AutoK执行，不输入CM标签、丰度、活动或结局。参数在真实CM空间效应前冻结，见[配置](niches/configuration_pre_effect.json)；具体HVG/epoch/K选择是本轮明确实现，不冒称CoVarNet未报告的精确参数。S002因所选HVG内存在零计数spot在训练前失败；有限恢复只对该类失败建立有逐spot记录的niche合格视图，低信息spot仍在原轴和主模型，成功拟合不重复。53请求分母和失败均保留，niche最终以新接受清单和独立审核为准。

| 来源 | 本轮准入切片 | 现有spots | 说明 |
|---|---:|---:|---|
| GSE203612 | 10 | 19,934 | 另2个PDX不作人类原生TME；catalog更大队列不冒充已取得数据 |
| GSE175540 | 24 | 74,220 | FFPE与frozen技术背景保留，不能把24张视作24个独立患者 |
| GSE273952 | 19 | 71,398 | 已有计数/坐标接受凭据；旧去卷积不改名复用 |
| **合计** | **53** | **165,552** | 后续in-tissue/正RNA/坐标规则的实际计数按切片记录 |

GSE238264/GSE301720在本轮已知canonical QC对象位置尚未取得，明确暂不准入，不因此声明全服务器不存在其他档案。legacy-ST、PDX和Xenium不混作原生人类Visium。空间source accession不在V4单细胞来源中，患者跨研究重叠与技术重复并未重新完整认证。

参考保持原50标签，45,385训练细胞＋9,000不重叠留出，41个关键身份基因完整实测来源、13,548共同非MT基因。54来源直接求交集会失去IGHD/PECAM1，故按空间效应前冻结的测量规则选表达参考；这不把被原CM发现排除的来源重新纳入频率矩阵。45个原CM背景节点不变，5个上皮标签仅作混合spot背景。

## 生产与审核状态

- 2026-11-07 01:10:31 UTC完成参考/输入准备，01:10:50独立核验通过：同50标签、行/来源、全部参考非负整数与500原始CSR行分层读回一致。
- 01:16:35 UTC（Melbourne 12:16:35）快照为1个活跃真实参考拟合进程、0个空间模型完成，GPU占用约2.6GiB；仅表示真实任务在跑，不是训练收敛或空间支持证据。
- 有期限队列接续参考250epoch → 每张cell2location 30,000epoch、N=5、alpha=20、1,000posterior（2并行）→ 原文共定位/Moran → 每切片独立公式/权重/epoch审核 → 全14/全来源汇总与独立读回。失败或不可估计不删分母。
- 首次参考训练前的devices重复关键字错误已保留；当前v2脚本仅移除冲突参数，没有重跑已接受输入或改变科学参数。

主posterior为expected mean abundance，q05完整作为敏感性。CM活动用原网络成员和W；主指标直接取固定作者函数，**60近邻、有向、inverse-square距离、行标准W**，不是旧CT75的6邻居版本。背景固定原45非上皮节点，空CM/近零载荷/断开图原样保留。

[明确合成的数学测试](audit/spatial_math_fixture_audit.json)已通过：100个合成spot，45个真实标签名，全部22个非空CM主分数与独立公式误差≤4.44e-16。**该测试不是任何真实Visium切片的验证结果。**

## 结果解释的预设边界

同父类随机成员对照、参考谱混淆控制、各基因面板留出同标签重放将与原文分数分开。累计≤对照比例包含并列，不当作P值或严格超过95%的成功率。若类型参考不可区分，即使空间原始分数高也不能认定精确亚型共处。

不把spot当单细胞接触、不按spot扩充患者做检验、不从共变推出通讯机制。原文CellCharter正交niche步骤已经启动，尚未全范围验收；两项主指标或niche拟合完成也不自动等于生物学验证成功。模型、矩阵、逐spot/细胞/患者表及原始日志私有，公开限汇总、代码和脱敏收据。
