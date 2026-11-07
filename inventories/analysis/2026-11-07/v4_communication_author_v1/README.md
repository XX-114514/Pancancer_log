# 当前V4通讯分析与空间队列进度

本页为带时间检查点：**空间主模型4/53已完整拟合、原文评分及独立审核通过，另2张运行；CellCharter全部53张拟合及逐张审核通过。当前V4通讯流程已真实启动输入准备，还没有正式通讯结论。**

[空间实际状态](spatial_runtime_checkpoint_v1.json) · [耗时与预计](spatial_runtime_ETA_v1.json) · [全部53张niche来源汇总](all53_niche_audit_by_source.tsv) · [通讯实际状态/代码SHA](runtime_checkpoint_v1.json) · [冻结配置](method_freeze.json) · [方法](../../../../methods/V4-CoVarNet-communication-plan-20261107.md) · [run](../../../../runs/20261107_V4_CoVarNet_communication_v1.md) · [代码](../../../../scripts/v4_communication_author_20261107/)。

05:19–05:20 UTC／墨尔本16:19–16:20取证，已完成4张均为GSE203612，每张完整30,000轮/1,000posterior，耗时80.39/82.54/84.68/86.79分钟。四张均通过作者共定位/60近邻有向Moran与独立核验。按两并行估算剩余49张名义约34小时；考虑后续来源/spot数差异，规划区间为**35–50小时，墨尔本11月9日凌晨至傍晚**。仅4张同来源时间标定，尚不能承诺精确完成时刻；故障、后续来源速度和资源变化会影响估计，不为压缩时间改变epoch或GPU并发。

CellCharter于03:10:48 UTC完成53张接受清单：GSE20361210、GSE17554024、GSE27395219，独立邻域特征误差≤2.38e-7。仅恢复视图的3个HVG零信息spot对niche不评估；主空间轴/模型不改。最终CM—niche联合解释仍依赖完整主去卷积分数，niche技术完成不自动等于CM生物学验证。

通讯输入严格复用当前V4主v2：13来源、164样本、560,871细胞、45非上皮标签和全部14因子。05:28 UTC实际输入准备已扫描原轴640,000行、抽取124,984个目标细胞，有限控制器真实存活；尚未完成输入审核或技术小试，不能称CellPhoneDB已跑完。[CM比较资格](audit/CM_comparison_eligibility.tsv)和[CMT来源样本数](audit/CMT_source_sample_counts.tsv)在通讯效应前固定。[来源实测范围](audit/native_measurement_scope.tsv)避免缺测当0。

执行顺序：原始输入及独立审核 → 明确非科学的20置换小试 → 每标签≤1,000的全局1,000置换CPDB → 14个CMT High＋11个非空CM pooled-Low共25组正式推断；并行CellChat替代算法及Immune Dictionary细胞因子响应 → 全家族比较/独立审核。High为CMT原H最高，不是中位数分组；正文pooled-Low与代码median-CMT-Low分别保留，缺失对NA。CM03 High19/Low145，部分CMT强来源依赖；细胞对P值不是患者环境效应证据。

通讯时间目前只能作**暂定4–12小时的CPU计算窗口**，从输入准备开始计，未用本轮正式任务实测标定；全量Low组和CellChat可能更慢。小试/全局正式结束后以实际耗时重估，出现软件/基因测量失败则保留真实阻断。空间最终对齐和图件审阅的完成时间另算，不承诺本窗口完成最终论文结果。
