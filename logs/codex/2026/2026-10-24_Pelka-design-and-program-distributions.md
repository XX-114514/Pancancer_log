# Pelka设计、现有程序读出审计与后续路线

## 用户要求与起点

继续仔细研究Pelka与其他研究，考虑重新处理或在现有进展上深入调整，组织后续工作。上一目标轮有真实来源权重敏感性完成及验收，是progress；本轮不是等待旧句柄。开始核实际git HEAD、工作区和最新run，确认该实验已完成，没有重复拟合。

## 检查与取证

读取AGENTS、README/STATUS/ROADMAP/TODO、RUNTIME_SYNC及论文叙事规则，最新来源权重/CRC DIALOGUE/状态DIALOGUE的run与methods，以及当日daily。保留全部并行atlas代码、状态和未跟踪文件。

源码审计确认两个现有DIALOGUE投影入口以患者/分析块均值进入后续检验。现有均值家族的完成不能充当分布尾部家族的检验。本轮未读取任何新原生表达或个体记录。

通过原始文献正常搜索、作者工具文档和GitHub连接器阅读Pelka代码及DIALOGUE、mc-ASTRA资料，核muscat与MuVIcell、NSCLC空间研究的公开摘要。PMC入口CAPTCHA及Nature认证跳转后停止，未取受限正文。Pelka网络代码两个默认tau不同且允许覆盖，实际论文调用未核，记录为精确复现的未决项。

## 本次实际工作

新增通用人工读出模拟、两份完整人工汇总及[run](../../../runs/20261024_program_distribution_design_simulation_v1.md)。先720个排列逐列独立核数值及反例，再执行全部场景300重复；因Monte Carlo精度不足，另固定完整3000重复并保留旧轮。两次actual exit0，详细时间、环境、参数与代码SHA在run/JSON；不存在真实人群新P或新模型。

新增[后续路线](../../../methods/Pelka-informed-program-distribution-roadmap-20261024.md)：P0覆盖/字典资格→P1分布测量→P2来源内状态网络→P3组成CM/来源迁移→P4组织/临床/通讯。明确来源依赖的价值、临床引导训练的循环风险、亚型与状态分层、尾部统计量精度、预定完整家族和新探索属性。链接加入STATUS/ROADMAP/TODO。

## 失败与修正

- 工作区根不存在AGENTS.md，仓库AGENTS已正常读取；未影响执行。
- 两次shell通配路径未匹配，只造成检索exit1，随后按实际文件名读取。
- GitHub首次误用master返回404；查看公开目录确认main后，连接器正常读到README及两份代码。
- 本机SciPy不可用；人工实验仅依赖既有NumPy，未安装或修改环境。
- muscat官方vignette超过工具内容长度，未读取全文；仅按原始摘要阐述方法范围，不声称方法细节已全文核验。

## 已知边界与后续

人工零场景部分Wilson区间不含5%，不报告全面校准通过；完整模拟限制在run中列明。尚未在真实数据上证明稀有状态是旧阴性的原因。新路线P0–P4未全部执行；下一先覆盖/测量审计，不立即全库重跑，也不重开旧家族。

本轮没有轮询服务器科学任务、更新上传或更改定时任务；这些状态不随本方法记录刷新。公共push此前自动审批拒绝保持，本轮不重试，云端待同步。并行atlas修改不进入本任务提交。因共享当日日志存在未提交并行修改，本任务使用独立daily节点，保留共享文件。

建议commit：Audit programme summaries and design a Pelka-informed state network study。


## 最终检查

仓库全量validate_repository实际exit0，0 errors/1既有warning。两轮各120位置、代码SHA、重复数和拒绝比例算术独立检查通过；全部新文档相对链接通过。额外Wilson区间包含点估计的严格比较首次遇到100%时浮点上界0.9999999999999998而失败；确认仅2.22e-16舍入，复核使用1e-14算术容差后通过，结果与模拟代码未修改。该数值容差不改变统计门槛或名义零场景偏差结论。

人工diff/路径/敏感内容审查仅纳入本任务通用方法、代码、人工结果和日志；并行atlas文件保持未暂存。公共push旧拒绝保持，未执行。

全量检查日志SHA256：`ae19d66d087e8ad4b547094098602beff0234bedeb894e4cc23dd6715762e67e`。
