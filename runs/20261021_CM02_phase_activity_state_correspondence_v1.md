# CM02分阶段活动与实测状态对应

## 问题及历史纠正

保留2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）。完整固定网络状态分解已接受，但原H的独立单位内平均包含不同阶段，不能直接当治疗前活动。

上一节点将“阶段投影”全部列为未来工作，表述过宽。[20261016原run](20261016_HNSCC_new_CM_baseline_clinical_v1.md)已经完成四scope原minmax/NMFns核对、真正基线投影及完整临床家族，原结论保持。当前不重复临床检验；仅补当前冻结marker状态与原模块各阶段活动的来源内描述。

实际原discovery_v3启动manifest使用downstream recovery v2脚本，SHA为4cf9af42279924f071ec1cffdbf49f1467cd4693aa85732bad6de9ad08485feb；此前引用的material v2通用脚本不是该产物的直接producer。真实代码的独立单位内平均规则一致。默认R缺原NMF；使用manifest指向的scrna_r/R4.3.3、NMF0.28。

## 效果前合同

- 既有原模型、冻结导出、输入频率/元数据、已审作者材料键只读；原W/S/H、网络及权威标签保持。
- 实际RDS的W/S/H/fitted、原freq_normalize、发现训练单位顺序与既有冻结导出逐项复核。有效A=W S满列秩且condition≤1e8。
- 治疗前/普通治疗后各患者先平均归一化材料组成，再对固定A作NNLS；复发仅活动/覆盖描述，不纳入两个阶段的状态对应。NNLS全量BVLS/KKT另核，原拟合F重演H，基线分数与20261016全部严格回放。
- L2 NNLS是冻结组件坐标上的描述性投影，不是重新求解原nsNMF的KL目标，也不要求在真实观测组成上逐值重现原H。原发现缩放包含治疗后数据，不能当外部训练或预测验证。
- 已冻结双调节parent定义、六个观察marker状态、全来源及原七成员并集、pre/post，共48描述位置。模块活动始终用原完整模型，不因网络仅七成员而截断basis。
- 本parent至少20观察且当前T注释至少20才计该病例；相关描述至少20已审病例。报告状态/本parent、状态/当前T、parent/当前T三口径，后两者与主要描述用同一病例。marker parent和原T注释并不严格一致，不能把这些操作性比率当真实细胞比例。
- 2000病例bootstrap，逐次重排平均秩；至少1800有效给点式区间。全部随机流用另一manual-tie rank实现核对，留一仍需至少20剩余病例。不足/常量保留，无新P/BH或按区间挑选晋级。
- 状态与节点/模块共享RNA和发现资料，此为内部对应和组成解释，不能充当独立细状态验证、直接互作或临床增量。

## 实际运行

```text
python -B scripts/describe_CM02_phase_activity_states_v1.py
```

环境由PANCANCER_PREFLIGHT_ROOT、PANCANCER_RAM_OUTPUT、PANCANCER_FINE_LOADER/SCOPE_OUTPUT、PANCANCER_PARENT_OUTPUT、PANCANCER_NETWORK_MEMBER_FILE/SHA256提供。复用已审核原作者桥、原生CSR与marker独立计数/双parent逐细胞scalar函数。新材料/患者/细胞/表达/学习数组仅RAM，写出只含来源聚合和技术收据；没有读取新的临床结果或序列。

2026-10-21T02:57:45.075207+00:00（UTC+8 10:57:45.075207）实际派发，精确PID/startticks及控制目录留私有。producer SHA为7f85348da0a6a2b5600e5d8f6bc01cdb1cfabbb60af396a9530bd374407ed11f。02:58前观察同作业仍运行，实际原模型R核对已通过；完整科学结果尚待终态，不能将本合同当完成。

## 发布与下一步

技术执行、数值复核、完整来源科学审阅分别验收。实际产物SHA与真实事件时间终态后追加；科学效果、数量及图仅授权私有，不提交公共Git。原自动拒绝的公开push与新随访表无特定许可不重试。NeoTRIP已有许可问题未答，未访问受限文件。

下一决定依据完整48位置及来源适配，不只看有利区间：保留CM02的多状态图谱解释；独立同癌种、同阶段、同病理终点和联合空间测量需另做资格，现有空间阴性及原基线临床家族不改。


## 2026-10-21T03:03:02.448435+00:00 UTC /2026-10-21T11:03:02.448435+08:00 UTC+8 分阶段对应真实完成与完整私有审阅

新作业02:57:45.075207 UTC启动，03:01:06.560518 UTC真实结束；转存时原PID已不存在。原模型/训练单位与真实producer、原冻结minmax/A=W S、全部原F→H与20261016真正基线逐值回放通过，新阶段NNLS/BVLS/KKT一致。作者Pre/Post/Recurrence与CM pre/post/post2及审阅患者单位双向唯一另核；该证明没有读取新的疗效或序列。

全原生RNA CSR/bincount与双parent逐细胞scalar、病例阶段状态计数的独立字典实现通过；全部48来源/阶段/范围/parent/状态位置及全部2000 bootstrap手工并列秩重放、留一/分母边界接受。全部来源科学结果与反例私有审阅，改变后续状态解释；不将点式区间当多重检验阳性或临床晋级。

输出资格v2发现原活动profile对小阶段也给分布；保留原私有历史，独立新资格版本将不足5病例的活动/残差分布置空，48主要描述与全部效果未变。公开仅技术证明，源科学值、图和学习参数不公开。原“全部阶段投影待定”已纠正：旧基线分析早已执行，本轮只补当前marker对应。现有空间/原临床完整家族结论均保持。

本轮公开网页原文/GEO核查遇recaptcha/出版入口不可访问，未绕过，也未将新外队列/同终点或受限文件资格写成接受。NeoTRIP许可仍待人类回复。整体goal active，暂无满足明确状态组合—组织位置—同终点增量的已确认互作；下一优先外部联合资料资格及效果前固定的新状态问题。具体科学分解见私有review，技术收据见inventories/analysis/2026-10-21/CM02_phase_activity_state_v1/technical_complete_acceptance.json。原公开push与新随访表被拒动作无特定新许可，不重试；普通本地commit在验证后执行，云端待同步。

[技术验收](../inventories/analysis/2026-10-21/CM02_phase_activity_state_v1/technical_complete_acceptance.json)。


## 2026-10-21T03:04:05.618660+00:00 UTC /2026-10-21T11:04:05.618660+08:00 UTC+8 本轮验证收束

完整仓库validate_repository退出0、0错误/1既有路径警告；原作业已结束、全部来源转存SHA、两新源码与实际运行/资格收据SHA、Python3.10 AST、两inventory列宽及diff检查通过。实际原模型/R/阶段键、冻结基线回放、阶段数值和完整私有科学审阅分别验收；本轮归类progress，不宣称完整互作或临床增量目标完成。

仅显式本轮13文件普通本地commit，保留原三份无关未跟踪文件。原自动审批拒绝公开push仍无特定新许可，不重试，云端待同步。验证日志SHA 85c14096012e60bf04a8e0ef69943ae93c44fecf7110c4710a5b10c9adaee92b；建议commit: Characterize frozen CM02 by treatment phase and observed T states。
