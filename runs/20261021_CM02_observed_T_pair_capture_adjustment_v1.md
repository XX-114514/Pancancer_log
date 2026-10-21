# CM02观察T状态组合与RNA捕获背景

## 前一节点与本轮目标

上一轮GSE288199全RNA/96测量范围实际完成，数值及私有来源审阅接受，为progress。临床TableS1用户稍后补充；其病例/药物桥未接受，沿原门槛不在该源计算关系效果。本轮在已有已审病例桥的GSE286827推进，保留2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）优先。

## 效果前冻结合同

- 使用原权威Zarr、已审材料/病例/药物HPV列、原生RNA及原完整七网络成员；不读取新疗效/TCR，不改W/S/H、注释或旧家族。
- 同原六观察recipe，全部3调节×3效应配对；pre为主要描述，post为同队列次要描述。不是独立来源复现。
- 完整来源为主，CM02原七成员并集为敏感范围；严格/扩展调节parent、全部原生保留RNA/原文库内低中高捕获，共2phase×2context×2parent×4depth×9=288位置全部保留。并集不是全部NMFbasis。
- 每患者/阶段先合全部文库；两parent各≥20观察、marker-T≥20，药物/HPV确定且阶段一致。至少20已审病例、完整设计满秩、condition<1e6和至少10残差df；缺失/共线/常量不删协变量或换范围补救。
- 状态/本parent为主。基础partial Spearman调整两parent log计数、药物和HPV；捕获调整再加入两parent各自非recipe人类RNA log1pUMI及log1p检测基因的细胞中位数，共9列控制设计。仅从完整human计数中减原22个recipe基因。这四个代理也可能携带生物激活，调整前后并报，不能称纯技术真值。
- 连续控制变量及六状态在每次拟合重排平均秩、标准化连续控制。病例在药物×HPV内抽样2000次，同一draw用于完整六状态及两个模型；有效≥1800才给点式CI。全部留一仍须20病例。所有主要/敏感性机会均描述，无新P/BH或按区间晋级。
- 同一病例子集另报未调整Spearman、state/marker-T与parent/marker-T相关；改变分母是改变估计对象，不能挑有利结果。基础和捕获模型各自资格及有效抽样并列。
- 每原生cell计数/marker规则沿原独立CSR/scalar核；患者阶段计数用另一字典求和，捕获中位由另一数组索引核。统计SciPy秩+lstsq/QR对照独立Pandas秩+SVD，每个实际draw/留一及敏感性均核，所有不足位置保留。
- 仅输出GSE级聚合、技术证明；全部cell/case/count/quality/设计/随机draw数组仅服务器RAM。不以cell/ROI扩大患者n，不把RNA共变当直接LR或空间证据。

本合同及实际启动脚本先于新效果；精确冻结/派发时刻、源码SHA与终态后续追加。当前尚未执行新效果。

## 实际入口及输入环境

以下为实际调用的等价脱敏入口，精确控制目录/PID和源码SHA保留私有启动证明。使用既有scanpy环境，单线程；本轮没有包安装。

```text
PANCANCER_PREFLIGHT_ROOT=${PROJECT_ROOT}/Project_v5/runs/20261011_GSE243013_preflight_v1
PANCANCER_PARENT_OUTPUT=${ACCEPTED_LIBRARY_BRIDGE}/private_output
PANCANCER_FINE_LOADER=${ACCEPTED_PATIENT_SCOPE}/qualify_CM02_fine_state_patient_scope_v1.py
PANCANCER_FINE_SCOPE_OUTPUT=${ACCEPTED_PATIENT_SCOPE}/private_output
PANCANCER_NETWORK_MEMBER_FILE=${ATLAS_ROOT}/CM/T_whole_TME_primary_tumor/discovery_v3/each_node.tsv
PANCANCER_NETWORK_MEMBER_SHA256=5731d7a3f32f1988358f3482fa103630a386dbdc8b24160e74fa2029f5ae0351
PANCANCER_RAM_OUTPUT=${TASK_ROOT}/private_output
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
${EXISTING_SCANPY_PYTHON} -B ${TASK_CODE}/describe_CM02_observed_T_pair_capture_v1.py
${EXISTING_SCANPY_PYTHON} -B ${TASK_CODE}/audit_CM02_observed_T_pair_capture_source_v1.py ${TASK_ROOT}/private_output ${TASK_ROOT}/source_audit
${EXISTING_SCANPY_PYTHON} -B ${TASK_CODE}/plot_CM02_observed_T_pair_capture_v1.py ${TASK_ROOT}/private_output ${TASK_ROOT}/figures_v1
```

终态技术证明：[acceptance](../inventories/analysis/2026-10-21/CM02_observed_T_pair_capture_v1/technical_complete_acceptance.json)。来源科学结果、完整源分母资格及图只留授权私有。

## 后续依赖

本轮回答观察状态关系能否超出数量及RNA捕获解释。稳定组合仍需独立队列、给定丰度/区室的组织关系及同终点临床增量。原空间/临床阴性家族保持。GSE288199临床缺表deferred；服务器来源结果回传/原公共push/新随访表受拒动作无特定新许可，不重试。


## 观察状态关系实际派发

2026-10-21T04:21:52.671219+00:00 UTC（UTC+8 12:21:52.671219）服务器真实启动；通用源码五份SHA与本地一致，完整2000合成重采样/留一、并列秩、常量和20病例/共线门槛已通过。仅通用源码上行，没有本机GSE288199派生结果上行，不涉及前述受拒回传。新源内分析仍只输出来源汇总，原数据/个体数组均服务器RAM。

随后同一PID/startticks取证仍活跃，真实输入/全统计终态尚待，不能由派发声称效果接受。源码main SHA6df0aea6104cf7885fc1936b87818460e3132735db222d7b44f9d8218e37edb4，统计helper SHAab636d6664c13c82f45f5e947580079f48b9cd6c7992de510bf79ab5fe58d17e。shell既有locale警告不当科学作业失败。独立输出审计将明确原producer所记source_phase_cases是该context有观察的病例数，同时补完整源阶段病例分母；不把零观察漏成不存在的机会。

Recorded/accepted 2026-10-21T04:37:03.790941+00:00 UTC /2026-10-21T12:37:03.790941+08:00 UTC+8


## 原生观察状态关系的整批完成与范围接受

实际生产于2026-10-21T04:30:25.948612+00:00 UTC（UTC+8 12:30:25.948612）结束；同PID已不存在。全部预设288位置、原CSR/marker scalar和患者字典计数/两parent捕获中位独立核对、每个实际bootstrap/留一/敏感性SciPy+lstsq+QR对照Pandas秩+SVD完成。没有低于20已审病例的效应替代，所有不足位置保留，原父类/药物HPV/完整设计门槛未变。

完整输出来源审计及私有森林图实际生成并视觉检查；所有来源/阶段/范围/parent/深度和点式区间资格私有审阅结束。原producer的source_phase_cases仅指该context有观察的病例，审计资格副本另加完整authority阶段分母及无context观察数；原产物/效果保留，不把缺观察填阴性。14份仅来源汇总/技术/图转存与全部实际源码SHA逐项通过，图例、完整九行/两阶段及区间/来源边界可见。未读取新疗效或TCR值，不计算新P/BH或临床/空间家族，不以源内描述充当独立验证。

实际统计环境核为Python3.10.19/NumPy1.26.4/Pandas2.3.3/SciPy1.11.4/Zarr2.18.3/Matplotlib3.10.8；本轮inventory早期环境标签推记Python311，已在提交前依真实环境更正为Python31019，不改历史其他run。未安装/更新软件。具体科学值、数量、候选解释和图仅授权私有；当前公共只技术轨迹、通用代码和SHA。原CM02保留图谱优先，稳定细状态—组织关系—同终点增量目标尚未完成，整体goal active；本轮是progress。

本有限描述轮次收束，下一依完整私有解释推进独立临床/队列资格和同状态测量，不通过变阈值或挑敏感性挽救。GSE288199 TableS1用户稍后补，不重复请求，也不把未审coded来源单位作患者关系检验。来源派生结果回传/公共push/NSCLC新随访表三类受拒动作均无新特定许可，不重试或绕过；本轮上行只有通用代码，原被拒GSE288199派生结果仍仅本地。

Technical acceptance SHA256: 1707309832ee1d5f416b04edc472d9f12bbcb4abd0dec92afdbdb77f98fda5b6
