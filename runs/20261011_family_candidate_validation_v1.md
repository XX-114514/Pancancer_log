# Run 细胞家族候选与身份内部程序检查

- Run ID：`20261011_family_candidate_validation_v1`。
- 父运行：`${PREFLIGHT_RUN_ROOT}`，对应Project_v5/runs/20261011_GSE243013_preflight_v1。
- 基线：v20261010，1,290样本、172节点、9CM；当前扩展家族分析不改写原节点。
- 状态：本轮计算、输入/结果审计和数值QR复算已完成；空间接触、正式状态、机制与临床验证未完成。
- 当前证据登记：native_LR_program_evidence_register_v11；这不是新的注释发布。
- 执行解释器：scanpy环境Python3.10.19；本次观察NumPy1.26.4、pandas2.3.3、SciPy1.11.4。历史逐作业完整环境和运行前Git commit未全部固定；脚本SHA和开始/结束时间见证据摘要。

## 任务轨迹

| UTC+8时间，2026-10-11 | 节点 | 实际结果与下一步 |
|---|---|---|
| 17:11:19 | 肺癌作者身份/RNA深度模型数值审计结束 | 87条记录、435模型；身份混合与RNA深度影响保留为条件敏感性，继续检查逐细胞来源 |
| 17:28:28 | 原生逐细胞输入审计结束 | 3来源、285,360细胞、296基因的目标panel，缺基因明确记录；继续同细胞复合物阳性检查 |
| 17:35:34 | 身份/原生输入里程碑登记 | 接受来源和癌种限定关系；编号和状态轴分开；651记录审计当时仍运行 |
| 17:44:17 | 651记录稳健性独立复算结束 | 各1,302,000次bootstrap/残差置换、15,589块删除；258条稳定，其中112条参考支持 |
| 17:44:27 | 全336记录原生表达审计结束 | 4条频繁同细胞RNA、332条稀疏；稀疏不代表蛋白互作不存在 |
| 17:47:24 | 修正编号组成共线模型后审计结束 | 全336×5=1,680模型；按实际存在节点选基准，旧共线结果保留 |
| 17:53:04 | 全336记录分类审计结束 | 28条新增条件关联保留；不手选替代端点，不替换全发现BH |
| 17:58:23 | 14文件空间panel审计结束 | 4,704候选×文件检查；只验证基因覆盖和元数据，未读取空间X分析邻近 |
| 18:11:33 | 自动可行性审计结束 | 全336按同一规则，15条后续记录：13接触检查、2程序检查，11种LR组合 |
| 18:14:23 | 15记录身份内部程序QR复算结束 | 3,710逐样本/分区模型；136块级检查，80可估计，56覆盖不足 |
| 18:16:10 | 本轮里程碑和状态刷新 | 0登记冲突；已登记本轮作业均终止，下一步空间原生表达与样本映射预检 |

完整开始/结束时间、时区及脚本SHA见[时间线](../inventories/analysis/2026-10-11/task_timeline_v1.tsv)和[证据摘要](../inventories/analysis/2026-10-11/family_candidate_checkpoint_v1.json)。这是依据已有证据的补录，不宣称这些节点当时已经同步GitHub。

## 方法与选择理由

接受来源/癌种限定的共变；多个癌种分别存在的效应具有探索价值，分别报告癌种内关系与癌种之间均值差异。原始编号保留，来源和表达审阅后建立分析家族；28组同名/103节点的覆盖核查不等于全部身份等价或全部可通讯检验。

以本群原生整数count>0判阳性；复合物要求所有基因在同一细胞阳性。发送者和接收者在同一样本各至少20阳性细胞，至少10继承分析块，作为可行性门槛；接触支不要求发送者程序调整后保留，程序特异支要求保留条件关联及足够R+/R−细胞。

15条记录在每个原始样本分别比较受体阳性/阴性细胞。程序排除LR和接收者命名基因；每组至少20细胞。原编号内部、肺癌作者CD4/CD8/未分配及全家族分别报告；全家族深度模型还加入实际存在原节点的固定效应。RNA深度模型包含logRNA和lognGenes，非变异或共线明确失败，不自动删项。

不报告逐细胞p值；各继承分析块内先平均文库系数，至少10块后做2,000次块bootstrap和2,000次符号翻转，跨所有可估计分区/模型作诊断性BH。翻转依赖块效应零假设下的对称性；这是筛选后诊断，不替换原发现FDR。来源未审阅的样本替代块不能称独立患者。身份、恶性、状态仍分开。

## 命令与代码证据

以下命令按阶段顺序执行；审计可等待已启动生产者，但不得在观察超时后自动重启：

```bash
cd "${PREFLIGHT_RUN_ROOT}"
python run_family_native_axis_jobs_v1.py native_axes
python run_family_native_axis_jobs_v1.py native_axes_audit
python run_family_conditioning_jobs_v1.py depth_member
python run_family_conditioning_jobs_v1.py depth_member_audit
python run_family_priority_jobs_v1.py priority
python run_family_priority_jobs_v1.py priority_audit
python run_family_spatial_panel_jobs_v1.py spatial_panel
python run_family_spatial_panel_jobs_v1.py spatial_panel_audit
python run_family_feasibility_jobs_v1.py feasibility
python run_family_feasibility_jobs_v1.py feasibility_audit
python run_family_receiver_programme_jobs_v1.py receiver_programme
python run_family_receiver_programme_jobs_v1.py receiver_programme_audit
python register_family_candidate_validation_milestone_v1.py
```

命令是执行结构记录，代码仍在外部运行目录；公共仓库未包含全量输入或所有历史runner。独立输入/结果审计与最后15条数值QR复算的独立程度分别说明，不能一概称独立生物验证。

## 失败与修正

旧all_fixed_member_mixtures模型的基准编号在部分来源不存在，存在编号比例与截距共线。按固定样本框的计数选择最后实际存在编号作基准，保留所有其他存在编号比例。新结构编号模型321可估计、15高影响不可估计；共同深度加发送者程序模型218可估计、118高影响不可估计。旧失败历史与不可估计案例保留。

部分SSH观察/上传超时；通过真实进程、终止manifest、审计文件和传输SHA分别确认，未把网络观察失败当科学任务失败或重启理由。一次公共日志导出误用本地汇总文件名，纠正为服务器实际覆盖表并独立计数后成功；不影响科学输出。

## 后续完成标准

对全部13接触记录覆盖的5份既有空间文件检查实际RNA、基因名、坐标单位、FOV/core/样本/患者映射与QC；不能把FOV或重复TMA材料视为独立患者。通过预检后做样本级邻近与保留密度/谱系/QC的零模型。参考支持仍需逐轴物理复合物、方向、糖基化/异构体等证据审阅。本家族扩展目前0直接原CM边匹配；尚未建立CM驱动轴。

输出和审计链接：[关键节点证据目录](../inventories/analysis/2026-10-11/README.md)。

## 2026-10-11 18:36:08 UTC+8 输入预检观察

按[只读计划](../inventories/analysis/2026-10-11/spatial_metadata_checkpoint_plan_v1.json)检查全部5份CosMx，采集执行时间10:34:46.528707—10:34:48.478582 UTC、exit0；这是元数据观察，没有新增空间模型或科学登记作业。采集器及转换器SHA见[JSON](../inventories/analysis/2026-10-11/spatial_metadata_checkpoint_v1.json)。

原CSV全行核对确认GSE299786肺腺癌TMA2索引列误入矩阵，故障及修复门槛见[incident](../logs/incidents/2026-10-11_spatial-csv-index-column.md)。尚未修改转换器/对象或进行邻近测试。下一有界任务为新版本输入修复及独立核验，然后继续全部13固定接触记录；原172节点、9CM、v20261010、发现BH不改变。
