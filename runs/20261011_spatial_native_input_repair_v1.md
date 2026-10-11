# 空间原生计数输入修复里程碑

取证：2026-10-11 10:52:59.350665 UTC / 18:52:59.350665 UTC+8。正式登记：10:52:17.951000 UTC。当前结果是输入修复完成；通讯轴验证仍未完成。

## 问题、范围和版本

GSE299786肺腺癌TMA2的Unnamed: 0列是CSV行索引，却进入表达矩阵。筛查全部14份已接受空间输入的基因轴，只发现这份对象含已知无名导出索引；这不等于所有基因别名均获生物学验证。生成独立新对象39123×1000；旧39123×1001对象和原转换器保留。

5份拟检验CosMx文件共254,458细胞、13条固定接触记录/9种LR组合保持不变；2条程序支记录独立保留。证据登记v12仅接受修复输入，不是新的注释V12。原v20261010、172节点、9CM、发现BH及权威Zarr未作修改。

输入：[原故障证据](../inventories/analysis/2026-10-11/spatial_index_column_incident_v1.json)、接受14文件库存、原H5AD、原GSE299786归档、自动15记录表及SHA锁定的原转换器。计划分别见[v1](../inventories/analysis/2026-10-11/spatial_native_input_repair_plan_v1.json)、[v2](../inventories/analysis/2026-10-11/spatial_native_input_repair_plan_v2.json)。全部汇总/输入SHA/代码SHA在[终态证据](../inventories/analysis/2026-10-11/spatial_input_repair_checkpoint_v1.json)。环境沿用scanpy Python3.10.19/NumPy1.26.4/pandas2.3.3/SciPy1.11.4；没有新增软件环境。

## 实际轨迹

| UTC+8时间 | 节点 | 结果 |
|---|---|---|
| 18:42:29 | 预声明原计划 | 14轴筛查、1已确认对象、新版本目录、原CSV全量比对 |
| 18:46:03.881285—18:46:06.905067 | v1生产者 | exit1；我错误假定无counts层，未生成修复对象 |
| 18:46:34.008617 | v1依赖审计结束 | exit1，未运行计算；开始时间未持久记录，不猜测 |
| 18:48:06.641921 | v2计划 | 源counts与X已核实完全相同；追加两层保留与修复 |
| 18:49:21.291248—18:49:39.034809 | v2生产者 | exit0，新对象生成；旧对象不变 |
| 18:49:51.378339—18:50:05.782303 | v2独立审计 | exit0，原CSV全39123×1000计数一致 |
| 18:52:17.951000 | 交接/登记/状态刷新 | v12证据登记；0冲突；本里程碑无活跃作业 |

v1计划提交a0a2d59于18:46:05.047925 UTC+8普通push，启动后1.17秒收到回执；v2失败/计划提交7da773f于18:49:04.168076推送并独立核验，在v2运行前完成。准确UTC事件、退出状态和脚本SHA见终态证据；不伪造当时未保存的时间。

## 修复和独立验收

新[索引守卫](../scripts/cosmx_export_index_guard_v2.py)在原转换器的列选择阶段排除无名导出索引，保留全零基因和控制探针。使用该守卫重新读原归档，再构建保留旧obs/var/坐标的新对象；X及counts同时修复，追加显式重算RNA总量/检测基因数，无细胞过滤或新身份赋值。

审计用另一套CSV分块读取及FOV＋cell-ID键，未调用原转换器barcode或列选择函数。全39,123细胞/1,000基因（39,123,000个值）与原CSV逐值一致；两层也与旧对象保留基因列一致。已有obs/var/坐标及原文件/归档/转换器SHA核查通过，13候选键保持。

修复后RNA总量中位数487，27零RNA细胞、701少于30计数；少于30只是描述，不是本次过滤规则。原计数中位数20,179包含索引，不能解释为RNA总量。

执行结构（路径均为外部逻辑路径）：

```bash
python "${PREFLIGHT_RUN_ROOT}/spatial_native_input_repair_code_v1/run_spatial_repair_jobs_v1.py" spatial_repair
python "${PREFLIGHT_RUN_ROOT}/spatial_native_input_repair_code_v1/run_spatial_repair_jobs_v1.py" spatial_repair_audit
python "${PREFLIGHT_RUN_ROOT}/spatial_native_input_repair_code_v2/run_spatial_repair_jobs_v2.py" spatial_repair_v2
python "${PREFLIGHT_RUN_ROOT}/spatial_native_input_repair_code_v2/run_spatial_repair_jobs_v2.py" spatial_repair_v2_audit
python "${PREFLIGHT_RUN_ROOT}/register_spatial_input_repair_milestone_v1.py"
```

作业、终止manifest、审计分别核查，没有因观察超时重启。运行代码在服务器，公开仓库只保存守卫和SHA/轻量汇总；生产运行前Git锁定不覆盖所有外部脚本，不能用事后日志commit冒充生产代码commit。

## 接下来如何进入验证

登记v12接受`${FAMILY_RUN_ROOT}/spatial_native_input_repair_v2/analysis_input_manifest.tsv`的analysis_path作为后续输入；其中原库存基因数/统计字段及初始pending状态属于历史输入/生产记录，不能覆盖独立审计和最新登记；必须重新读取analysis_path的var/X。

下一阶段覆盖13×5=65项，不手选替代记录；检查独立谱系/命名基因及LR同细胞表达、背景、QC和可测程序，保留零值、缺失及交叉谱系情况。随后核实组织/core/患者、空间单位和边界，才拟定邻近推断。文件/FOV不等于患者。

作者公开资料提供QC后细胞标签和TMA前缀，尚不能据此填患者映射。[作者说明](https://github.com/KChen-lab/ST_Comparison/blob/main/Annotations/READ.me)、[处理代码](https://github.com/KChen-lab/ST_Comparison/blob/main/Seurat/CosMx_ICON.R)。当前未下载或分析这些RDS。

此次没有证明配受体结合、空间接触、接收程序机制、CM驱动或跨患者复现；全目标继续推进。
