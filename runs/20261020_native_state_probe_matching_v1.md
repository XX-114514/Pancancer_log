# 原生空间相同探针口径及作者谱系RNA捕获

## 2026-10-19T16:12:33.984619+00:00 UTC /2026-10-20T00:12:33.984619+08:00 UTC+8 原生空间相同探针口径接受，RNA捕获运行

全部14输入、作者/几何4输入及原66代理位置保留；1556参考口径/48gene轴/全部目标映射和154技术视野计数独立接受。目标原生RNA任务已于16:11:02 UTC启动，无目标科学接受/新空间效果。[方法](../methods/native-spatial-state-probe-matching-20261020.md)、[run](../runs/20261020_native_state_probe_matching_v1.md)。goal active。

## 合同与环境

[匹配计划](../inventories/analysis/2026-10-20/native_state_probe_matching_v1/frozen_plan_before_matching_v1.json)、[原生捕获计划](../inventories/analysis/2026-10-20/native_state_probe_matching_v1/frozen_native_capture_plan_v1.json)、[当前技术证据](../inventories/analysis/2026-10-20/native_state_probe_matching_v1/technical_matching_acceptance_capture_running_v1.json)。实际输入以原清单SHA及源码SHA固定，commit在发布形成不倒填。既有scanpy环境NumPy1.26.4/SciPy1.11.4/anndata/h5py，单BLAS；没有安装或修改软件。

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/qualify_native_spatial_state_probe_contract_v1.py
python scripts/audit_native_spatial_state_probe_contract_v1.py
python scripts/run_native_spatial_state_capture_v1.py
```

执行前提供PANCANCER_PREFLIGHT_ROOT、任务tmpfs的PANCANCER_RAM_OUTPUT及对应PANCANCER_STATE_NATIVE_PROBE_PLAN/PANCANCER_NATIVE_CAPTURE_PLAN。匹配16:01:22.720906–16:01:37.717439 UTC /次日00:01 UTC+8；独立接受16:06:17.554065 UTC /次日00:06 UTC+8。新捕获16:11:02.115749 UTC /次日00:11 UTC+8启动，当前未最终接受；公式已通过合成fixture和原参考缓存实际比较，后续真实范围/退出及独立输出审核仍待。新空间/临床模型尚未运行，原身份/CM/Zarr保持。

## 已确认失败与输入域修正

capture v1在完成首图后、第二输入raw metadata全域相等断言失败，无完整终态。独立键集合和坐标/面积诊断支持唯一完整超集；未改变native表达对象或原注释。reader v3仅接受唯一完整超集并逐FOV保存额外记录数量，缺失/重复仍拒绝。真实读取器四类合成归档与原数值fixture通过；[capture v2计划](../inventories/analysis/2026-10-20/native_state_probe_matching_v1/frozen_native_capture_plan_v2.json)源和helper SHA锁定，16:16:43.001601 UTC /次日00:16 UTC+8一次启动。生物门槛、原4输入与预期190038原生细胞不变，尚无全量接受。

## 原生生产完成，独立审核启动

v2于16:18:42.987874 UTC完整终止，技术终态efa79e09ea535391339ec81d4ecdbb0088500cea23d9ef003814d8c020ea2a1d；全190038细胞/190038000原生计数位置、负探针和metadata、31代理91参考值及189目标标量接受。16:21:13.658166 UTC独立审核启动；2426视野代理/4846细亚型及source汇总尚待审核/科学分级，未运行空间或临床模型。

## 最终输出/科学与视觉接受

独立审核16:21:32.195722 UTC完成95693e…0999c，全2426视野代理/4846细亚型与source统计重建，最大误差5.55e-16；原始表达归档未在审核中再流式逐行重放，范围与生产完整raw核验分开。来源报告16:26:41.495194 UTC完成a645f0…4b32f；科学/PNG实际像素审阅2026-10-19T16:34:00.993016+00:00 UTC /2026-10-20T00:34:00.993016+08:00 UTC+8。详见[完整技术证据](../inventories/analysis/2026-10-20/native_state_probe_matching_v1/technical_native_capture_final_acceptance_v2.json)。原生捕获阶段已结束，下一为真实有限空间预测；目前无新空间/临床效应。现有官方MISTy全合成软件检查通过，默认KFold/OOB不是本项目外层空间验证。
