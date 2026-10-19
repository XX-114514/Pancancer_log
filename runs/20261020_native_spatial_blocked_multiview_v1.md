# 原生空间分区多视图：真实有限模型

## 启动

事件 2026-10-19T16:54:59.613156+00:00 UTC；取证 2026-10-19T16:55:55.968495+00:00 UTC /2026-10-20T00:55:55.968495+08:00 UTC+8。

[方法](../methods/native-spatial-blocked-multiview-20261020.md)、[效果前计划](../inventories/analysis/2026-10-20/native_spatial_blocked_multiview_v1/frozen_plan_before_real_models_v1.json)、[真实fixture与启动](../inventories/analysis/2026-10-20/native_spatial_blocked_multiview_v1/technical_fixture_and_launch_receipt_v1.json)。原4输入/63目标和两读数全部126位置。当前真实生产运行，结果/独立审计/科学接受待；不是空间互作机制或临床验证。

## 执行与资源

原scanpy Python3.10、LIANA1.7.1、sklearn1.7.2、NumPy1.26.4/SciPy1.11.4；无需安装。已审通用源码复制任务tmpfs后绑定SHA，当前发布commit另记、不倒填。提供任务目录PANCANCER_MV_CONTROL和其private_output为PANCANCER_RAM_OUTPUT；其plan含私有合同路径/完整SHA。原项目只读。

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python scripts/run_native_spatial_blocked_multiview_v1.py
```

两image worker，每RF单线程。tmpfs至少512MiB空闲/汇总32MiB上限；RNA/评分/图/模型/逐细胞预测只RAM，私有persist仅模型/视野/fold误差、覆盖、技术和重要性汇总。stdout不公开。跟随原launch receipt PID/start ticks，执行故障不当科学阴性，不因SSH观察超时重启。

## 验收待办

全部126家族守恒、所有fold/输入与状态；独立逐视野误差重聚合、raw/负探针对照及实际模型复核；来源科学分级/实际图像。所有科学效果/程序/个体保留私有，完整输入和不可估计不隐藏。


## 2026-10-19T16:58:33.336863+00:00 UTC /2026-10-20T00:58:33.336863+08:00 UTC+8 真实生产完整结束、独立核算启动

全部4输入/126真实位置已于2026-10-19T16:57:19.409377+00:00 UTC结束，完整终态e187616ee043701294421c002444973803684be0a1fd34d4a050c4fa4d5e0bbe，无执行失败；这不表示全部可估计或效果阳性。2026-10-19T16:58:00.409977+00:00 UTC独立审核一次启动，范围为全家族/fold/视野误差与固定第一输入第一目标fold0的两个读数，直接sklearn RF/显式中心化Ridge复跑。尚无科学接受。[方法](../methods/native-spatial-blocked-multiview-20261020.md)、[run](../runs/20261020_native_spatial_blocked_multiview_v1.md)。目标active。

审核源8e3519543dac425957d9386c21192c63229f57c635471b3a37854cfb7d734aee；独立汇总另一pandas实现逐视野重聚合，直接RF与显式Ridge正规方程复跑原固定第0图/第0目标/fold0、两读数，选择时未审效果。输入变换复用已SHA绑定/原生捕获接受代码，不冒称第二次全原始RNA重放或全部RF重复拟合。当前不导出科学效应。
