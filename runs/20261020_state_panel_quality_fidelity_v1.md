# 现有面板代理质量调整保真度 v1

## 2026-10-19T15:34:43.173206+00:00 UTC /2026-10-19T23:34:43.173206+08:00 UTC+8 面板外参考的深度调整测量问题实际运行

全部91447细胞原分数/metadata，709原生来源大类文件、2836面板分数文件；固定质量控制与块等权QR/Schur，原完整空间门槛保持。正负/常量/共线/不等权fixture及3.10 AST通过，15:33:21 UTC一次启动，无新科学接受。[方法](../methods/state-panel-quality-adjusted-disjoint-fidelity-20261020.md)，[run](../runs/20261020_state_panel_quality_fidelity_v1.md)。goal active。

## 输入、环境与复现

计划及全部输入清单SHA见[冻结合同](../inventories/analysis/2026-10-20/state_panel_quality_fidelity_v1/frozen_plan_before_adjusted_fidelity_v1.json)。复用scanpy环境NumPy1.26.4/SciPy1.11.4及anndata/h5py，只读既有metadata和分数，不读原RNA值或临床/空间效果。生产源SHA86bf00a4b4a05eeff0d20d51b7032219995fee055a12175767873a8c8788a4f0；代码commit在发布节点形成，不倒填运行时Git版本。

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python scripts/qualify_state_panel_quality_adjusted_fidelity_v1.py
```

调用前提供PANCANCER_STATE_MODEL_OUTPUT_ROOT、PANCANCER_RAM_OUTPUT和PANCANCER_QUALITY_PROXY_PLAN，后两者为独占任务tmpfs新路径。15:33:21.031006 UTC /23:33:21.031006 UTC+8启动，进程真实存活；退出/完整结果仍待。正式标签/原2026年10月CM（v20261010）/Zarr及已闭合检验保持。个体数组、学习基因及科学系数不公开。

## 数值修正及当前实际执行

v1常量块合成反例失败，未接受其块分流；历史生产继续留存。相同统计规则的v2新增块均值相对方差保护，15:36:48.071332 UTC /23:36:48.071332 UTC+8一次启动；新源SHA47a8774f645b0a759c31dce3626c1dfdff0fa9d1cf6e74d5c632dcac431692b2及[新计划](../inventories/analysis/2026-10-20/state_panel_quality_fidelity_v1/frozen_plan_before_adjusted_fidelity_v2.json)保留完整修复轨迹。命令改用scripts/qualify_state_panel_quality_adjusted_fidelity_v2.py；原科学门槛不变。源/计划SHA绑定实际运行，不从观察超时重启。

## 完成与接受

生产15:40:37.518077 UTC、全来源/分流审核15:43:34.251854 UTC、来源图15:45:01.275094 UTC完成；真实科学与像素审阅2026-10-19T15:49:38.017844+00:00 UTC /2026-10-19T23:49:38.017844+08:00 UTC+8。全范围和SHA见[技术终态](../inventories/analysis/2026-10-20/state_panel_quality_fidelity_v1/technical_complete_quality_fidelity_acceptance_v2.json)。科学系数/分流/图保持私有；原完整门槛/正式身份和旧家族不变。下一原生目标捕获/代理去重后另立空间问题，本轮没有新空间效果或通讯发现。
