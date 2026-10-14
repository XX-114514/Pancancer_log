# Spectra多初始化共享及身份内状态发现 v1

UTC 2026-10-14T01:09:01.038759+00:00/UTC+8 2026-10-14T09:09:01.038759+08:00。[固定合同](../inventories/analysis/2026-10-14/Spectra_multiseed_state_discovery_plan_v1.json)。输入独立全240533raw rows/108251626gene-unit聚合值PASS；4源/12000source-major/block均衡发现细胞，胃source为内部已用来源投影。共同轴16625，full-native RNA分母。

官方commit9c079...隔离GPU实现，3seed/每最多1000epochs/1000minibatch，lambda.008/delta.001/rho.05，固定143global/18T/5Mac参考，加3/2/2free，共173factor slots；discovery-only每major2000HVG并集+完整prior基因。不删除全LR、不使用响应或通讯p选特征/超参数。候选factor数目不是正式状态数量。

运行入口：隔离venv `python train_spectra_state_candidates_v1.py`，PANCANCER_PREFLIGHT_ROOT/PANCANCER_SPECTRA_PLAN提供根及冻结。模型/全训练scores/markers/checkpoints私有，新目录拒绝覆盖。状态：predeclared_not_started，code语法通过；多seed输出需独立轴/数值/稳定性及原身份/来源/RNA诊断，随后cNMF/UCell/9CM/DIALOGUE/空间，不能因fit成功完成goal。

## 实际终态：2026-10-14 09:37 UTC+8取证

3seed真实01:14:25—01:31:39 UTC exit0，两个PID已不存在；每173factors/6512features/12000cells，完整输出SHA/finite/轴检查通过。146shared/20T/7Mac是candidate slots，不是正式state。状态completed_fit_pending_crossseed_source_identity_diagnostics；cNMF实际同时运行，UCell/所有9CM/状态后DIALOGUE/空间待。[checkpoint](../inventories/analysis/2026-10-14/state_models_progress_v1/RESULTS.md)。
