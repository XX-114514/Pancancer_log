# 正式运行记录

当前结果入口：[原生全家族解释准备](20261018_NSCLC_native_interface_result_review_v1.md)，绑定原全数值审核，等待通过后才读取效果。

当前新增入口：[原生接口外部测量机会预检](20261018_NSCLC_native_external_opportunities_v1.md)，与仍运行的原生完整关联家族分开登记。

本目录保存正式运行的轻量级可复现记录，不保存分析对象或完整程序日志。

最新CM复现支线：[20261012 CM定义与免疫来源恢复](20261012_CM_immune_recovery_v1.md)，4源131库/6分层，54捕获检查与432重连/分母敏感性；内部骨架支持不等于完整或独立CM验证。

创建方式：

```bash
bash scripts/new_run_record.sh "phase07 annotation"
```

记录文件按年份组织，命名为 `YYYYMMDD_HHMMSS_run-name.md`。创建后应：

1. 在运行前填写输入、代码 commit、环境、参数、命令和完成标准；
2. 在运行后填写结束时间、退出状态、输出验证、错误和解释；
3. 在 `inventories/runs.tsv` 增加一行；
4. 重要产物加入 `inventories/artifacts.tsv`；
5. 不复制 `${RUN_ROOT}` 中的完整日志或大型对象。

已有的外部历史运行可先在 inventory 建索引；只有在需要深入审计时才补建对应 Markdown 记录。


## 2026-10-12 新公共来源完整技术合集

[稳定基因ID预处理Zarr v2](20261012_public_preprocessed_zarr_v2.md)：9来源/792文库/3180973细胞，全量独立输出审计、标准读回及持久证据PASS；来源正式AC与生物学整合仍待。已登记runs/artifacts索引。

## 2026-10-14 全谱系细胞状态输入与模型

[全100766原生输入验收](20261014_all_lineage_state_native_inputs_v1.md)、[6block Spectra/18major cNMF真实运行](20261014_all_lineage_state_models_v1.md)。原172/9CM保留；状态资格、全数据投影、全CM/DIALOGUE/空间仍待。


## 2026-10-16 原CM全LR与接收功能

[全LR验收/TF功能启动](20261016_DLG_LR_receiver_function_v1.md)：4012估计、14BH/6区间，5448固定功能槽运行；原身份/CM保留，机制仍待。


## 2026-10-16 原CM接收功能与空间技术校准

[5448TF联合证据](20261016_DLG_TF_joint_evidence_v1.md)、[全52来源空间代理](20261016_state_spatial_proxy_calibration_v1.md)。TF全量独立通过；空间块均值精度修复后独立审计仍运行，不能当空间生物学验证。


空间状态程序最新：[全来源代理校准终态](20261016_state_spatial_proxy_calibration_v1.md)；[GSE273952全转录组冻结投射](20261016_GSE273952_frozen_programme_spatial_v1.md)。
