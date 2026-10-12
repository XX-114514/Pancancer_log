# 公共来源统一预处理 Zarr v1

## 当前状态

用户已授权临时加工及最终Zarr；代码/合成测试/实际压缩预检已准备，真实全量合并**尚未启动或完成**。等待既有9源792文库全部处理及独立审计，不能以输入3338473细胞当最终保留数。

- 既有生产冻结`cdcea29b45e989cb8d03235b9a92521611b6dd17`；全量QC审计v3 SHA`37eaefef9d920c156e30c33e1f0a8c84d619c5565e80c74fe122bfaf305fe65f`。
- 暂存`${PREPROCESS_ROOT}`在个人备用项目区，最终`${FINAL_ROOT}`在原Project_v5独立run `20261012_public_preprocessed_zarr_v1`。
- 科学环境Python3.10.19/AnnData0.11.4/Zarr2/Scanpy1.11.5；BLAS/Numba1线程，压缩2线程。
- 新组装代码Git SHA以实际`assembly_start.json`冻结；启动前必须普通push并核验。组件不重跑/改变既有双细胞和DecontX生产。

## 实测容量预检

04:12:41 UTC/12:12:41 UTC+8，498完成文库/9来源/2115465保留细胞的实际压缩率，按完整原输入非零数上界缩放：约46.88GB矩阵，含35%+12GB余量约75.29GB，原卷配额可用593.64GB，保留500GB后可容纳。此为上界估算，非最终文件大小；最终写前全792再估，写入持续检查。

证据：[容量](../inventories/runtime/2026-10-12/combined_zarr_compression_preflight_v1.json)，[合成读回](../inventories/runtime/2026-10-12/combined_zarr_synthetic_contract_test_v1.json)。[方法](../methods/combined-public-preprocessed-zarr-v1.md)。

## 冻结命令结构

```bash
${PYTHON} -u ${FINAL_ROOT}/scripts/run_public_zarr_completion.py \
  --preprocess-root ${PREPROCESS_ROOT} --collection ${COLLECTION_ROOT} \
  --final-root ${FINAL_ROOT} --scratch-root ${TEMP_SCRATCH_ROOT} \
  --code-git-sha ${FROZEN_CODE_SHA}
```

独立进程/锁；每阶段写进度、真实PID、exit code和组件SHA，终态需审计PASS。所有attempt保留，不覆盖已存在目标；Git日志须审阅后普通push，非自动定时发布。矩阵、患者/cell行和内部路径留服务器。

## 验收/用途边界

最终`outputs/public_preprocessed_union_v1.zarr`；原count、实际采用校正差值、重新归一化X、元数据、gene mask齐全。全行独立核验后才可称技术对象完成。作者预处理槽/脑谱系污染可辨识性/来源别名的限制单独保留；技术合并不增加正式新增AC、不声称所有污染或双细胞完美清除，不重新注释、不作临床合并推断。

追加合成v2于04:15:52 UTC全部PASS：旧基线精确命中1行、新合集重复1行隔离，108→106，完整AnnData读回/全行独立核验/超过2^53整数ID原值保留通过。真实临床细胞0；v1保留。代码SHA与v2报告一致。Cloud轻量5测试通过，0错误、1既有路径warning。


## 2026-10-12T12:18:34+08:00 统一Zarr组装组件已启动

- 3ed9acb7b9b79f3b17a3639374081a5e0c6db0fc代码/合同先普通push并独立HTTPS核验一致；新组件04:18:07 UTC实际启动PID39340，3个源码SHA与冻结一致、进程存活/等待manifest核实。最终目录在原项目卷的新独立run，暂存仅加工用途。
- 04:18:34 UTC真实528/792、0失败、2236179保留细胞、82773额外预测双细胞移除、1012967严格技术准备候选；6来源步骤结束，余3仍处理。完整9源审计通过后自动进入确切指纹/流式合并/全行独立验收，未声称当前已得到最终Zarr。
- 合成v2增加旧基线匹配与大整数ID保护均PASS。组装launch只含执行汇总公开；真正count/患者行/内部路径不上传。R交换/分片不自动清理。新下载仍暂停，来源最终AC保持0。
- 既有记录finalizer v2仅维护准备索引，本次新组装执行授权的物理技术合集；二者不同角色，旧记录中未整合/user scope门槛是当时状态，当前用户授权已覆盖物理组装。研究别名/污染可辨识性门槛仍保留。

## 2026-10-12T12:23:12+08:00 当前attempt被稳定ID修订替代

等待器39340经命令/独立进程组核实后停止，仅组装等待角色；尚未写科学矩阵、原预处理继续。真实基因别名拆分需稳定ID统一，后续转[独立v2](20261012_public_preprocessed_zarr_v2.md)。原start/stop/代码历史不删除。此状态取代上方“正在等待”历史。
