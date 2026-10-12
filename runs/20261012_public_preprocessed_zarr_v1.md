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
