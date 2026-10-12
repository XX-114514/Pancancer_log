# 稳定基因ID的公共来源统一预处理Zarr v2

## 输入与当前状态

沿用已授权9来源/792原始文库及同一QC/DD/DecontX生产；用户允许临时备用卷、最终原项目卷。初始v1按名称设计存在真实基因别名拆分，04:23:12 UTC停止尚在等待的组装器自身，**没有科学矩阵写入、没有停止前处理**。旧attempt保留。[v1历史](20261012_public_preprocessed_zarr_v1.md)。

新final run `20261012_public_preprocessed_zarr_v2`，`${FINAL_ROOT}/outputs/public_preprocessed_union_v1.zarr`是该次完整对象；`union_v1`为存储schema名字，不是旧等待器产物。代码普通push/冻结后才启动，目前未启动或完成真实合并。

## 冻结与门槛

- 原科学生产Git `cdcea29b45e989cb8d03235b9a92521611b6dd17`、完整审计v3 `37eaefef9d920c156e30c33e1f0a8c84d619c5565e80c74fe122bfaf305fe65f`不变。
- 组装包含merger、independent auditor、orchestrator、稳定gene identity helper四个SHA；以实际start冻结。全9源audit PASS前不开始科学写入。
- Python3.10.19/AnnData0.11.4/Zarr2；BLAS1、压缩2线程，CSR分块、实测配额/500GB十进制留空。
- 原始输入基因名称72554→稳定特征49985；400歧义特征、1244仅symbol特征显式区分。int64别名求和后验收int32，原count总量逐cell守恒，合并后再normalize。
- 纯人工v3 108行，基线及新源重复各1隔离→106行/18特征；别名求和/歧义/守恒/规范化/大整数/完整AnnData读回及独立全行PASS，真实临床细胞0。不能冒充真实9源验收。

[当前方法](../methods/combined-public-preprocessed-zarr-v1.md)。命令结构同v1，`${FINAL_ROOT}`替换v2独立新run；scratch独立v2。所有failed/superseded attempt保留，无删除、自动定时push或新增下载。

## 最终验收

输出counts/采用去除量/X/元数据/gene mask；按源上下文限制用途。逐cell独立COO投影求和/总量/归一化/ID/元数据检查，完整输出后才称完成。正式来源别名AC仍0，污染不可辨识/作者槽未知/CAR-T产品/富集snRNA等限制保留；物理技术合并不代表新批次模型或新身份。


## 2026-10-12T12:28:42+08:00 稳定ID组装v2已真实启动

- d8a53ee4feecff2129b70cc1564c90f706e667dd先普通push、独立HTTPS SHA一致；四组件源码SHA与合成v3相符。04:28:42 UTC观察560/792文库、0失败、2380251保留细胞，剩余3源仍处理。v2启动记录与manifest核实，原卷新独立run，临时scratch独立v2。
- 新组件等待原全9源审计PASS，然后精确指纹隔离→稳定ID CSR合并→全行COO独立验收。旧v1等待器已停止保留，不是两个科学组装同时执行；科学队列与DD/污染参数不变。
- 新流程仍不下载、不批次模型、不做生物学身份；最终独立可读对象与来源最终AC/污染限制分开汇报。R交换等历史不自动删除。
