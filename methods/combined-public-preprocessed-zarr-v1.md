# 已下载公共来源统一预处理 Zarr v1

## 授权与范围

用户2026-10-12允许备用卷作临时加工，要求继续至Zarr；最终对象在原Project_v5卷的新独立run。原4676787细胞/54GSE权威Zarr不改，下载仍暂停。本阶段是全量基因/元数据技术合并；没有新的批次校正、生物学身份或临床推断。前置为9来源逐文库全量独立audit PASS，不能以抽样或部分source代替。

## 存储与层

- 最终`public_preprocessed_union_v1.zarr`使用AnnData/Zarr2标准CSR编码，分块追加；不会构造全细胞×全基因稠密矩阵，也不复制R交换文本。
- `layers/counts`：来源发布计数，经本轮有审计账本的QC与双细胞选择后保留；整数int32，非整数/溢出停止。作者是否已校正的exact slot未知时，不称它未经处理的液滴计数。
- `layers/ambient_removed_counts`：float64稀疏**实际采用**的去除量；不采用的模型估计不填入。分析计数=`counts.astype(float64)-ambient_removed_counts`，避免全量重复保存未改变的校正计数。
- `X`：上述分析计数每细胞总量10000、log1p、float32；不复用选层前旧X。
- filtered-only DecontX须通过收敛/本谱系marker保存/UMI损失门槛，且单细胞去除≤50%才采用。LBCL、乳腺有作者污染处理记录但exact slot待对应，pHGG当前神经/胶质面板覆盖不足，均保留来源counts而将新模型作为诊断。高损失/不可辨识不是已证明无污染。
- Blosc zstd level5 bitshuffle，固定1Mi元素数组块；容量按每文库实测压缩率×完整输入非零数上界预估，另35%+12GB缓冲；原卷以物理余量/组硬配额较小者验收，保留**500GB十进制**空闲。写入时持续检查，不能用物理空闲代替配额。

## ID、缺失与排重

- 稳定cell ID、原始文库、技术、组织、富集和产品分支保留。source-local cluster仅为原文库诊断，不作全局细身份。
- 各来源基因符号outer union；`varm/measured_in_source`与`uns/source_order`给出测量覆盖。补零只表示未测，不能解释为生物学零表达。
- 元数据保留列覆盖账本；布尔采用可空类型（missing≠false），字符串缺失显式sentinel。数值通常float64；超过2^53的整数转字符串防止ID精度丢失，独立审计原值。
- 与冻结54GSE基线及新合集执行完整整数计数BLAKE2b192+同原始barcode指纹；保留基线基因词表顺序，新增基因稳定追加。barcode提取与原基线一致，DNA片段不含文库`-1`后缀；非匹配时保留原字符串，不编造barcode。
- 精确命中行进入私有隔离账本，源输入不删除。计数指纹不相同仍可能是同一文库的不同处理版本；研究/样本/文库别名复核未闭合前`source_final_AC=false`、正式新增AC=0，技术合并不解除此条件。

## 交付与审计

最终Zarr在原项目卷独立可读；其自身数组无临时卷软链接。小型逐文库QC诊断、原feature映射、重复排除账本随最终run保留；原输入/分片/交换文本仍在原位，用户本次没有授权删除。

独立审计逐个保留细胞检查：原count精确一致、元数据/顺序/ID唯一、选层/校正非负及上界、独立重算总量和归一化、未测基因从不非零、各CSR指针/shape、AnnData稀疏读接口。合成对象完整`read_zarr`验证格式与缺失/排重/大整数边界；不能把它冒充真实临床验收。真实输出以`independent_final_audit.json`和终态为准。

源码：[合并](../scripts/merge_public_preprocessed_zarr.py)、[独立验收](../scripts/audit_combined_public_zarr.py)、[等待/执行组件](../scripts/run_public_zarr_completion.py)。执行合同见[run](../runs/20261012_public_preprocessed_zarr_v1.md)。
