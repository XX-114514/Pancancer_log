# IFNGR1/2同细胞RNA证据 v1

## 目标

在原五源212库的冻结T/巨噬范围中，直接扫描权威Zarr的原表达X，统计IFNG、IFNGR1、IFNGR2及同一细胞R1/R2均>0。当前单细胞metadata已核240533细胞、1832样本节点，所有计数与缓存一致，无新身份/状态或原9CM改写。

## 范围与验收

[冻结输入/代码SHA](../inventories/analysis/2026-10-13/IFNG_receptor_single_cell_plan_v1.json)。原输入合同、节点字典、五源cache manifest固定；逐源gene→core_column精确连接，拒绝符号猜测；16384行CSR窗口，只读含目标细胞的范围。完整三基因边际阳性与原缓存逐1832单位核对，再独立重建全部联合聚合并新读各GSE/大类/阳性模式代表行。不声称第二遍全体细胞独立联合表达扫描。private cellbits/sample-node/sample-major留服务器，仅发布GSE/亚群聚合。

RNA共同阳性不保证蛋白受体复合物；未检测不表示蛋白缺失。独立来源、空间或功能机制尚未验证。

## 执行

服务器scanpy环境/BLAS2线程；PANCANCER_PREFLIGHT_ROOT/ANNOTATION_ROOT、PANCANCER_MAC_T_POOL_PLAN、PANCANCER_SINGLE_CELL_PLAN提供入口，分别执行mac_T_IFNG_receptor_single_cell_v1.py与audit_mac_T_IFNG_receptor_single_cell_v1.py，独立execution manifest。真实表达扫描尚未启动。
