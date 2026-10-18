# CM02空间分子测量与新增肺癌来源资格

## 执行前范围

- 2026-10-14免疫图谱CM（v1，T_whole_TME_primary_tumor）CM02，旧312临床及18条件空间槽关闭保持。
- 新83232位置＝3原MERFISH技术输入×6方向×完整4624资源，不新增P/FDR，不择强接口。
- 输入前提：原3输入全部RNA/几何旁表验收、完整原功能面板；LIANA consensus的既有参考证据表。
- 原输入/身份/W/S/H只读；新目录`${PREFLIGHT_ROOT}/CM02_spatial_molecular_opportunities_20261018_v1`，64MiB上限、500GB+1GiB预留。
- [效果前合同](../inventories/analysis/2026-10-18/CM02_spatial_communication_extension_v1/frozen_molecular_measurement_plan_before_effects.json)，[方法](../methods/CM02-spatial-communication-extension-20261018.md)。

## 可复现命令

```text
python measure_CM02_spatial_full_LR_opportunities_v1.py --root ${PREFLIGHT_ROOT} --plan FROZEN_PLAN --output ${NEW_OUTPUT_ROOT}
```

原scanpy Python3.10环境，单BLAS线程。生产完整计数后仍需独立原RNA/复合体/空间边接受。所有科学表和细胞/图数组保留服务器；公共仅技术收据、代码、合同和来源。

## 新公开数据资格

GSE271689/GSE221733官方SOFT已实际获取，原作者仓库4个小型处理表共9.41MB，逐文件Git blob与自算SHA接受。全固定34基因覆盖检查表明矩阵有缺失，不直接套用原临床代理。作者完整蛋白h5ad只给Google Drive链接，未下载未知大小对象；公开原始DCC schema已读，不假定RTS/probe等同gene symbol。Xenium TLS新研究需reasonable request，未联系作者或提交申请。

本节点尚未登记新83232位置终态或直接通讯证据；失败/不可估计位置全部保留。后续事件与精确源码/数据SHA写同任务[日志](../logs/codex/2026/2026-10-18_CM02-spatial-communication-extension.md)。
