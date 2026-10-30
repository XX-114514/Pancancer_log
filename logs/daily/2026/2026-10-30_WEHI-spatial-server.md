# WEHI数据服务器建立与空间数据下载

2026-10-30 13:33 Melbourne / 02:33 UTC：WEHI数据服务器接入本仓库（deploy key、普通push）。WEHI代码移入`scripts/wehi_xenium/`，路径改由本地配置提供。仓库检查环境建好，620个测试只剩1个已知环境缺口。
Xenium共31单元入库，其中新增S18 HNSCC 10位患者。GEO的CosMx、MERFISH、Stereo-seq、Visium和三个单细胞参考共415个文件下载完成；IMC归档仍在下载。
发现S16、S20、S22是同一批TMA的连续切片，不能作为独立患者；另外纠正了两条服务器笔记（GSE328481的平台；GSE273952/GSE316195不是配对数据）。
详见[任务日志](../../codex/2026/2026-10-30_WEHI-server-setup-and-spatial-data.md)、[数据汇总](../../../inventories/runtime/2026-10-30/wehi_spatial/DATA_SUMMARY.md)、[规划](../../../methods/WEHI-spatial-line-plan-20261030.md)。

2026-10-30 18:30 Melbourne / 07:30 UTC：测量了WEHI与旧服务器之间的通道（SSH约4–30 KB/s），按下一阶段依赖规划分层迁移（M0–M3），约200 GB公开或历史数据不迁移。见[规划](../../../methods/WEHI-data-migration-plan-20261030.md)和[任务日志](../../codex/2026/2026-10-30_WEHI-data-migration-plan.md)。
