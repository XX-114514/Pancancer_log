# V4原生输入失败下载恢复

- 用户要求：重新下载失败数据，并编写应对后续失败的代码。
- 事件UTC：2026-11-08T12:11:39起重新核验/恢复；取证2026-11-08T12:19:30.401105+00:00，UTC+8 2026-11-08T20:19:30.401105+08:00，Melbourne 2026-11-08T23:19:30.401105+11:00。
- 输入：同6原生/18映射；18映射和首原生复用，只重下5失败成员。预期SHA/CRC/归档长度来自既有独立完整输入及官方ZIP库存，私有清单保留。
- 基线：050f4f80并发注释提交已快进保留；下载修复源代码SHA见恢复合同。
- 环境：本机Python3.9/NumPy2.0.2/AnnData0.10.9；服务器Python3.10.19/NumPy1.26.4/AnnData0.11.4；传输核心为Python标准库、POSIX锁。
- 参数：256KiB块，6次范围尝试/60秒，指数退避+抖动≤30秒，3轮/18小时，保守8GiB尝试字节/reader；全CRC/长度/SHA和原生坐标。
- 服务器逻辑位置：`${PROJECT_ROOT}/runs/20261108_V4_native_coordinate_fix_v1/download_recovery_v1`。大型cache、残片和矩阵留服务器；不删除、不改其他活跃上传。
- 实际状态：已部署并通过服务器故障测试，旧输入重新核验，五份恢复实际运行；最新完整数以[检查点](../inventories/analysis/2026-11-08/v4_download_recovery_v1/runtime_checkpoint_v1.json)为准。
- 完成：六原生通过且写ready才exit0；其余exit2且保留可恢复状态，不报成功。分析/科学验收另记。

[方法与重启命令](../methods/V4-Xenium-resumable-download-recovery-20261108.md) · [测试](../inventories/analysis/2026-11-08/v4_download_recovery_v1/software_verification_v1.json)。

有限接续控制已部署：48小时窗口/最多3次追加worker，只在原writer退出且未全部接受时运行，不停止活动下载；已接受则不重复。独立服务器进程以网络预算0复用真实缓存块通过。原分析消费者的期限仍单独核对，控制器不改分析参数或假装完成分析。
