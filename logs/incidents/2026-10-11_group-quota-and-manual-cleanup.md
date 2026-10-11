# 数据卷组配额耗尽与手动清理准备

## 已观察事实

取证UTC 2026-10-11T14:09:31.866016+00:00 /UTC+8 2026-10-11T22:09:31.866016+08:00。data4卷df仍有约3.7T可用，但分析组quota已达到1500G硬上限；用户个人quota显示none不能排除组限制。资源aria2于21:51:30 UTC+8写cache flush返回errno122/Disk quota exceeded；worker实际13:51:38.494517 UTC以exit1终止，13:55:22.652731 UTC确认worker/child均不存在，v16/active_jobs空/状态0冲突。资源恢复尚未完整，未启动全基因生产或审计。

13:53观察只读了旧progress尾行79MiB/31%，没有同时确认manifest/PID；先前日志“仍为同worker且增长”的表述错误。本次明确更正，最后可靠活跃进程取证是13:50:09.711972 UTC；不能从文件大小或旧尾行推出正在运行。历史原日志保留。

Git同步原data4仓库在fetch/index-pack写临时pack时同样配额失败；本机HTTPS普通push无凭据失败。保留所有原仓库与失败字节，在/tmp新建no-hardlinks Git副本，正常快进/普通SSH push，4fcceccc50f0462ea9b6592343193bfcf9265b56于13:56:46.217152 UTC推送，13:56:53.175751 UTC独立本机HTTPS核验同SHA。[回执](../../inventories/analysis/2026-10-11/quota_recovery_git_receipt_v1.json)。原data4 Git未改写，临时Git不是数据长期存储。

## 用户指定的处置方式

用户要求只识别冗余、提供代码，由用户手动删除。代理没有执行任何清理候选的unlink、删除模式或去重替换。仅读大小/硬链接/SHA/自身进程FD及maps，生成清单及默认preview脚本；删除必须用户交互输入对应APPLY确认词。

| 候选组 | 文件数 | 预计可释放GiB | 证据与影响 |
| --- | ---: | ---: | --- |
| Conda压缩包归档 | 602 | 1.033 | 全部有解包目录、原包SHA冻结；仅删.conda/.tar.bz2，环境/解包目录保留，后续重装可能需下载 |
| 失败Range数据块 | 158 | 0.148 | 任务已终止、无使用；只删.bin，保留48份failure.json及其他收据；损失旧Range续传字节 |
| TF失败截断副本 | 1 | 0.006 | 完整校验TF另存且SHA检查；完整LR/TF及当前aria2 target部分/控制文件保留 |
| 重复GPU传输文件 | 21 | 3.774 | 与Project_v5保留副本完整SHA逐字节一致、独立inode/nlink1；可选去重将旧副本换相对symlink，保留旧路径，之后须避免原地修改两处共享文件 |

合计预计4.961GiB/5.327GB，实际释放以quota复查为准；没有确认更大历史目录可删。Project_v1约394GiB/Project_v2约327GiB/Project_v3约530GiB，含原始数据/参考/历史结果。最大两个v3run约180.7GiB与147.8GiB，仅占用定位，不列整目录删除；权威Zarr/172身份/9CM及新下载原始来源继续保留。

## 代码及检查

[通用代码](../../scripts/manual_cleanup_reviewed_v1.py)、[脱敏计数/SHA](../../inventories/analysis/2026-10-11/manual_cleanup_reviewed_summary_v1.json)。逐文件完整路径/空间文件名与白名单只留服务器/tmp和用户本机，不上传公共仓库。手动入口/tmp/pancancer_manual_cleanup_20261011_v2.sh：preview只读；cache清理前三组；dedup可选处理21重复副本。默认preview不改候选。

首次生成器将206个Range目录文件都当.bin，断言拒绝；进一步查明158数据块+48错误JSON，仅158列入。首次只读preview在device/dev字段差异处拒绝，未进入删除模式；保留v1白名单/入口，新v2规范dev字段。最终真实服务器preview exit0，782个候选SHA/保留副本/进程占用核对通过，输出“no file deleted or replaced”；终止事件精确时间unknown，观察2026-10-11T14:09:31.866016+00:00。未用真实候选运行--apply，删除行为仍由用户执行。

下一步：用户手动处置并复查组quota，再以新的资源版本恢复；不重启已失败v3或用旧31%日志称健康。正式分析目标未完成，v16保持；全45源header合同未实际执行。
