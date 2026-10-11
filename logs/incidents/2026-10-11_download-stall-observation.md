# Incident：第二批部分下载停滞观察

- 时间：2026-10-11 03:22:53 UTC
- 状态：observed_pending_diagnosis；不是已修复或全队列故障
- 证据：[快照](../../inventories/runtime/2026-10-11/runtime_snapshot_v2.json)

GSE315928 RAW tar为12%，最近传输日志超过10分钟未更新。GSE294273为68%、GSE236581 counts为10%，采样瞬间0B/s。批次下载Python队列及部分aria2进程仍存在，Zenodo为36%且24KiB/s；进程存活不等于每个文件健康。

可能原因包括来源限速、连接重试或调度等待，当前证据不能定位根因。不要直接改端口或宣布代理故障。下一步对该文件核对队列重试记录、当前子进程、文件大小时间与HTTP响应；有断点续传恢复/官方校验后才关闭。保留partial及失败历史，不删除或把其计作完成。本次记录工作未重启/改动下载器。

## 2026-10-11 10:56:11 UTC /18:56:11 UTC+8 状态纠正

状态改为`transfer_workers_exited_pending_recovery`。无下载进程、无前处理worker，仅v6可读监控和基础设施。第二批8/14校验，未完成6文件；胃癌12%、黑色素瘤68%、CRC10%、Zenodo39%。这不是仍在后台健康下载。

退出记录：胃癌03:10:14 UTC/11:10:14 UTC+8 aria2 exit1；黑色素瘤和CRC03:24:50 UTC/11:24:50 UTC+8 aria2 exit2，黑色素瘤传输末尾明确Timeout；TNBC/前列腺癌HEAD ReadTimeout；Zenodo03:53:34 UTC/11:53:34 UTC+8重试耗尽保留partial。事件时间取自failure文件mtime，不冒充精确进程退出时刻。证据SHA与取证时间见[汇总](../../inventories/runtime/2026-10-11/download_failure_evidence_v1.json)。

v2 guarded吞并逐文件失败后仍写queue_terminal；v6仅凭日志存在显示“断点续传中”。监控刷新不是传输进度。可用空间约3.70TiB；现有证据支持超时/断流与重试耗尽，未定位唯一网络/代理责任端。恢复需新版本、端点/Range验证和真实进程状态；本次检查没有执行恢复，旧partial/失败日志全部保留。
