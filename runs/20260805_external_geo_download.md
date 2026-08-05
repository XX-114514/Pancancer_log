# 运行记录：20260805 external GEO download and retry v3

## 运行身份

| 字段 | 值 |
| --- | --- |
| Scope | `Project_v3` 新数据下载扩展 |
| Run type | production download + versioned recovery retry |
| Initial run | `20260805_060316_slurm_download_now24` |
| Active retry | `20260806_025419_slurm_retry_v3` |
| Snapshot | 2026-08-06 05:29 +08:00 |
| Status | `in_progress` |
| Scheduler job | 26090；RUNNING；1 CPU；8 GiB RAM；14-day time limit |
| Source commit | `not_recorded`；Project_v3 路径未检测到 Git 元数据 |
| Environment | `scanpy` Python environment |

逻辑根 `${DOWNLOAD_RUN_ROOT}` 指向：
`${RUN_ROOT}/11_external_download_20260805`。

## 目标与范围

下载决策工作簿中严格标记为 `Download_now` 的 24 个 GSE：

- 下载 GEO 官方 supplementary files。
- 下载 GEO family MINiML 与 SOFT metadata。
- 只下载 SRA RunInfo metadata；FASTQ 在工作簿中为 optional，本运行不下载 FASTQ。
- multi-omic、PBMC、TCR、ATAC、VDJ、spatial 文件保留为 `sidecar_noncore`，不得混入核心 GEX。
- 按工作簿要求排除 GSE161529 的 PDX/xenograft reanalysis files。
- 不处理 `Local_no_redownload`、`Reserve_core`、`Validation_optional`、`Do_not_download_core` 行。

源工作簿 SHA256：`0fa62206e70940c62294675669aaa4e9e08e08593011e14848684d1832b1c35f`。

## 下载方法与参数

### 初始下载

| Parameter | Value |
| --- | --- |
| workers | 1 |
| retries | 6 |
| timeout | 300 s |
| chunk size | 8 MiB |
| scope | 24 `Download_now` datasets |
| manifest tasks | 122 |

由于数据盘接近容量上限，固定单流下载。下载器使用：

- GEO supplementary index 的严格 `ftp.ncbi.nlm.nih.gov` URL 过滤。
- HTTP Range 断点续传。
- retry 与指数退避。
- gzip 完整性、tar/zip 结构、文件大小和可用 MD5 检查。
- append-only events、进程锁、timestamped manifest/report。
- 已存在且验证通过的文件标记 `skipped_verified`，不覆盖。
- 失败或不完整重试使用版本化文件/attempt，不删除既有证据。

### active retry v3

| Parameter | Value |
| --- | --- |
| workers / Slurm CPUs | 1 / 1 |
| memory | 8 GiB |
| minimum request interval | 12 s |
| circuit breaker cooldown | 300 s |
| timeout | 300 s |
| chunk size | 8 MiB |
| recovery tasks | 141 |

重试启动前对 GEO 与 SRA endpoint 各执行 3 次 preflight，6/6 次返回 HTTP 200。recovery manifest 重新读取官方 supplementary indexes，因此任务数从初始固定 manifest 的 122 扩展为 141；二者分母不可直接混用。

## 初始运行结果

初始运行从 2026-08-05 06:03 至 12:16：

- 24 datasets。
- 122 tasks。
- 已知 expected bytes 65,328,877,836。
- 4 verified，118 failed。
- 主要失败证据为与 NCBI FTP 的 `Connection refused`；状态是 `complete_with_failures`，不是成功完成。

这次失败触发 versioned retry v2/v3；原 report 和 events 保留。

## active retry v3 当前进度

### 已观察事实

- Slurm job 26090 自 2026-08-06 02:54 起运行，监控快照时约 2 小时 35 分钟。
- 141 个 recovery tasks 中 74 个已有 terminal progress row，即 52.5%。
- terminal breakdown：61 `verified`、7 `skipped_verified`、6 `unavailable_upstream`。
- 已 terminal 的 verified/skipped_verified bytes 合计约 10.17 GB。
- 6 个 `unavailable_upstream` 均为 SRA RunInfo endpoint 返回空内容，属于上游无可用响应；不应写成已验证下载。
- 当前正在处理第 75 个任务：GSE274934 `GSE274934_RAW.tar`，expected size 15,310,264,320 bytes。
- 20 秒双时点采样中，该文件从 7,959,707,163 增长到 8,001,650,203 bytes，增加约 42 MB；说明传输仍在前进而非挂起。
- event stream 约每 60 秒写入 transfer heartbeat。
- 05:28 heartbeat 记录当前大文件为 9,931,030,043 bytes，约为 expected size 的 64.9%。
- 当前任务约完成 52.3%，但在 tar 完整性检查通过前不能计入 `verified`。

## 下载后的验证与数据边界

每个任务完成后按其类型执行：

- metadata text/CSV：非空、可读和必要字段检查。
- gzip：完整解压校验。
- tar/zip：容器可读与成员结构检查。
- 有 expected size 时要求精确大小一致。
- 有 expected MD5 时要求 checksum 一致。
- 只有验证通过才标记 `verified`；部分文件、网络失败或空 SRA response 保留独立状态。

下载成功不等于进入核心泛癌表达分析。每个 GSE 仍必须按照 decision manifest 中的 guardrail 区分 tumor/normal、GEX/sidecar、human/xenograft、baseline/repeated measure，并在预处理前核验 integer count matrix、样本身份和配对关系。

## 输出与证据

| Artifact | Logical path | Snapshot status |
| --- | --- | --- |
| initial decision manifest | `${DOWNLOAD_RUN_ROOT}/manifests/download_decisions_20260805_060316_slurm_download_now24.tsv` | completed |
| initial file manifest | `${DOWNLOAD_RUN_ROOT}/manifests/download_file_manifest_20260805_060316_slurm_download_now24.tsv` | 122 tasks |
| initial summary | `${DOWNLOAD_RUN_ROOT}/reports/download_summary_20260805_060316_slurm_download_now24.json` | complete_with_failures |
| recovery manifest | `${DOWNLOAD_RUN_ROOT}/manifests/recovery_manifest_20260806_025419_slurm_retry_v3.json` | 141 tasks |
| retry progress | `${DOWNLOAD_RUN_ROOT}/reports/retry_progress_20260806_025419_slurm_retry_v3.tsv` | 74 terminal rows; task 75 active |
| event stream | `${DOWNLOAD_RUN_ROOT}/state/download_events.jsonl` | active / append-only |
| retry stdout | `${DOWNLOAD_RUN_ROOT}/logs/slurm_20260806_025419_slurm_retry_v3_26090.out` | active |
| final retry summary | `${DOWNLOAD_RUN_ROOT}/reports/retry_summary_20260806_025419_slurm_retry_v3.json` | pending |

## 风险与下一步

- 数据盘在启动时记录为约 98% used；单流下载降低并发压力，但不能消除容量耗尽风险。
- 大型 tar 的“文件正在增长”只证明网络传输活跃，不等于文件完整或可用。
- SRA RunInfo 的空响应需要在最终报告中保留 `unavailable_upstream`，不得伪造 metadata。
- 下一步继续运行剩余 recovery tasks，完成每文件完整性检查，并生成最终 retry summary；随后对 core GEX 与 sidecar 做人工范围复核，再进入预处理。
