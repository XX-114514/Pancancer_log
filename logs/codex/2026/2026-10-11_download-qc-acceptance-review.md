# Codex任务：下载、QC与新增验收复核

## 用户要求与起点

检查目前下载和AC情况，判断能否进入下一步；每一次任务及时push。AC本次按正式新增验收统计，未见该缩写在仓库定义；不宣称已审计其他特定AC软件。

初始本轮读取HEAD为7da773f，工作期间其他聊天正常提交到1aef8d2；其空间/家族工作和未跟踪研究综述日志均保留，不代为stage。来源基线4,676,787/54GSE；此前源库存5,341,204、QC608,684、正式新增0。本次下载线复核不修改科学身份/空间机制结论。

## 计划与规则

读取AGENTS、README/STATUS/ROADMAP/TODO、RUNTIME_SYNC、相关run/daily/incident与上次任务日志；已应用服务器pan-cancer-large-cohort-scrna-pipeline、single-sample-scrna-initial-preprocess、pan-cancer-geo-supplement-audit：原始文库QC/doublet、保留计数与来源、不推断结局，去重验收后才计新增。保留下载失败与partial，跨研究整合仍待用户决定。

## 执行与真实证据

- 10:56:11 UTC/18:56:11 UTC+8：服务器只读取证，使用既有exporter读取已知证据/当前进程，保存新的download_readiness_snapshot_v1；公开路径使用变量。源库存11/5,341,204。
- 核查accepted凭证、threshold TSV和后台进程；4份QC分片，GSE286827完成03:42:39 UTC/11:42:39 UTC+8，新增384,628；总993,312。
- 新脚本audit_qc_readiness.py部署run-local独立版本，scanpy Python3.10真实执行：19:02:09 UTC+8开始流式读取全部QC账本，只核汇总/结构/证据SHA。1,034,264原始账本行、993,312保留行、149文库；原始/保留数、文库原始数合计、X/counts形状全部一致；所有doublet computed/mt measured，保留行均predicted_singlet。输出只含汇总，拒绝覆盖旧审计。
- 四QC凭证仍pending_final_source_acceptance，clinical_sidecar未登记，正式新增0。原始计数5份/2,345,161待QC策略与来源验收。
- 无下载或前处理进程，仅monitor/infra。第二批8/14文件校验，6未完成；Zenodo39%无官方MD5完成证明。查看aria2黑色素瘤日志末尾Timeout和失败JSON；汇总最新failure时间、类型及SHA，不上传完整日志。
- RMC原始23,880细胞全指纹比对基线0、内部0；来源/文库排重仍需闭环，不以精确匹配0证明无重叠。

复现命令（脱敏）：

```bash
python scripts/collect_runtime_snapshot.py --project-root "${PANCANCER_ROOT}" --output "${NEW_SNAPSHOT}"
python scripts/audit_qc_readiness.py --run "${RUN_ROOT}" --output "${NEW_AUDIT_JSON}"
python scripts/validate_runtime_snapshot.py inventories/runtime/2026-10-11/download_readiness_snapshot_v1.json
PYTHON=python bash scripts/cloud_check.sh
```

本次exporter部署版本以runpy collect调用取证；日志/进程命令只做白名单检查，真实主机地址/路径不提交。

## 失败尝试与纠正

export完成后附加JSON探查遇到list没有items，退出1；已生成的snapshot仍完整，复制并独立校验。误探查不存在的qc/目录后依据实际worker代码定位manifests；服务器rg缺失改用find/grep，旧日志名猜测未找到后按实际文件读取。这些是取证尝试，不能记为科学流水线失败。

真正传输失败为HEAD ReadTimeout、aria2退出及Zenodo重试耗尽；旧queue_terminal只代表尝试结束。v6监控依赖文件存在产生“续传中”误报，本次以进程/失败/日志新鲜度纠正当前状态，没有重启下载或改动v6。磁盘3.70TiB，不据此断言网络/代理唯一根因。

## 修改与处理判断

新增下载/故障/全QC账本汇总、可读DOWNLOAD_READINESS、只读QC脚本；更新STATUS/README/TODO、正式run、daily、incident与同步规则。保留全部旧快照和科学记录。4份QC可进入逐来源验收，5原始计数可逐文库审查和QC；TNBC/CTCL先解来源，未完成counts只审sidecar。不授予全量整合或科学身份晋升。

## 验证与边界

服务器真实账本审计all_checked_conditions_pass=true。其证据SHA包含每份proof、ledger、threshold。仅结构读回，不重新加载全部数值矩阵；数值全读依赖既有producer证明。doublet完成不保证灵敏度，低检出结合作者历史QC解释。来源别名、文库排重、临床sidecar、final acceptance及网络唯一责任端未完成。

提交前检查与远端发布收据在下方追加。建议commit：records: audit stopped downloads and QC readiness。下一步来源验收优先，下载恢复用新版本/端点与Range验证，不重复过滤、不计metadata凑规模。

## 提交前检查结果

新runtime快照validator 0 errors；服务器全QC账本真实核对通过。cloud_check的5项合同测试、仓库0errors/1既有逻辑路径warning、release 0errors均通过；首次最后diff --check因新增EOF空行退出2，已修正文档末尾后重新检查。既有路径warning涉及annotation_rebuild@SHA逻辑引用/旧收据，未发现本次真实内网地址、患者行或凭据。未跟踪的其他聊天研究综述/方法说明/汇总均不stage。本次新增artifact索引含三份证据SHA。

修正EOF后完整cloud_check exit0；5测试/两个runtime快照校验/仓库与release/diff均通过，保留1项旧逻辑路径warning。

普通push使用已有GitHub身份临时remote-agent签名、本机Git transport，密钥10分钟过期，无私钥复制。首次ControlPath含长名字和哈希超过macOS Unix socket限制，未形成转发；改短独立控制路径重试，未切换下载代理端口。

第二次短ControlPath重试发现首轮失败留下本地签名socket，绑定失败；未删除，使用新短签名socket重试。同期其他聊天已提交36d1e06方法说明，正常保留。

## GitHub发布已核实

2026-10-11T11:07:10.755665+00:00（2026-10-11T19:07:10.755665+08:00）：普通push成功，独立HTTPS fetch/ls-remote及本地HEAD一致，完整内容SHA `c0c3d12447f3b5b718406e1411f95dcbbcf9e852`。[发布收据](../../../inventories/runtime/2026-10-11/download_qc_git_receipt_v1.json)。

共享Git index发生并发提交：本任务15文件stage后，其他聊天验证提交c0c3d12已包含这些文件及其日志两行；本次拟执行commit显示nothing added，未改写历史或重复提交。核对c0c3d12确已包含全部本任务内容后，正常push连续历史36d1e06/3cd9b8d/c0c3d12；3cd9b8d所述待同步是当时历史状态，现已一并发布。保留其他聊天未跟踪研究综述日志。本收据另正常commit/push。
