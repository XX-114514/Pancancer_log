# 完整基因接收响应输入与先验覆盖

## 问题与范围

预声明UTC 2026-10-11T12:45:16.867675+00:00 /UTC+8 2026-10-11T20:45:16.867675+08:00。父证据v16，上一轮70模型计算验收通过但0名义及主BH显著；当前目标未完成。先确认哪些固定通讯具备可解释的靶基因先验和完整接收表达，再设计新的机制检验。

固定336记录（209LR）和回映15记录，3来源/164文库、原13子群/4家族，选定285360细胞，原生完整源809784细胞。数字不同的同名群在表达单位中分别保留；分析家族用于并行汇总，不代替原身份。[详细计划](../inventories/analysis/2026-10-11/fullgene_receiver_target_plan_v1.json)。

## 输入和方法

作者官方[先验](https://zenodo.org/records/7074291)的人类LR、ligand-target和ligand-TF三文件，固定官方大小/MD5；其中LR20KB已经下载校验，另两份正在获取。资源下载未完成不称完整先验已到库。

按原子群及文库保留完整原生基因轴，汇总原始计数、检测细胞数、完整RNA归一化均值及log均值；空单位与未测基因分别记录。先验覆盖检查全部336/15及1008跨来源组合，top200靶按先验权重和符号确定，排除LR/家族marker自关联；不根据当前表达挑靶。

[2025NicheNet协议](https://www.nature.com/articles/s41596-024-01121-9)给出从接收转录响应排序配体的路线；[作者FAQ](https://github.com/saeyslab/nichenetr/blob/master/vignettes/faq.md)明确排名不是显著性判定，PPI来源也不同于已知LR结合。当前数据只有多来源观察关系；先验不是本数据的受体特异作用证据，亦没有正式赋予细胞状态。

## 验收与状态

尚未启动科学生产或审计。必须核验原生全基因/原编号身份/样本轴及全部聚合值，独立审计通过后登记新证据。3来源已经用于发现/探索，不称独立验证。原9CM、172子群、Zarr和历史结果不变；当前没有新增统计p值。

## 2026-10-11T21:06:33.599638+08:00 资源与合成检查点

[全部事件时间线](../inventories/analysis/2026-10-11/fullgene_receiver_target_resource_timeline_v1.tsv)、[脱敏取证与SHA](../inventories/analysis/2026-10-11/fullgene_receiver_target_resource_checkpoint_v1.json)、[预声明云端核验](../inventories/analysis/2026-10-11/fullgene_receiver_target_plan_git_receipt_v1.json)。3份官方先验本机校验完整；服务器下载exit1、TF部分、靶尚未到库。分块探针通过，下一下载版本未启动。合成编码故障已修正且fixture v2通过；真实完整基因生产和独立审计均未开始。

失败部分文件、旧脚本及终态manifest保留；详情在[持续任务日志](../logs/codex/2026/2026-10-11_milestone-log-and-cloud-sync.md)。科学证据保持v16，资源准备不晋升通讯验证。运行命令将使用固定scanpy Python3.10.19与环境变量PANCANCER_PREFLIGHT_ROOT；完成资源及代码冻结后记录实际生产起止和exit。

21:07:05.928423 UTC+8（UTC13:07:05.928423）：资源恢复及失败轨迹内容d435622253a0fa575328a2b1ef39aa1c73c61b55普通push；21:07:12.642061 UTC+8独立HTTPS ls-remote/fetch核对同SHA。 [发布核验](../inventories/analysis/2026-10-11/fullgene_receiver_target_resource_git_receipt_v1.json)。

## 2026-10-11T21:11:02.854074+08:00 资源恢复v2

[固定下载合同](../inventories/analysis/2026-10-11/fullgene_receiver_target_resource_v2_spec.json)和[脚本](../scripts/fetch_server_priors_fullgene_v2.py)为资源阶段；8路1MiB精确Range，全文件大小/官方MD5确认后才能进入科学生产。真实3先验/原全基因输入未变，旧失败和部分文件保留。启动尚未发生。

## 2026-10-11T21:23:02.879031+08:00 固定代码与实际覆盖

[代码/输入冻结](../inventories/analysis/2026-10-11/fullgene_receiver_target_code_freeze_v2.json)、[资源运行取证](../inventories/analysis/2026-10-11/fullgene_receiver_target_resource_running_v2.json)。当前3来源实际节点4/6/9、合计9/目录13，原可测panel289/295/296；缺捕获家族表达覆盖NA，保留所有336/15/1008行。TF文件已完整，21680非零一般基因行不是TF目录，字段已纠正。源计数及细胞身份不变；合成检查PASS，真实计算未开始。

复现命令（固定scanpy Python3.10.19，两个根目录由环境提供）：`python fullgene_execution_worker_v1.py production`，真实终态后 `python fullgene_execution_worker_v1.py audit`。worker验证全部冻结代码SHA/原输入/父交接和官方资源终态；生产与审计分别记录新manifest/PID/起止/exit，最后手动验收科学交接。
