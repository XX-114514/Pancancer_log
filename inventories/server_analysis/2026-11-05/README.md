# 2026-11-05服务器分析证据快照

[分析报告](../../../docs/reports/20261105_server_analysis_audit/README.md)组织科学解释；本目录保存可核对的聚合材料。取证时间见`snapshot.json`，补充证据各行另有时间。原始事件时间只采用原run记录；文件mtime不是执行或完成时间。

| 文件/目录 | 含义 |
|---|---|
| `source_evidence.tsv` | 206项allowlist证据，原始SHA、大小、mtime、取证时间及公开副本hash |
| `supplementary_evidence.tsv` | 69项补证，包括配置、原W的hash及迁移问题列表；不公开W本体 |
| `imported_code.tsv` | 43个代码原件和公开副本SHA，以及路径/保护修改 |
| `analysis_lines.tsv`、`version_crosswalk.tsv` | 16工作线和不可互换的数据/CM版本 |
| `state_corrections.tsv` | 9项早期pending状态或解释限制的更新 |
| `tables/` | V5/V7聚合转移、V7 rank/执行摘要、跨单位空间汇总、NMF匹配、去批次和早期QC |
| `annotation_pilots.json`、`annotation_preflight.json` | 真实单来源候选试点及历史预检摘要，不是最终身份 |
| `annotation_current_presence.json` | 72侧表存在检查；未逐细胞重验，不证明全量分流完成 |
| `spatial_completion.tsv`、`nmf_tasks.tsv` | bin摘要存在及45个NMF运行收据摘要；不替代退出码和科学验收 |
| `nmf_scale_sensitivity.tsv` | 本次读取已有W做的尺度敏感性核查，没有重新拟合 |
| `migration_receipt_summary.tsv`、`data_presence.tsv` | 15份旧迁移问题表均无问题行；指定对象存在情况，本次未重hash大对象 |
| `run_directory_counts.tsv`、`recent_run_directories.tsv` | 1,179个运行目录的名称日期汇总及近期列表，不能据此判断在哪台服务器执行 |

公开副本将物理根替换为逻辑变量、统一LF，并将TSV末列空缺写作`NA`以保留列数及缺测语义；source SHA仍针对原字节。hash_only/NA表示只发布hash或无适用值，不能解释为零。NMF源码commit收据沿用原3列无表头格式（项目、commit、源URL），没有完整环境导出。

离线核验：

```bash
python scripts/server_audit/verify_snapshot.py
```

该入口核对64个公开证据引用、43代码hash、45任务、4,676,787细胞转移总数、路由计数、rank覆盖及符号检验算术。它验证记录一致性，不能验证细胞身份或机制。主要采集入口为`scripts/server_audit/collect_snapshot.py`；补充取证操作及环境见任务日志，复采须使用新目录。
