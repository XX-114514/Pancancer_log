# Run：CM来源内验证首轮

- Run ID：`20261010_CM_source_replication_v1`；范围：Project_v3 annotation_rebuild
- 日期：2026-10-10；本次重新核查时间2026-10-11 03:22:53 UTC
- 状态：completed_artifacts_and_audits；科学身份/临床/因果发布不据此完成
- 代码commit、准确首尾时刻和完整环境：本次未全部锁定，不能声称从公共仓库一键重现全量分析
- 目录：`${PANCANCER_ROOT}/Project_v3/annotation_rebuild/20261010_CM_source_replication_v1`

## 已观察结果

冻结1,290样本/172节点/9CM；121已审样本/45患者单位，60组成pairs/120CM边记录。13可推断候选、2争议描述性候选；357组成可估计分层、6主模型可估计分层、0 BH显著。患者/组织/癌种/分选映射不全，不能声称独立泛癌机制或临床关联成立。

## 证据与边界

[本次报告](../inventories/runtime/2026-10-11/REPORT.md)列出当前register、independent audit、final handoff和读回/图像审计；[JSON](../inventories/runtime/2026-10-11/runtime_snapshot_v2.json)保存逻辑文件、大小与SHA，患者行和完整日志留在服务器。本次没有重算该运行，也没有启动后续科学任务。

后续最终交接状态优先于生产者早期pending检查点。GSE29962340群/54,099细胞的最终身份决定独立保留。下一步补齐来源/患者与技术映射，再做敏感性和独立临床验证。
