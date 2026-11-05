# V7全量doublet历史核验

用户要求确认此前是否真正去除双细胞。06:04:56 UTC / 14:04:56 UTC+8开始，只读核对冻结V7全部54来源/4676787细胞、原QC对象和账本。

原31来源按doublet规则排除5518，新23来源排除13821，合计19339；旧其中341与其他QC原因重叠，两批cap均未触发。全4676787细胞与原过滤后集合及分数精确对应，未发现旧工具阳性残留或已剔除细胞重入。28999个not_run记录是真正无分数，不是有效singlet判定。

K-09的245新候选全部原先运行过Scrublet但未判阳性，属于新增谱系复查；K-09另7个未运行，T6总786。独立复核06:16:51 UTC / 14:16:51 UTC+8 PASS，输入SHA及分数完整性06:19:50 / 14:19:50 PASS。原V7未改，没有新增模型或排除。

[完整结果](../../../inventories/analysis/2026-11-05/v7_doublet_provenance_audit_v1/README.md)、[run](../../../runs/20261105_v7_doublet_provenance_audit_v1.md)、[任务日志](../../codex/2026/2026-11-05_v7-doublet-provenance-audit.md)。本结果对应冻结V7/core，和新TME的来源/QC版本分开。
