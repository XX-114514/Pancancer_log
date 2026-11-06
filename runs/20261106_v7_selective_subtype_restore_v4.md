# V7保留既有亚型的选择性修正运行

[结果](../inventories/analysis/2026-11-06/v7_selective_subtype_restore_v4/README.md) · [方法](../methods/V7-selective-subtype-preservation-20261106.md)

本运行是用户明确纠正后的命名恢复，不重新拟合或改写RNA/QC。输入为此前完整V7当前视图及已审计v3大类诊断；环境scanpy Python3.10.19，科学代码以代码SHA登记。服务器逻辑根沿用`${PROJECT_ROOT}/runs/20261106_v7_publication_marker_annotation_v1`。

| 阶段 | UTC起止 | 退出 |
|---|---|---:|
| restore_subtypes_selectively_v4.py | 2026-11-06T12:00:24.260430+00:00 — 2026-11-06T12:02:44.593973+00:00 | 0 |
| audit_selective_restore_v4.py | 2026-11-06T12:02:44.594211+00:00 — 2026-11-06T12:04:14.565483+00:00 | 0 |
| aggregate_selective_markers_v4.py | 2026-11-06T12:04:17.801559+00:00 — 2026-11-06T12:05:30.838259+00:00 | 0 |
| audit_selective_markers_v4.py | 2026-11-06T12:05:30.838704+00:00 — 2026-11-06T12:06:00.792483+00:00 | 0 |

实际结果：4,061,939个细胞、50标签，2,843,051原名保留，1,218,888大类marker回退；606,373QC保留细胞仍在审计层。30个原非占位标签逐项登记；2,022个此前clean refit父类例外不恢复成错误旧名。完整轴、QC、原字段/CNV、全部保留/回退标签及私有TSV全流核验通过。继承原名不是新增发表验证。

UTC+8为UTC加8小时，Melbourne为UTC加11小时。实际文件SHA和终态以上述收据为准。v3的整体合并是上一轮理解错误，保留历史并明确取代其完整展示角色。
