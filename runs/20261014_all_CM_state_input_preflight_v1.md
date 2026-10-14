# 原全172节点/1290样本状态输入预检 v1

2026-10-14T03:45:21.383967+00:00 /UTC+8 2026-10-14T11:45:21.383967+08:00冻结未执行。 [计划](../inventories/analysis/2026-10-14/all_CM_state_input_preflight_plan_v1.json)。当前T/M只覆盖参考core weight0—26.98%，本阶段回全部冻结sample/node，读canonical sample_uid/expression_v1_node_id及native units/genes/current counts/provenance，只建私有row/unit选择，不复制RNA矩阵。所有键、每单位及总cell边际必须精确一致；metadata flags不自动变推断资格或患者独立性，原172/9CM保持。

`python preflight_all_CM_state_inputs_v1.py`，PANCANCER_PREFLIGHT_ROOT/PANCANCER_ANNOTATION_ROOT/PANCANCER_CM_STATE_INPUT_PLAN。只输出新的版本目录；不读取临床结局、不用p选数据；随后独立audit/每lineage来源细胞量及native轴可行性，冻结缺谱系发现与流式全gene投影。

03:49:13—03:49:34 UTC真实exit0；52source/1290libraries/172nodes/4532933cells/35567units，所有counts键/单位行边际精确，RNA未复制；独立raw categorical/row audit待。输入currentCatalogue确为172、counts35811 rows，冻结scope过滤后35567，不把未进入1290的单位混入。

04:00:26独立sorted key/canonical row重建4532933全部行、35567单位/count/catalogue/source边际PASS。18majors完整可准备发现，原RNA未复制，推断/患者资格不提升。 [结果](../inventories/analysis/2026-10-14/all_CM_state_input_preflight_v1/RESULTS.md)。
