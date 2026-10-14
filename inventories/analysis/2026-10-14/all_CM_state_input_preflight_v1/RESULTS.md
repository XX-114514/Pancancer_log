# 全部原CoVarNet状态输入预检：已核所有行

03:49:13—03:49:34 UTC producer exit0，04:00:26独立audit完成。原172节点/1290样本、52sources、4532933cells、35567 sample-node units；所有私有canonical global row选择、单位key/count和原major/subcluster一一匹配，独立sorted-key membership重建全部453万行通过。[审计](independent_audit.json)、[source范围](all_source_native_scope.tsv)、[18major/source范围](all_source_major_native_scope.tsv)。原Counts35811行中scope外样本单位未混入35567。

Major范围包括T1391840、Mac398575、B482248、Monocyte125704、Myeloid57302、T/NK214073、Lymphoid128717、Dendritic16266、Mast42396、Neutrophil56628，以及Epi1032027、Fib147160、Endo122858、Mural69292、Melanocytic97127、Glial93663、Mesenchymal50002、Hematopoietic7055；总数精确4532933。所有major已存在原数据，初期pilot缺口是输入范围，并非原数据没有这些cells。

原生gene轴16497—58336，median29580.5；不能把不同轴不在签名的genes填实测0。只保存私有行/单位索引及继承provenance，不复制任何全量RNA矩阵。原身份/编号/W/H/9CM不改，1290库不称1290独立患者；分析资格仍另外检查，不因metadata行匹配就全部可推断。

下一按固定source/major/block均衡准备有限发现cells，保留所有18major（single-source群标来源限定）；限定训练matrix和流式原生投影，检查500GiB余量，不复制453万全gene对象。Spectra过程先验与cNMF无先验对照/既有programme候选分开，全部CM成员/edge映射还需实际state scores；不以本row preflight完成完整goal。
