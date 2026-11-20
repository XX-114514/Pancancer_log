# WEHI R4跨平台综合及目标证据核对

用户持续目标：下载/验证CM12与其他模块、high/low通讯、免疫及临床指标、空间验证。上一轮全24组模型比较和上传为progress。本轮从现有代码/表/回执核对，不以旧pending文案推断进程，不重启模型。

检查：本地旧综合和builder、source_measurement/group_qualification审计、三来源患者LR完成回执、全15highlow关联、当前全24与临床/成像小表，Git隔离检出干净且HEAD为2f48812。沿用已读AGENTS/RUNTIME_SYNC与相关方法、run、日志。

执行build_CM_evidence_review_v2.py和audit_CM_goal_coverage.py均exit0：新版20个源SHA与3输出重读匹配；来源表核心身份覆盖68个来源×癌种组合，严格保留两类分别high/low>=3仅2个来源。边际计数不能当共同患者集合，后续资格来自既有三来源完整审计。未新增统计检验、重跑模型或患者映射猜测。输入文件证明的是当前分析范围，不声称穷尽外部数据。

修改：新增v2中文综合、68行来源表、六项需求证据登记、source manifest及两份代码模板；旧综合仍保留。科学结论、限制和下一步详见[报告](../../../inventories/analysis/2026-11-20/WEHI_R4_integrated_evidence_v2/README.md)。未解决：病理/分割核对、患者桥接、外部临床验证。无本轮科学运行失败。

取证记录时间：2026-11-20T01:49:45.085194+00:00；UTC+8以相同瞬间加8小时。原模型/临床执行时刻沿用原回执，本轮只做综合审查。

发布：路径/隐私、SHA、diff与仓库轻量检查后普通push，独立远端核验收据保存在本地/服务器。建议提交：Integrate all24 spatial evidence and audit CM goal coverage。

发布准备取证（UTC+8）：2026-11-20T09:49:56.548589+08:00。新增文件SHA/大小/真实路径扫描全部通过。

补充下载节点：官方Zenodo8215682页面确认942MB count/H&E/aligned-spot归档，公布MD5 a7aa24ddaa699f75c19b0c96a6a33a8a。服务器直接连接Proxy CONNECT aborted；未改代理。本机官方HEAD200，按已有SSH流式写入服务器新.part，无大文件落本地；实际已写18,857,984 bytes，传输仍进行，不能称下载或组织学验收完成。成功后需官方MD5及完整ZIP CRC，患者桥接仍未知。

仓库检查完成0错误/1既有历史路径警告；新目录路径/大小/SHA扫描及diff通过。一次索引记录命令因工作目录/Python编码失败未执行，改默认本地环境后完成；没有重启下载。
