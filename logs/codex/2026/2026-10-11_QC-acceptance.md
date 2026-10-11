# Codex任务：四份QC全量验收

## 用户请求与起点

用户要求验收QC情况，并延续每次任务及时普通push授权。起点HEAD5f3fc18，保留其他聊天未跟踪研究综述日志。先读AGENTS、README/STATUS/ROADMAP/TODO、RUNTIME_SYNC、daily、下载run/incident、前轮QC任务与证据。

## 方法与范围

遵循服务器逐样本前处理skill、此前已读大型队列/GEO审计skill。四份分片993,312细胞；不重复过滤或重算doublet，不启动整合。此前只核ledger和Zarr结构，本次升级为全保留细胞/全部稀疏计数与原始H5AD独立比对。读原始源SHA、官方metadata SHA、完整ledger、149文库阈值/registry；重算floor/逐文库median+3MAD阈值、稳定ID、原始ID、指标、X/counts/原始计数和临床映射一致性。按5000行块扫描，单线程，输入只读、输出不可覆盖，不加载全量矩阵。

## 运行节点

2026-10-11T13:11:23.365945+00:00 /21:11:23 UTC+8启动后台audit v1，script SHA690b7b2e3d709f253a4295229249990263a0cf499cc09a63c05edb1ff6dced72，执行环境scanpy Python3.10，进度保存在${RUN_ROOT}/logs/independent_qc_acceptance_20261011_v1.log。真实结果尚待终态；不把启动作为验收成功。

## 原研究核查与边界

三GEO83个RNA文库精确GSM；CAR-T59文库公开H5AD患者产品关联。发现GSE292046是T细胞富集snRNA-seq，8文库来自2患者的前/后和重复，不能算8独立患者或无偏全组织组成；保留技术/富集字段。CAR-T原研究别名为GSE241783，59患者输注产品，不能视为原位肿瘤生态。另检查GSE294109为人类RNA队列而非同论文动物/空间部分。

官方/原研究来源：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292046 、GSE241783同站、https://pmc.ncbi.nlm.nih.gov/articles/PMC12629806/ 、https://doi.org/10.1016/j.ccell.2026.05.007 。GEO/PMC open部分受recaptcha、出版社open失败；使用搜索返回的官方记录/原研究数据声明，不引用第三方方法推断。未取得的作者doublet方法细节保持unknown。

原proof中clinical_sidecar_sha256=null是external TSV未使用，clinical_source明确为official_published_H5AD_obs；不能仅凭null认定缺少临床来源。本轮核对实际registry与obs，不推断患者结局。

## 待终态补录

验收结果、异常、输出SHA、技术/正式来源资格、验证、提交与独立远端核验随后补录。建议commit：records: accept full-cell QC audit with source qualifications。

节点检查：cloud_check exit0，5测试通过、仓库0errors/1既有逻辑路径warning、release0errors；diff检查通过。原研究官方URL和脱敏context已保存；真实数值审计仍运行。源码在目标环境编译通过，启动后真实进度与资源检查正常。
