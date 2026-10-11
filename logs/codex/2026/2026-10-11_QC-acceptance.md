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

## 全量技术终态与确定的映射问题

2026-10-11T13:15:27.440923+00:00 /21:15:27 UTC+8：四份全部36条件PASS，共993,312细胞；逐保留行全部非零计数对原始来源一致，X/counts一致，源/官方SHA一致，稳定ID与原始顺序/临床registry一致。各总counts：CAR-T4716938125、头颈癌3332069517、脑772248028、肾842459456。

13:13:04 UTC /21:13:04 UTC+8 context审计：CAR-T59患者、HNS29已关联患者+1未关联库7463QC细胞、脑2患者、肾15患者。doublet删除103/289/1352/186；低检出不能证明残留不存在。GSE292046有富集/重复设计；CAR-T产品单列。没有把临床TSV缺失当作官方来源缺失。

启动节点1417788bc00bbf0cfe7c704eb4dfba65dcace772普通push，并独立HTTPS ls-remote/fetch核对SHA一致。随后定位头颈癌官方GEX标题post2导致旧patient解析遗漏；同患者pre/post已指向唯一patient_uid。建立只修patient/timepoint的新registry和stable-ID overlay v2，结局/层级保持不变、原QC/旧registry保留。实际修复/独立核对结果随后追加；正式来源验收仍需研究别名/基线文库审查，不因技术PASS自动增加。

## 元數據修復與登記

13:18:17.691405 UTC /21:18:17 UTC+8：GSE286827 patientfix_v2完成，依据官方标题token及同患者pre/post唯一source映射，保留post2字面code。新registry和稳定ID overlay只留服务器；原文件SHA保留。独立审查用split/token而非producer regex，7463ID一一精确对应，所有结局/tier等原字段未改动，未映射库变0，患者仍29。13:18:52.908477 UTC独立PASS，13:18:52.909137 UTC技术验收登记。相关JSON保存汇总与SHA，无患者或细胞行。正式source acceptance仍0。

## 变更与验收判断

加入QC_ACCEPTANCE可读报告、全数值审计/技术登记、修复/独立覆核JSON和版本化修复代码；更新README/STATUS/TODO、run/daily及证据索引。993312技术验收通过，可来源内前处理；576512肿瘤组织与416800产品分开。富集snRNA/两患者重复和低doublet检出边界保留。原研究别名GSE241783已查明，但54基线原始文库/别名完整对照未完成，不晋升正式新增或整合。

未验证：Scrublet敏感性与原作者doublet完整方法、全部基线研究别名/文库对照、唯一患者结局、跨研究整合。没有删除数据、重跑过滤或改动科学身份。Web open遇到recaptcha/出版社失败，采用官方搜索记录及本地已冻结metadata，未把无结果视为无信息。

终态发布前：完整cloud_check exit0，5合同测试/仓库0errors(1旧逻辑路径warning)/release0errors/diff均通过，所有新JSON人工核对只有聚合和逻辑路径。原研究方法open未取得者明确unknown；真实全量计数和元数据独立复核已通过。专用签名转发在空闲后自动退出，使用新短socket恢复，不换下载代理端口。
