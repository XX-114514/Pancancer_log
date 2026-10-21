# 较大独立乳腺空间来源：版本、端点和状态测量资格

## 前态与范围

前轮Hoch固定作者代码/小表资格已接受，日期按月输出的历史每日/同日误读已纠正，来源事件支持科学值私有，本地提交e6f5c10。总体目标及原2026-10-14免疫图谱CM02 v1保持；本轮尝试更大临床/空间来源，不把任何蛋白代理直接称原CM02投射或免疫治疗预测。

## 公开来源与版本

Danenberg et al 2022 Nature Genetics，DOI 10.1038/s41588-022-01041-y，原数据[记录](https://zenodo.org/records/5850952)。只读公开版本API实际接受，固定处理数据/代码v2.1.0 [记录7324285](https://zenodo.org/records/7324285)，MBTMEStrPublicProcessedDataCode.zip，910666172字节，公开MD5 f69910a57b0ea4957db12e18aa064263。未下载整个归档。图像来自v2.0.0另档；[更正mask记录15615709](https://zenodo.org/records/15615709)说明修复预处理/后处理cell标签对应，后续不得混旧mask。

原始文献/官方说明公开的病例总数只作来源事实，不作为我方有效分母。原文网页本轮重定向/PMC验证页未读到全文，不绕过、不据无匹配文字推断marker缺测。只用公开作者数据目录/面板/通用代码和表头建立资格。

## 已接受的小范围证据

ZIP目录实际接受，目录SHA 17911e718220292268a8cddfac32fb4bfa6c021a615ead0e011464c4a3658814。作者面板AbPanel.csv、类型字典SI_metaclusters.csv及五份通用R代码完整成员CRC/自身SHA接受；细胞、临床、类型表只读表头，未解析细胞/临床行或权重数组。全部请求级原域名/默认TLS/现有Clash，不改全局网络；单步4MiB预算，实际每步轻量字节见私有收据。目录/表头资格不是全归档MD5或真实数据接受。

SingleCells.csv同时有ImageNumber/ObjectNumber/metabric_id、作者cellPhenotype、组织/肿瘤/界面/血管等is_字段、蛋白字段、中心坐标及面积。公开说明明确每cell以ImageNumber+ObjectNumber唯一；from/to邻居边在同图像ObjectNumber内对应且是无向，不作发送/接收因果。ImageNumber在其他METABRIC研究不能作跨研究相同编号连接。

面板/表头共同支持CD3/CD4/CD8/FOXP3、Ki-67、PD-1、ICOS/OX40/GITR，以及Histone H3/DNA等测量背景候选。TCF7/IL7R/原RNA参考marker缺测，不得把单个蛋白当原细程序完整测量。CD45RO字段与作者代码将panel CD45重命名为cd45ro相关，尚不能据此把该字段作独立验证的memory isoform；不进入当前memory状态定义。ER别名由Rabbit secondary一行处理，需如实记录来源别名而非自作测量推断。

作者类型字典明确包含混合Treg & Tex、CD4+ T cells & APCs；不得把前者全当纯Treg、后者全当纯CD4。另一精确CD4+ T及CD8+ T标签可作候选父类，真实每cell映射/比例和同cell多marker仍待。R子群代码也对Reg使用CD4/CD8/FOXP3/PD-1及costim markers，不执行作者聚类或自行命名纯Treg。

IMCClinical.csv表头为metabric_id、ERStatus、LymphNodesOrdinal、sizeOrdinal、Grade、ERBB2_pos、yearsToStatus、DeathBreast、isValidation、PAM50、IntClust。作者SurvivalPlots.R使用Surv(yearsToStatus, DeathBreast)，属于乳腺癌死亡相关端点；不能改称ICI应答。时长起点/删失语义、真实患者唯一/事件编码/缺失、临床子集与预定义validation标志仍待具体核。未读临床行、细胞行或运行效应。

## 后续有限动作

先在原CSV小成员完整CRC/新鲜目录和字段SHA基础上，核患者唯一、随访/事件编码、预定义validation及基本临床完整性，仅元数据来源上限，不看任何细胞状态效果。符合后另冻原生细胞测量/父类/标记变化/坐标与键资格，源压缩成员仅RAM流读，不取完整ZIP或图像、不保存个体或学习数组。只有真实分母和标记/身份资格后另定完整状态/空间/临床机会，不能按效果改父类、阈值或终点。原CM02/身份/家族保持，GSE288199表待用户补、三个受拒动作不重试，公共push待同步。

观察 2026-10-21T17:39:33.605529+00:00 UTC / 2026-10-22T01:39:33.605529+08:00 UTC+8。


## 小型临床元数据资格先冻结

观察 2026-10-21T17:41:19.613051+00:00 UTC / 2026-10-22T01:41:19.613051+08:00 UTC+8。在任何新来源细胞/临床关联前，固定IMCClinical.csv小成员CRC/大小/新鲜目录，核作者患者键、重复/冲突、正有限yearsToStatus、DeathBreast及isValidation明确二进制编码和五项既有临床字段非缺失。未知不补，0时间不加偏移，不挑父类/亚型或结局。个体只在RAM，来源分母/事件/validation完整性私有，独立pandas分组及另一字典转换复算；三项合成检查通过，真实元数据分母待。时长起点/治疗背景仍待原始定义，不把作者flag自动称独立队列，不把当前metadata当原生状态捕获或原CM02投影。大型细胞表/影像/模型未读、无效应。


## 临床小表与字段语义资格实际完成

小表终态 2026-10-21T17:41:28.417220+00:00 UTC；观察 2026-10-21T17:50:13.775095+00:00 UTC / 2026-10-22T01:50:13.775095+08:00 UTC+8。完整小成员CRC、新鲜目录、患者唯一/事件/作者validation标志、基本字段完整性与独立pandas分组/另一字典编码复算接受。来源病例/事件/临床完整性及报告私有；三项合成检查通过。它们不是细胞捕获后分母、来源独立性或临床效果。

spreadsheets技能以bundled artifact-tool只读导入两份公开字段说明，未修改/导出。is_tumour是样本范围而非单细胞恶性，is_hotAggregate是成像伪影，不能作免疫热区；邻居边无向、复合键/像素标尺及已spillover校正的MFI语义记录。TCF7/IL7R缺测、Treg/Tex及CD4/APC混合类、CD45RO别名限制保留，未改权威身份。没有读取原生细胞行或运行效应。

下一步固定原版本单个CSV压缩成员的全行/CRC、复合键、伪影/样本标志、marker变异、坐标及真实患者/父类覆盖，仅RAM；源资格后另冻完整科学机会。当前大成员worker未实现/未启动，不能把预定协议当执行。普通及bundled本机Python缺SciPy已如实记技术限制，未安装包；元数据/计数资格可用现有运行时，不因此停止有效工作。原CM02/闭合家族保持，GSE288199表待补、三个受拒动作不重试、公共push待同步，goal active。


## 原生细胞资格实现冻结，真实成员尚未启动

观察 2026-10-21T18:00:14.137579+00:00 UTC / 2026-10-22T02:00:14.137579+08:00 UTC+8。新独立RAM范围读取器已实现，不修改闭合小对象传输代码：现有loopback代理请求级CONNECT、原Zenodo SNI/默认TLS，206/精确范围/长度/对象身份，1MiB四块RAM缓存；含失败请求768MiB、总2h、单请求45s、临时失败最多三试，不跟重定向、不用200全档回退、不存原始payload。

固定原v2.1.0目录/成员大小、压缩大小/CRC、表头和字典SHA。全行复合键、image单病例/样本范围、明确main-tumour且非hotAggregate，固定marker与坐标有限、面积正；补偿负MFI保留。作者混合label及通用代码明确的显示别名才映射，未知编码留待解决，不猜numericcluster。每病例父类至少20细胞及marker变异和临床键覆盖仅汇总；不划状态阈值、不算空间/临床效应、不投射CM02。

八项合成检查实际通过（缓存/seek、200拒绝、对象变更、长度/范围、预算/域名、CRC损坏、sample而非malignancy过滤及负值、混合名称）。源真实全行/CRC尚未执行，frozen证据与源码SHA独立保存；先验证/本地提交再启动。GSE288199 Table S1仍待用户补、公共push拒绝不重试。


### 分段响应validator缺失的技术纠正

观察 2026-10-21T18:02:04.615803+00:00 UTC / 2026-10-22T02:02:04.615803+08:00 UTC+8：首读及诊断重试在任何细胞行前因新增的missing_object_identity检查退出，非科学阴性。该固定公开端点206无ETag/Last-Modified；旧小对象协议允许两者均缺失。新传输明确允许缺失但逐请求tuple一致，保留原域名/default TLS/严格206范围/长度、固定原版本前后目录与全成员CRC，不能声称有服务端validator或全档MD5。诊断错误仅自身通用代码，不输出库错误内的行值；进度文件按独立run命名避免覆盖失败历史。九项合成检查通过，v2实现/SHA先冻结，待重新启动实际全量。公共push不重试。
