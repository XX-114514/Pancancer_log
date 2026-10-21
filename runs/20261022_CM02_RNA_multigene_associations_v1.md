# CM02固定RNA多基因关联

## 范围与前态

仅2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）原GSE286827。完整输入资格见[前置run](20261022_CM02_RNA_multigene_preflight_v1.md)。GSE288199 TableS1按用户指示缺失待补，不重复索要、不猜临床桥；三个受拒动作不重试。上一输入检查已普通本地提交8146790，公共push仍未同步。

## 固定实现与执行边界

全部32槽与先前私有设计逐项完全一致：两左父类、两端S/G2M、pre/post、全参考或原CM02成员范围。每父类至少20细胞、至少20已审病例、残差df至少10、完整设计条件数不超过1e8，连续协变量不自动删。评分排除原三marker，保留全测量稳定人类特征宇宙平均秩/maxRank1500；程序不直接赋细胞周期相位或细身份。

每病例阶段父类平均程序分数，调整双方在全部当前T中的比例、log1p细胞数、每细胞log1p原生UMI与检测基因数的父类平均，以及已审治疗/HPV。参考缺失与双文件重叠不补为无该谱系。pre/post含重复人，不能称独立复现。

2000 arm×HPV层内病例bootstrap、2000秩结局Freedman–Lane残差置换，pseudoY不重排秩；这是条件残差交换近似，不是因果识别。BH32包含不可估位置P1，缺失不当阴性。另全24父类/程序/阶段/范围与原完整basis NNLS活动的对应仅描述rho，无新P/CI或临床家族。

独立检查固定在看效果前：每CSR块首末所选行用稠密pandas秩重算；全部32主效应用另一套pandas秩/最小二乘核对；每个可估槽全部bootstrap和置换再独立计算。原完整A=W S、minmax和NNLS/BVLS/KKT、既有baseline仍需真实重放，learned/个体数组只RAM。

## 检查与当前状态

新增九项合成检查通过，含正/负方向、完整随机重放、不可估/共线、分母/参考缺失、原CSR重复索引与真实块选择。首轮本机CSR测试因macOS与Linux峰值RSS单位不同而误触预算断言，已只在合成测试提供Linux单位stub；服务器8GiB实际门槛不变、真实作业未受影响。代码与原私有设计已SHA冻结，科学效果未读取。

启动使用既有scanpy Python与scrna R；两个原公开NCBI参考仍以严格原URL/默认TLS/完整SHA从既有请求路径流送，仅RAM，不提取源文件，不读TCR字段。代码/公开程序可写自有RAM目录，个体表达/分数/病例记录不得落盘。传输观察和真实科学子进程启动分开记录；等候超时不代表失败，不因此另启重复作业。

冻结观察 2026-10-21T15:40:45.469666+00:00 UTC / 2026-10-21T23:40:45.469666+08:00 UTC+8。技术SHA见[固定实现](../inventories/analysis/2026-10-21/CM02_external_joint_resource_v1/technical_RNA_fixed32_implementation_v1.json)。真实执行尚待，源级科学结果留私有。


## 实际传输观测

观察 2026-10-21T15:44:00.334176+00:00 UTC / 2026-10-21T23:44:00.334176+08:00 UTC+8：原源完整字节本机已验，服务器接收进度真实变化。科学启动尚未观测，状态为传输中；网络等待不是统计失败，不因此改端点/端口或重跑。通用帧工具此前已测，当前bootstrap生成语法检查通过。完整仓库验证0 error/1既有warning，diff空白检查通过；后续终态另记录。


## 科学真实启动

科学子进程 2026-10-21T15:46:51.035069+00:00 UTC启动，通用实现普通本地提交669652f。公开源帧完整接收后再从RAM交给原四字段读取器；不以15:42传输开始冒充科学启动。启动收据含原代码SHA及精确PID/启动标识，科学终态及独立随机复算仍待。


## 完整终态与验收

真实科学结束 2026-10-21T15:53:13.356446+00:00 UTC；本次观察 2026-10-21T15:56:41.581668+00:00 UTC / 2026-10-21T23:56:41.581668+08:00 UTC+8。源8份白名单汇总/终态全部精确字节复制并核SHA；同PID已退出。原全部32槽、完整随机向量独立数值重放、源分母/调整设计/BH、原完整模型A及H/原baseline和阶段NNLS/BVLS/KKT均接受；源级效果与区间、分母和图表只在私有目录。本仓库只留代码、完整方法、时间线和技术SHA，不公开实际科学值。

独立源复算是数值重放，不等同新生物队列验证；pre/post和全参考/CM02范围共享病例，不能当独立重复。活动对应的24项只描述，无新P/区间或临床模型，不自动把状态同步等同原模块活性、细胞身份或直接通讯。原权威/原注释/原模型及既有临床/空间家族保持。

森林图实际查看已接受；本机matplotlib默认缓存目录无写权限后自动使用允许的临时目录，未影响图或科学结果。代码同冻结SHA未改。官方Hoch代码链本轮公开核对仅到release，具体Rmd网页cache miss未读到，临床时间语义不晋级、无大型人类对象读取；GSE288199临床表仍待用户稍后补充。

下一步是既有作者细类的状态归属/覆盖、完整原槽病例留一和有限测量敏感性先资格，再冻结实际代码与完整机会。规划留私有，未执行、不先开新临床家族。普通本地提交后公共push仍因先前自动审批拒绝而未同步；本轮不重试三个受拒动作。

复现入口：已审公共源tar仅RAM经stdin输入run_CM02_RNA_from_RAM_sources_v1.py，代码包按frozen_code_manifest.json逐文件SHA核；PANCANCER_RNA_INPUT_OUTPUT/TERMINAL_SHA256接前置资格、PANCANCER_ELIGIBILITY_SOURCE接固定原映射、PANCANCER_ORIGINAL_PHASE_SOURCE/SHA256接原R模型合同，另同前置programme/loader/network/root/RAM环境。Python=既有scanpy3.10.19、R=既有scrna4.3.3；未安装包。receiver与科学子进程均有精确PID/启动标识私有记录，当前已真实结束。通用bootstrap和本机sender为本任务私有启动包装，公开仓库提供固定receiver/scorer/aggregator/inference/auditor/projection；不能把一行stdin示例当完整基础设施安装器。

[技术完成收据](../inventories/analysis/2026-10-21/CM02_external_joint_resource_v1/technical_RNA_fixed32_complete_v1.json)。
