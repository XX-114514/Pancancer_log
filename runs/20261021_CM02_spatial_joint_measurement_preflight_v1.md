# CM02同癌种空间联合测量资格

[CM02主线](../methods/CM02-specific-states-spatial-clinical-line-20261021.md)、[上轮状态计数诊断](20261021_CM02_fine_reference_patient_scope_v1.md)。原2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）继续优先；已关闭的168描述不重开。

## 问题与范围

GSE300147仅同癌种组织结构候选，不具备当前新辅助ICI同阶段/同终点资格。原文17切片/10HNSCC肿瘤含技术重复，非HNSCC病例需另排；实际样本/切片映射、面板测量和坐标/分割必须核验。当前不把宽Treg/CD8或某一RNA标志替代完整细状态。

## 首阶段：公共归档目录

官方tar可精确HTTP206；首512字节实测通过tar校验。脚本index_public_spatial_tar_headers_v1.py因源无ETag在首头失败，未读取常规成员载荷；原失败和代码保留。v2显式使用源可用的Last-Modified/If-Unmodified-Since、固定总大小与每响应头一致，并在结束复读首头SHA；有ETag时仍优先If-Match。不把这称整个载荷checksum验收。

每个512字节tar头以手算checksum/size和标准TarInfo双核；只允许有限PAX/长名目录数据，归档路径永不解压。最多2000成员/4000HTTP请求/4MiB目录字节，精确Content-Range和双零块结束才接受。所有表达、坐标、TCR、临床与图像载荷跳过，只写本机小型私有索引和技术收据。PANCANCER_ARCHIVE_URL只接受NCBI官方HTTPS，PANCANCER_ARCHIVE_OUTPUT是新任务私有目录。

```text
python -B scripts/index_public_spatial_tar_headers_v2.py
```

目录成功不等于数据已获取/空间身份通过；下一按真实成员类型决定可单独获取的面板/feature轴、坐标/分割和处理矩阵。实际schema和范围在任何空间效应前固定，原标签/CM保持，源级科学汇总私有。公开push仍受原拒绝后的特定许可限制，无重试。


## 2026-10-21T01:01:43.399617+00:00 UTC /2026-10-21T09:01:43.399617+08:00 UTC+8 目录完整索引与处理数据资格冻结

v1因官方源无ETag于首头停止，未读常规载荷；v2显式Last-Modified/条件请求和固定总大小、首头末尾复核替代，仅目录元数据。v2于00:58:02.045269 UTC完成：122成员、125HTTP请求/64000元数据字节，双零结束及逐头手算/标准库校验通过，所有常规载荷跳过。不能称整个归档或载荷checksum通过。

真实索引显示18份矩阵与18份细胞坐标合计113362324字节，可选择性读取，无需整体52G归档。新资格阶段事前冻结：全部18机会先核（包括后续科学排除的非HNSCC），只范围读取这两类处理文件；在服务器既有spatial环境RAM内解析，所有完整载荷/特征/细胞数组不落地。只输出预定公开marker白名单的面板存在、坐标schema/有限性与表达—坐标字面细胞键一对一匹配；不取图像、逐转录本大文件或分析TCR/临床数据。不计算空间/临床效应，也不据某一成功样本缩小全范围。原scanpy无Parquet读取器，已发现既有spatial环境含Arrow/HDF5，先核实际运行版本，不安装新软件。


## 2026-10-21T01:11:15.860466+00:00 UTC /2026-10-21T09:11:15.860466+08:00 UTC+8 实际面板/坐标初检与精确传输恢复

全18联合资格作业01:05:30.331858 UTC结束；全部机会保留。一个已终态失败定位于原range读取长度断言（69行），不是表达、坐标或生物学判定；新独立审核仅将网络读取改为固定源版本的4MiB精确分块，短响应/网络问题最多3次同块尝试，解析与身份错误不重试。原成功载荷SHA均重核，18份矩阵/坐标字面键由Counter及有限坐标另复算；原失败槽用同一完整RAM字节回放原解析器，不缩小范围或替换数据。官方Ensembl已核CXCL13稳定ID及别名，在全部面板同时查标准名、ID、别名和前缀疑似项，排除命名造成的假缺失；不借用预测或近似标志直接称该细状态。

GEO来源单位v1已明确全18归档与SOFT样本键相等、来源病例标题分组完整；诊断并非characteristics字段，原限定选择未取得。v2明确补读source_name及recurrence/p16来源类别，仅GSE级输出，源字段分类与其他研究编号无关，不读取疗效/随访。原两版源资格记录保留。


## 2026-10-21T01:21:11.863442+00:00 UTC /2026-10-21T09:21:11.863442+08:00 UTC+8 同癌种空间联合测量与来源范围完整接受

全部18处理矩阵/坐标已实际选择性读取，源字节和个体数组只RAM。原17成功载荷SHA/面板/schema复现一致；原单一短响应失败经4MiB固定源分块恢复，同完整RAM载荷回放原解析器通过。全18细胞键以Counter多重集、坐标以另一有限性实现复算；所有预设marker的标准名/官方别名/Ensembl/前缀匹配及原始特征类别已审，不将面板缺项替换成原状态。具体科学面板/数量与适配结论仅来源私有汇总，逐样本明细保持服务器。

SOFT source_name诊断与明确标题重复分组由独立MINiML XML重放；非HNSCC排除、重复切片和来源病例分组吻合。复发史、解剖取样部位与ICI时点分开；来源recurrence类别可核，但不能当免疫治疗前后时间，解剖部位仍待。p16未识别字面类别不填阴性/HPV真值，源病例也不等于新跨队列患者身份审阅。来源结构入口不替代同阶段/同终点临床验证。

完整技术验收及来源科学范围审阅结束；没有细胞状态正式分配、空间效应或新临床推断。下一问题在看到表达/空间效果前固定为实际marker共表达与捕获质量的全机会资格：多marker正证据、谱系交叠未决、非TCR人类基因计数质量背景、原发/复发与解剖部位范围、技术重复归病例。预定状态组合与敏感性仅私有计划，不能把可测的一般功能障碍标志直接命名成未测细状态。原CM02仍优先，原两版CM/正式标签不改，全部本轮任务均已终态。原拒绝公开push及随访表无重试，goal active未完成。

[完整技术验收](../inventories/analysis/2026-10-21/CM02_spatial_measurement_preflight_v1/technical_complete_acceptance.json)绑定全部实际源码/输入/审核SHA与运行版本；源级私有审阅与下一测量合同SHA分别记录。当前没有新增空间效应，临床终点验证尚待。


## 源码与重放入口

index_public_spatial_tar_headers_v1/v2.py保留严格ETag失败与明确Last-Modified替代；preflight_public_xenium_joint_measurement_v1.py在服务器现有spatial Python3.11运行。audit_public_xenium_joint_measurement_v1.py以PANCANCER_JOINT_PREFLIGHT读取原输出，重读全部选择成员并用Counter/Python有限性与官方标识核查，原失败槽用相同完整RAM字节回放原解析器。qualify_public_xenium_source_units_v1/v2.py的source_name扩展由audit_public_xenium_source_units_v1.py独立MINiML重放；review_public_xenium_measurement_scope_v1.py只读两份通过回执，完成全部来源测量范围审阅。所有实际源码、输入/审核与私有科学审阅SHA和运行版本绑定[技术回执](../inventories/analysis/2026-10-21/CM02_spatial_measurement_preflight_v1/technical_complete_acceptance.json)。下一marker测量计划已冻结但尚未执行；当前没有空间P、效应或临床增量结果。
