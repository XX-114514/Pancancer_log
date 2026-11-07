# 2026-11-06 V4 CM03：乳腺癌血管相关免疫组成方法

本阶段按用户新方向，将组成、空间组织、募集联系和治疗关联作为主线，状态分析可并行但不作前置。原CM03仍是IGHD B—PECAM1内皮—XCL1 NK及原两边；新乳腺癌代理与空间统计另立版本，不重命名原CM。

[权威完整报告](../inventories/analysis/2026-11-07/v4_breast_vascular_v1/README.md)含资源、代码算法表、全家族结果和未验证内容；[冻结配置](../inventories/analysis/2026-11-07/v4_breast_vascular_v1/frozen_breast_plan_v1.json)在本轮可估计效应前保存，原型严格测量与v2分开。

原生空间：选择作者公开6份Xenium和全18层级映射，按源文件字典序。通过HTTP精确范围读取ZIP成员；仅在完整CRC/SHA读回后可入分析，部分传输不冒充完整MD5。最大下载6GB、输出30GB、磁盘保留500GB，单CPU有限控制器，不改共享环境或旧GPU/CPDB任务。AnnData中心坐标与SpatialData边界用途分开；图像像素不可当微米。

作者Tangram三层标签为外部参考证据，不是独立真值，最大映射质量不是校准分类概率。复核原生log1p计数、原生微米与barcode覆盖后，以保守多marker检验B、T、NK、血内皮与淋巴管。ACKR1检出血内皮不直接叫HEV；CCL19/21来源分别核PVL、成纤维和淋巴管。IGHD/XCL1缺测时不硬给原亚型。空间首阶段只描述；正式患者级/区室零模型必须另冻结，并控制TLS和局部淋巴细胞密度，不能将细胞当患者重复。

I-SPY2：复用官方完整处理表达及元数据SHA，精确原始患者键；重复键全排、无模糊匹配，HER2阴性治疗前、真实pCR及HR明确，只取pembro+paclitaxel和paclitaxel对照。19,134基因/988表达列，合格248位，69/179两臂。全部marker完整才得分；患者百分位以合格全队列为分母，统一池内SD。严格panT含TRAC缺测，原24槽均不可估计；未生成可估计效应前增加CD3D/CD3E控制v2，原签名不补缺。CM03/09仍NA。

主21槽及增量3槽分别BH/BY，NA以P1留在分母、展示NA；模型、完整病例、rank/conditioning、收敛/分离、CI均显式审查。生产用本机R4.5.2 base glm未惩罚logistic，maxit100/epsilon1e-10；独立Python NumPy从原SOFT及原矩阵重新读键、明确平均并列秩、重新建全部设计，再直接似然Newton/线搜索重拟全部18设计与完整FDR。无共享环境安装。source、score、设计和患者逐行文件只私有；公共仅通用代码/汇总/图/脱敏收据。

当前18设计通过，但主/增量两家族均0 BH。对照同期资格及适应性随机化限制保留，治疗交互不称已证明的治疗特异获益。原CM相关RNA代理、组成代理与真实空间结构分开，不因阴性切换marker、来源或阈值。


## 关键组成尺度扩展

已按用户允许的大类尺度在14张审核通过的Visium上真实计算两条原CM03边及一条探索边，保留原CM结果。方法、来源汇总和Xenium通讯边界见[新增报告](../inventories/analysis/2026-11-07/v4_breast_vascular_v1/keymodule_v1/README.md)。
