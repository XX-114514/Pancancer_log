# 既有NSCLC治疗后临床队列：冻结程序测量入口

当前终态：全基因测量、独立原矩阵与固定程序投影验收完成，后续完整病例级临床家族也已执行验收，见[临床run](20261016_NSCLC_frozen_signature_clinical_v1.md)。下方前半保留05:34的历史入口判断；09:37后的章节记录后续实际完成，不能以早期pending覆盖新终态。

## 问题与实际范围

承接HNSCC完整自动家族未晋级后的临床/队列方向，先核既有GSE243013输入。原始GEO对应治疗后NSCLC免疫专用队列，适合问治疗后生态/病理疗效关联。不能直接当治疗前预测；病例独立性、治疗/终点与当前端点捕获尚待另审。

已有全原生病例/作者状态提取缓存，及先前CD8—巨噬问题的固定200例选择集。其权重聚合和来源可测性已有旧独立审计，本轮没有重新分析疗效或复制患者行到本机/公共仓库。该来源不在当前免疫图谱来源合同，不能假定已有与原172或新图谱的权威身份映射。

## 方法和真实时间线

- 2026-10-16 05:31:17 UTC（13:31:17 UTC+8）：原CRC两端MCP1及新DIALOGUE全部56端点槽，共58固定计划位置的缓存基因/绝对权重覆盖预检完成，不按内部留出P或临床效果筛程序。
- 05:34:11 UTC（13:34:11 UTC+8）：逐基因集合与绝对权重标量求和独立复算PASS，未提取与原生缺失严格分开。病例/身份及临床关联未验收。
- 05:35:24 UTC（13:35:24 UTC+8）前实际核对既有count_matrix_download_v1完整压缩计数文件存在；历史完整性证明下一入场重核，本轮没有新下载。

缓存有真正每细胞log1pCP10K后的均值，可在充分测量/固定标准化时线性聚合程序贡献；仅均值CP10K取log不能替代。现有提取目标来自先前LR/footprint集合，尚不足当前完整/原95%权重门槛。这是提取范围限制，不能据此判断基因未表达或原矩阵无基因。

完整私有覆盖表`${PREFLIGHT_ROOT}/NSCLC_existing_frozen_programme_measurement_20261016_v1/all58_private_fixed_programme_cache_coverage.tsv`，SHA93b58c0ee6a1ae3a010ef93c850e8ff58f5d47b5bef9688cb20f060e12faaba7。原native提取审计SHA1c65ae928bc193353b277bfef47482172e6f981ced20bbb82fe2252700288967；原固定选择/家族审计SHAa25c659408a0f58f6958392346cdbbc96460a053b064a2c5081ac96e72836a16。独立覆盖审计状态PASS_ALL58_FIXED_GENE_WEIGHT_CACHE_COVERAGE_AND_UNEXTRACTED_NOT_ABSENT。

## 下一动作与解释边界

先重核计数文件、原条码/基因轴及源指纹，从既有原矩阵补提取全部冻结程序所需基因及质量统计；与旧缓存交集逐项一致、全RNA分母保留。无新数据下载，不降95%门槛或对未提取基因补零。只从基因测量覆盖不能推断相应细胞已捕获。

下一审核全病例临床/当前端点，保留旧200选择集作其原问题的对照；不能将CD8捕获门槛自动套给Treg/其他端点。作者免疫状态/粗谱系与原172细群明确区别，不能由名字或单marker直接对号；缺上皮等端点保持不可估计。治疗后关联、来源内复现、跨来源一致和真正独立验证分别报告。

主要效果前冻结完整临床家族、患者/病例块及必需丰度/质量/治疗/病理协变量；保留既往本来源临床探索事实，不称项目级独立验证。此run交付的是可追溯测量差距与可执行入口，尚未运行补提取/程序评分/疗效检验，也没有发现新轴。完整goal active。

## 代码、失败与发布

preflight_NSCLC_existing_frozen_programmes_v1.py及audit_NSCLC_existing_frozen_programmes_v1.py均真实执行完成；此前较大来源文件清单SSH观察超时，缩小摘要并压缩连接后取证恢复，无科学重复提交。纯读取代码/文件定位的多次本机未匹配glob已改为rg路径筛选，不改科学数据。

代码轻量compile、仓库/隐私/diff检查；只发布方法、执行轨迹和SHA，不发布患者行、科学统计或学习签名。与HNSCC修正节点同普通commit/push、独立HTTPS核云SHA。下一方法见[主方法](../methods/automatic-receiver-prior-and-signed-attribution-v1.md)。

发布检查实际完成：本轮6份新增/修正源码内存compile全部PASS；仓库验证0错误、1项既有路径警告；工作区和暂存diff检查通过。仅本任务文件及共享文件自有段落纳入，保留并发工作。普通push后的内容SHA与独立云核验收据另存本地工作区。

## 2026-10-16 09:02—09:29 UTC /17:02—17:29 UTC+8 全基因入口实际推进

08:54后重新取证既有完整矩阵和内存、真实文件系统空闲及组配额，满足本轮全基因单位数组和文件系统预留门槛。09:02:54启动原矩阵流式汇总，09:24:13 UTC真实终止并写出全部原生基因的整数计数、阳性细胞数、每细胞CP10K及log1pCP10K均值；完整原RNA分母、压缩矩阵SHA/CRC和坐标唯一顺序通过。旧缓存全部提供基因、全部单元的三个表达口径逐项相等，原作者标签/病例选择未改。原生检出基因统计与作者n_genes_by_counts一致；另一个n_genes字段存在差异，质量字段保留真实原矩阵口径，不混用。

独立验收09:19:15启动并等待原生产真实句柄，随后按预先固定、不读疗效的64单元复读原矩阵，覆盖每个作者细胞群全部原生基因；同时核全部文件SHA、所有单元守恒及全部旧缓存交集。09:28取证进程仍live，尚未据生产终态宣称独立全验收。它独立重放64单元，不冒称独立复读全部单元。

09:21:21 UTC全部作者病例临床字段及全部细胞类别捕获审核完成；保留原始字面治疗/结局和缺失，病例内冲突另表登记，无临床效应。病例编号仍不等于已审独立患者；连续病理比例的语义未凭列名指定。09:29:18启动全部58冻结槽的完整原生轴覆盖重核，仍不做程序分数或疗效检验。

私有目录：`${PREFLIGHT_ROOT}/NSCLC_full_native_gene_units_20261016_v1`（生产终态）、`NSCLC_full_native_gene_unit_audit_20261016_v1`（独立审计）、`NSCLC_all_case_clinical_scope_20261016_v1`（临床字段）、`NSCLC_full_native_programme_coverage_20261016_v1`（覆盖）。完整矩阵SHA `c96d918548ed857ee69c1a1f1e40669d072c07ba6d7a757b6c6051386ee7968e`；全基因log均值数组SHA `8f17f4920524cacdcc709a525d2dfda42d9c45f4ff9da0d518b745a4af3ce402`；单位/质量表SHA `a032ef126f10ce20a327eefbddad43ea9fcb565cdb6b1587ef7917846d1728f9`。原生生产源码SHA `1bc875fe88f43b99e083f2408b664257b1520cda900a2cfbbb9f3b06587259c9`，launch_receipt保存实际PID/startticks。

下一需核原上下文的major分母、训练质量均值/SD、冻结逐基因校正与作者细胞群的谱系相容性及外推范围。作者Treg/巨噬/CD8群不直接冒称原172细身份；缺上皮不得转投免疫细胞。临床全家族/处理及结局含义须在效果前明确。[GEO原记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013)及[原论文](https://pubmed.ncbi.nlm.nih.gov/40147443/)确认治疗后取材与免疫范围，不能当治疗前预测。

## 09:37—09:51 UTC /17:37—17:51 UTC+8 真实验收及临床运行

09:37:24 UTC，独立完整文件/所有旧交集和预固定64单元全基因复读通过，状态PASS_ALL_OUTPUTS_OLD_OVERLAP_AND_PRESPECIFIED64_RAW_UNITS_ALL_NATIVE_GENES；终态SHA `bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1`。明确并未独立复读全部单元，但生产全细胞RNA/坐标守恒均已验收。

09:36:38全部支持上下文的原RDS均值/SD/权重/逐基因五项校正导出与加权系数一致；全部原56槽含未支持项和来源可测部分分区校准保留。09:39:46向明确相容作者谱系的固定投影完成；09:44:24独立全部56范围、47,628项分数/谱系分母/训练范围标志通过，终态SHA `be340aaf5e0f69e6b0f8bc6a44ea143ffda24287d61c2c76021416a2fe6aed8a`。独立审计v1再次因index/column同名排序歧义停止，v2只改sort_index，原分数不改，旧失败保留。全source域检查表明亚群粒度/组成外推需说明；不据投影存在宣称细身份或状态已验证。

09:48:53在任何临床效应前写出完整合同及所有预定设计，09:48:55预检终止；09:50:05正式唯一全家族运行，09:51:08完成，独立设计/HC3/完整BH/全部LOO及预固定三例完整wild/bootstrap审计真实启动后待验收。原/残差口径完整控制后的推断不变性是强制项；不能因评分形式选择有利结果。该临床阶段的主要检验是wild bootstrap，不是原CM的残差置换。

方法见[预冻结临床家族](../methods/NSCLC-frozen-signature-clinical-family-v1.md)。私有`${PREFLIGHT_ROOT}/NSCLC_context_frozen_parameters_20261016_v1`（参数/校准）、`NSCLC_context_native_programme_projection_20261016_v1`（分数/独立审计）、`NSCLC_frozen_signature_clinical_preflight_20261016_v1`（临床合同/病例资格/设计）、`NSCLC_frozen_signature_clinical_results_20261016_v1`（完整统计/审计）。设计表SHA `4a07a989f0ba4ecb7a807455e720fbb53c535c92576ce3cede62ee6f8c698236`；临床合同SHA `60a0146fd102fd1057005ee27012fd6be95bb0d6ae3e746dac8df4e536754a4c`；完整临床表SHA `6ab08b1f21a6c4188e52c1f1487adcfbc0ec8526ee1753e1369429b604c767da`。科学结果不复制公共仓库。当前整体goal未完成。
