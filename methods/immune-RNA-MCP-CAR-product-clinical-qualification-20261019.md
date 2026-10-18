# CAR输注产品：原作者三个月结局的元数据资格

## 来源与问题

在当前RNA MCP学习完成后，核已存在的公开CAR产品来源。官方CELLxGENE集合与原H5AD DOI匹配原资源论文；[PubMed](https://pubmed.ncbi.nlm.nih.gov/37738975/)将语境描述为Axi-cel输注产品、三个月PET/CT随访。[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241783)对应同一发表资源；研究层面的对应不代替精确GSM/材料桥，也不把同一资源计为两个独立来源。

## 原始字段与完整键

只读取已缓存官方H5AD obs，不读取表达或程序分数。原下载H5AD与既有原生基因归并H5AD的二进制SHA不同是正常的转换关系；需保留两种SHA/转换manifest，绑定原官方metadata SHA，并严格比较全部原细胞ID及sample_id、donor_id、response3m、3mo PET/CT、cell_type五项字段。

以当前权威图谱original_cell_id与原H5AD index精确一对一桥接全部保留细胞，逐产品、原donor、当前患者块验证互相唯一、原取样/结局字段在产品内恒定。只建立新私有旁表，不改权威Zarr、患者标签或旧临床结果。

## 结局语义

保留原作者response3m的R/NR二元编码。当前有标签产品中R对应CR，NR含PR/PD；因此这是作者三个月二元结局，不能擅自写为CR+PR作为应答的ORR，也不凭摘要推断其他结局。其余未知值不按PET片段补标、不代阴性；二元编码、原PET类别与缺失保持可追溯。

原作者cell_type与项目source_author_label须单独核实。已声明的缺失词（空、nan、None、<NA>、NA、__SOURCE_VALUE_MISSING__、unknown/Unknown）可在审核比较时规范化；其他文字差异不得忽略。CSV解析保留字面缺失词，不让解析器静默改变索引。规范化是比较规则，不修改原注释或升格为正式细胞状态。

## 独立复核与下一阶段

完整原metadata列/ID、当前精确细胞桥、所有产品结局/病例关系、作者标签来源与逐产品计数由另实现join/group核对；只有独立通过后才用于新推断设计。GSM精确桥未核时明确保留，项目已用来源不称完全未触碰外部验证。

这是制造后的输注产品领域，不是肿瘤ICI微环境。当前RNA程序在其他来源学习，CAR队列内结局关联不能自动证明细胞间通讯或临床预测能力。下一在全原RNA测量/独立接受后再固定完整新家族、单位、背景、原发现尺度与不可估计位置；本资格不产生相关、疗效效应或P。

Python/scanpy单BLAS，既有文件只读，4/2MiB私有小产物预算和500GB+1GiB保留；不新下载表达，个体表/条码、学习基因/权重和科学效果不进公共Git。
