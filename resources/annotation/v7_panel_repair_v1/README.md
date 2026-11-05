# V7 panel repair v1

新版包含112个证据模块、409个不同基因，并逐条处理旧162项。它是供谱系内重新注释使用的证据注册表，不是已校准的自动分类器。旧V7标签保持原样。

- [全部112项及全部marker](ALL_PANELS.md)：直接阅读。
- [完整Excel工作簿](V7_panel_repair_v1.xlsx)：112项panel、162项迁移、634条逐基因依据及21项来源范围，4张表均可筛选。
- [筛选用完整表](all_new_panels.tsv)：正向、降低方向、核心组件、命名marker、鉴别对照、适用背景和来源。
- [旧162项逐条迁移](old162_to_new.tsv)：迁移规则，不是迁移细胞标签。
- [机器读取的权威registry](registry.json)：22个parent、43个branch、37个state、5个isotype、4个boundary、1个作者完整功能签名。
- [逐基因依据](gene_reference_trace.tsv)：634条模块—基因记录，其中266条在本次已核验作者材料中有出现证据；368条仍明确标记为项目整理、待进一步逐基因核实。出现证据也不等于作者使用该组合进行分类。
- [21项来源及适用范围](references.tsv)：包含论文、作者代码和官方方法文档，不是21个独立验证数据集。
- [方法与算法边界](../../../methods/V7-panel-repair-literature-workflow-20261105.md)。

每个模块的核心组件分别检查，不要求每个细胞表达全部marker；“组件中存在可测基因”“所有核心基因在全部细胞来源可测”和“同一细胞共表达”是不同字段。降低方向必须使用可测的匹配对照，未测或掉失不能当阴性。

`author_exact_signature`仅表示作者完整签名的保留，不表示身份分类器。`curated_literature_informed`为本项目依据论文整理的组合。`curated_reference_gap`仍用于Oligodendrocyte、Neuron及Epi_renal_tubular_like：已有相关参考背景，但尚缺本轮逐项核验的合适原始marker依据。后者只覆盖近端小管方向，不能泛指全部肾小管。

红系、HSPC、ILC、Mast等在大类证据通过后，可以正式停在大类；是否再细分取决于数据和研究问题。marker命名需要该marker本身的可测性及相对富集证据，不能靠同panel其他共享基因替代。
