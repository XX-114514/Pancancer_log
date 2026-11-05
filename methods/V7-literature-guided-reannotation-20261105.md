# V7逐类注释复查：文献做法与执行合同

本轮输入是冻结V7的4676787细胞，原标签和CNV结论只读保留；新增证据与修订另存。完整注释指每类有可解释的身份层级、状态层、证据和处理结论，不预设必须消灭所有coarse/unresolved。

## 文献与可复现边界

- [Lodi 2025作者代码](https://github.com/lambrechtslab/Pan-cancer_Lodietal)：初始样本QC/DoubletFinder之后得到大类，在同谱系集合中再查看QC、竞争谱系表达和已有双细胞证据；清理后重新做谱系HVG/PCA/批次处理/图/多分辨率聚类，以RNA marker和手工核对命名。T脚本将T与NK一起复查。
- [Qian 2020正文](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf)：逐癌种和谱系复查、marker支撑与适当分辨率选择；共享与组织特异群分别评价，不能以批次混合最好作为唯一目标。
- [Chu 2023作者代码](https://github.com/Coolgenome/TCM)：T细胞CD4/CD8/非常规T/增殖分支，类型与功能程序分开；参考投射作为旁证。

以上决定分析顺序和证据结构。Lodi公开代码有未定义regression_conditions、MarkersToAdd占位和版本依赖；不能宣称逐参数复现未知部分，也不将原论文某个cluster编号、分辨率或预期双细胞比例直接套到V7。实际参数与适配要写入每次run。

## 第一步：冻结与原因归因

同时读取V7 major/final/confidence、strict eligibility、requires_recluster、原cluster、父类修订原因和原始计分/扩展panel/描述性marker表。建立每个最终标签及其原cluster的crosswalk，计数守恒。历史原因与当前新证据分开存储。

| 情况 | 处理 |
|---|---|
| rare_lineage_no_forced_subtype | 主动保留粗类；复查身份后可作为终点 |
| major_refined_subtype_requires_recluster | 新归入父类的细胞未在对应谱系重新建模；纳入新谱系分析 |
| coarse_provisional或strict_major为false | 父身份仍需复查，不能只去掉标签后缀 |
| expanded_signal_without_canonical_anchor | 旧相对表达/亚型panel条件未满足，不等同父类marker不表达 |
| top_expanded_panel_is_contamination_conflict | 查看同细胞共表达、父类关系、参考与QC；不自动定为doublet |
| marker-only或no_robust_de_novo_marker | 描述性命名/无稳定谱系内差异，复查分辨率、状态和来源 |

## 第二步：父身份、QC和测量范围

1. 先查父类marker表达比例和完整panel；共享程序或单一异源marker不直接判双细胞。
2. 按来源基因可测性计算分母。未测基因不是阴性；不把共同10432基因作为后续HVG/marker默认上限。全基因对象须精确验证与V7细胞轴一致，只读取counts及证据，不继承后加标签为答案。
3. 按来源/样本和谱系看nCount、nFeature、mt、已有doublet/ambient指标，保留原QC记录。新的排除记录原因及敏感性；稀有谱系不能按全局阈值机械剔除。
4. CellTypist保留模型、范围、置信度及来源重叠审计，作为辅助。免疫模型不能否决非免疫身份；稀有类不因参考缺失而被否定。

## 第三步：谱系内重新建模与命名

有父类调整、混群或亚型分辨不足的集合，使用全基因counts重选HVG、PCA、邻接图，检查多分辨率marker/QC/来源结构。批次变量由真实技术/来源元数据决定，记录未校正和校正敏感性；不直接回归癌种、IFN或研究目标状态。

命名同时依赖父身份、正marker、竞争marker、谱系内DEGs、跨样本/来源复现及分辨率稳定性。增殖、IFN、应激、细胞毒单列状态。小群单独检查，不靠增加分辨率强造亚型；旧阈值命中不等于发表级验收。

红系、血小板、Mast和HSPC可在身份明确时保留大类或大类+稳定marker群；HSPC若涉及肿瘤/髓系邻近程序继续审阅。ILC先核对TCR/CD3阴性证据的可测性、ILC阳性程序和NK/肥大/髓系竞争信号，不能默认其父类已成立。

## 验收与交付

每一谱系交付旧标签→原因→证据→修订标签表、marker表达/DEG、QC和来源分布、必要嵌入图及分辨率对照。细胞级overlay留服务器；公共仓库只存审阅后汇总、通用代码和脱敏日志。原标签不覆盖，终端manifest及SHA后才能称运行完成；生物学验收另列。
