# 优先免疫CM的TCGA与治疗队列功能代理探索

## 对象和边界

2026-10-14免疫图谱CM（v1，T_whole_TME_primary_tumor）CM02/CM04。此前依据原W、网络和原生表达命名的两个T模块，分别研究调节相关/细胞毒，以及记忆/GZMK效应的功能读数。原W/S/H、注释和发现输入均只读；原10月16日HNSCC30槽阴性家族不重开。功能代理不是原CM精确投影，独立数值重拟不等于独立患者临床验证。

## 固定测量

[主要效果前合同](../inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/frozen_plan_before_effects.json)列出34个唯一基因、7个panel、全部槽和门槛。基因按队列内（TCGA癌种内）平均秩百分位，多基因等权均值，再取两臂min；所需基因100%测得且可变，否则不放宽签名。pan-T控制由CD3D/CD3E/TRAC组成。没有新训练、优化阈值、按效果替换基因或二分连续得分。

- pCM02：FOXP3/IL2RA/CTLA4/IKZF2/TNFRSF18与GZMB/PRF1/CTSW/NKG7/CD8A/CD8B。
- pCM04：TCF7/LEF1/CCR7/IL7R/SELL/LTB与GZMK/CTSW/CCL5/GZMA。
- 原生T节点表达从接受的全基因充分统计读取，不抽样细胞，不将表达展示当DEG或同细胞证据。

## 患者、终点和设计

TCGA预设33癌种，实际01型原发实体瘤读取涉及32癌种；LAML按固定材料规则无患者，仍保留槽。合并重复aliquot、冲突患者隔离；未知事件不当删失。PFI/OS分开，TCGA-CDR癌种/随访限制保留；TCGA不当ICI队列或随机未治疗对照。基础Cox加pan-T、年龄、性别；分期敏感性仅用明确整体I–IV（病理优先、临床回退），不以T/N/M冒充分期，不将未知填0。样本>=60、事件>=max(20,10p)、分期覆盖>=70%、满秩/条件数门槛按合同。

既有5个ICI来源按真实患者/临床key对接：Riaz只pre，Hugo排除on并平均重复治疗前患者，Rose只精确交集且保留时点限制。客观应答为CR/PR与SD/PD，NE/未知缺失，生存事件必须明确0/1。主要模型joint+pan-T；双臂增量另用joint+armA+armB+pan-T，并应用更严格样本/少数类门槛。增量是min给定两臂后的平衡/非线性问题，不是配体受体检验。

## 完整检验家族

TCGA 264（33×2×2终点×2模型），既有ICI应答10、生存20、增量10。新增SU2C在下载/效果前独立冻结[8槽合同](../inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/SU2C_external_plan_before_download_or_effects.json)：应答2、生存4、增量2。每个完整家族计算BH/BY，失败槽以1保留在分母、结果展示NA；不将不同家族q当合并312槽q。LRT P与Wald CI分开。

## 公开数据补充

SU2C-MARK原研究：[Nature Genetics](https://www.nature.com/articles/s41588-023-01355-5)，[作者Zenodo v3](https://zenodo.org/records/11179623)。393全队列不是RNA数；152 RNA按作者QC Keep保留136、Flag16排除。TableS1 header=2，按精确RNA/患者key连接，事件/天数使用作者明确字段。只下载公开TPM/临床成员，各成员ZIP CRC与SHA通过，未核完整ZIP MD5。原完整TPM本地保留，服务器输入为效果前固定34基因/全部152RNA精确子集及作者完整临床，源SHA和子集SHA绑定；不改变评分数学。受控dbGaP phs002822.v1.p1未申请，不代填PI/伦理或联系作者。

## 数值、图和敏感性验收

主162个可拟合设计全部以另一logistic/Cox Newton实现和显式Efron风险集似然/梯度/信息重拟，全部家族BH/BY核对。临床源字段和映射另查，未声称由独立采集患者验证。全150个Cox有描述性残差/时间诊断（不是正式PH证明）。既有10个IFN/增殖控制槽和SU2C8个年龄/性别/线数/联合治疗控制槽均保留不可估计项，不新设P家族、不替换主分析；后者内部数值核对不冒称主独立162设计范围。

4图保留全部33癌种和6治疗来源：PFI、OS、全ICI终点、16原生T节点全panel表达。v1因NA不进入自动范围导致行对齐风险，v2固定行边界；CI截断加箭头，平衡图扩展轴。全部v2实际查看，底层312槽TSV与v1逐字节相同。原版本不删除。

## 资源与输出

本轮读取既有bulk大矩阵，仅产生小panel、设计/审计/图表；100MB本轮服务器输出上限与500GB留空合同，不写data1或复建原Zarr。最终报告和科学系数仅授权工作区；公共仓库保存合同、代码、来源与技术收据。实际运行详见[run](../runs/20261018_priority_immune_CM_bulk_v1.md)。后续先同终点独立复现、跨研究患者重叠、临床/纯度与精确细胞身份测量，再研究治疗交互；没有非ICI比较不能确证治疗特异获益。
