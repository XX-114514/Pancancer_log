# Mast大类、候选亚型与状态的分层复核

输入是冻结V7的45036个Mast候选，全部保留；未从其他父类招募，未覆盖旧标签。H来源敏感性模型只包含36007个保守筛选细胞，res0.4得到15个探索群；图中H-编号对应这些群，不能当作15个已接受亚型。

图顶部同时标明marker对应的身份、候选亚型、状态或竞争谱系，底部保留基因名；左侧给工作注释、review标志和细胞数。点大小按来源实际可测细胞计算，颜色是各基因内部缩放，不能跨基因比较绝对表达。

- [身份、候选亚型和状态图](../../../../../assets/annotation/20261105_v7_review/Mast_H_v4/Mast_H_identity_subtype_state_v4.png)
- [竞争谱系与QC图](../../../../../assets/annotation/20261105_v7_review/Mast_H_v4/Mast_H_competing_lineage_QC_v4.png)
- [pan-MC作者参考程序图](../../../../../assets/annotation/20261105_v7_review/Mast_H_v4/Mast_H_reference_states_v4.png)
- 同名SVG可编辑；[marker与目标对应表](marker_to_type_state_map.tsv)、[15群工作注释](Mast_H15_display_annotations.tsv)、[来源可测性修正的表达](Mast_H15_expanded_measured_markers.tsv)、[RNA前10个DEGs](H_source_robust_top10_RNA_DEG.tsv)。

## 全部45036个候选的当前复核结论

| 复核路线 | 细胞数 | 当前解释 |
|---|---:|---|
| 支持保留Mast大类 | 38959 | 跨功能Mast marker支持，未触发扩展竞争规则；不要求硬分成熟亚型 |
| Mast支持且有竞争程序 | 2296 | 父身份支持保留，环境RNA/双细胞/共享程序待分辨 |
| 联合doublet候选 | 1323 | Mast与竞争程序共存且原记录背景内doublet分位≥95%，未确诊或自动排除 |
| 低Mast marker信号 | 2115 | 低于保守规则，不能由此否定Mast |
| 低Mast且有竞争程序 | 343 | 父类与质量都待复核 |

原20559个Mast_unresolved中18241在此轮支持停在Mast大类，另2318仍有明确复核原因。38959是大类支持终点，不是新增细胞或最终主分析QC纳入数。原V7字段不变，细亚型赋值数为0。[旧标签到本轮结论](old_V7_to_Mast_reviewed_transition.tsv)、[层级状态计数](Mast_hierarchical_status_counts.tsv)、[字段含义](annotation_schema.json)可追溯。

## 状态与候选亚型的证据边界

在上述38959个大类支持细胞中，非互斥RNA程序：MC_TC相关1127细胞/37来源，MC_T相关894/32来源，cycling169/29来源，IFN373/31来源，HLA-II4063/40来源。规则要求同细胞至少两个对应marker各≥2 raw UMI；低深度/未测基因不能据阴性赋相反亚型。MC_T/MC_TC参考[JCI2025](https://www.jci.org/articles/view/174981)，不声称已验证蛋白定义亚型。HLA-II程序不单独证明功能性抗原呈递。[完整程序统计](Mast_positive_state_programs_v1.tsv)。

群富集状态与逐细胞阳性数量不同：例如H-10的503个细胞是cycling-enriched群，并非503个都通过双marker门槛。pan-MC参考C0–C9来自[作者Fig1代码](https://github.com/xiaoyanger01/pan-MC/blob/main/Fig.%201.R)，不与H群编号等同，不直接迁移标签。共享CCL5/IL32/GZMB也不单独证明T污染。

## 实测与验算

I完成2026-11-05 04:39:35 UTC / 12:39:35 UTC+8；[独立验算](independent_audit.json)04:57:34 UTC / 12:57:34 UTC+8 PASS：45036个背景分位、8竞争program、1425 marker记录及2000 HVG资格均重算。层级overlay读回完成05:07:08 UTC / 13:07:08 UTC+8。父身份、亚型、状态与QC独立，细胞级对象留服务器。

H的来源约束降低来源集中度，但不证明更准确；扩展panel仍发现竞争信号，所以后续需比较来源/样本复现、环境RNA与双细胞证据。低信号或混合程序的6077个候选尚未完成最终归属。谱系外补招、全局T/NK/ILC边界及最终图谱发布仍待完成。方法和参数见[run](../../../../../runs/20261105_v7_literature_guided_reannotation_v1.md)及[方法合同](../../../../../methods/V7-literature-guided-reannotation-20261105.md)。
