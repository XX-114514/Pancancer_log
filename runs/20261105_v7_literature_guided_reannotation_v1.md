# V7文献指导的逐类注释复查 v1

范围：冻结V7的4676787细胞。输入SHA、阶段终端及逐标签原因表见[证据目录](../inventories/analysis/2026-11-05/v7_literature_guided_v1/README.md)，方法见[执行合同](../methods/V7-literature-guided-reannotation-20261105.md)。本run不替代十月新图谱，也不重算冻结CNV。

## 实际阶段与命令

逻辑运行根：`${PROJECT_ROOT}/annotation_rebuild/20261105_v7_literature_guided_reannotation_v1`。实际命令通过环境变量提供根目录，不公开服务器、账户或细胞路径。

```bash
python scripts/v7_literature_guided_20261105/audit_v7_stage_a.py
python scripts/v7_literature_guided_20261105/export_mast_counts.py
Rscript scripts/v7_literature_guided_20261105/mast_initial_diagnostic_v2.R
python scripts/v7_literature_guided_20261105/mast_raw_qc_and_t6_review.py
python scripts/v7_literature_guided_20261105/review_T6_Treg_replication.py
```

这是五个独立有依赖的实际阶段，不能无条件重复；Python新输出拒绝覆盖。R由专属detached supervisor启动，PID/创建时间/真实exit及stdout留私有运行根。没有新增定时监控。

| 阶段 | 当前实测状态 | 结果边界 |
|---|---|---|
| A 冻结原因及旧panel表达 | completed，02:32:35 UTC / 10:32:35 UTC+8 | 26父类、106标签、143交叉记录；481271细胞待谱系重建；无标签改写 |
| B Mast全基因输入 | completed，02:36:28 UTC / 10:36:28 UTC+8 | 45036×95038 raw counts；逐轴精确一致；不含谱系外补招 |
| C Mast初诊 | 首轮exit1，02:38:06 UTC / 10:38:06 UTC+8 | Seurat符号转换重名，尚未聚类，失败保留 |
| C2 Mast初诊修复 | completed，02:47:11 UTC / 10:47:11 UTC+8，exit0 | 显式唯一feature ID及原名映射；HVG2000/PCA50、初诊30PC、8分辨率、RNA markers；无批次校正/状态回归 |
| D raw marker及候选审阅 | completed | raw UMI 1/2/3门槛敏感性、T4/T6全panel候选；不是双细胞或成熟亚型分类器 |
| E T6 Treg复制支持 | completed，02:45:24 UTC / 10:45:24 UTC+8 | 841/842样本支持；CellTypist Regulatory T cells 142908/156631；不整群赋Treg |
| F Mast逐群review | completed，02:52:12 UTC / 10:52:12 UTC+8 | 27探索群，38350父身份支持候选、6686竞争程序/QC复查候选；49源header可测性修正 |
| F独立算术复核v2 | PASS，02:55:52 UTC / 10:55:52 UTC+8 | 全45036细胞、57129234非零来源测量掩码、837 marker记录及27群来源计数重算；不是生物学最终验收 |

## 关键观察

Mast_0/3分别有98.8%/99.2%细胞至少两个Mast marker阳性，旧亚型未决不等于父类缺失。1408个Mast_coarse_provisional中，556个在raw≥1 UMI时同时满足Mast和T各至少两个marker；提高到每基因≥2 UMI后仍279个，应进入混合程序复查，不能仅去掉coarse后缀。尚未确认是doublet、ambient还是其他原因，不直接删除。

Tcell_6旧表中6个Treg marker、3个经典锚点均通过。旧算法先按通过基因平均分选winner，TIGIT单基因CD8_exhausted候选0.47666高于6基因Treg均分0.32886；该winner缺经典锚点后被拒，Treg候选未进入最终支持检查。应修正“整群缺锚点”的解释，并优先复查CD4/Treg成分和内部混合。

## 验收与下一步

冻结SHA不变、计数守恒、panel与原生全基因精确细胞轴、阶段输出读回和独立算术复核已完成。C2/E/F有真实成功终端，工作性review overlay已经产生；尚未发布新最终注释。下一步复查Mast的7个竞争程序/QC群，核对源特征命名/可测性并完成清理后的谱系敏感性，再进入T/NK/ILC。允许身份明确的Mast大类终点，状态与亚型分开。

初诊resolution0.4的27群不能视作27亚型：有明显来源/基因版本差异。cluster10的CD3D/E/TRBC1与Mast共存、doublet分数中位0.110；cluster15有LYZ/LST1/CSF1R及较高doublet/decontX；cluster13/25/26有间质或上皮竞争程序。均保留为待复查，不自动排除。cluster16共748细胞，Mast身份与TOP2A/BIRC5/RRM2增殖程序共存，作为cycling-enriched状态候选。

## 继续阶段：对应类型图与逐细胞复查

用户补充要求：marker和对应亚型同时展示。新版图列上方写身份/候选亚型/状态或竞争谱系，行写工作注释与细胞数；[三张v3图及panel映射](../inventories/analysis/2026-11-05/v7_literature_guided_v1/marker_targets_v3/README.md)分别展示身份/候选亚型/状态、竞争谱系QC、pan-MC文献参考panel。原图和F阶段不覆盖。

G于04:12:11 UTC / 12:12:11 UTC+8完成，60基因/1620表达记录。decontX只有指标，raw RNA没有全基因去污染层。使用1322个原样本/方法背景组校准doublet分位，44818个Mast细胞可评估。H的跨功能Mast支持和原竞争panel工作性分流为41191无强竞争信号、559联合doublet候选、828竞争程序、2365低Mast marker、93弱Mast加竞争程序；均不是最终doublet/排除判决。

H来源约束敏感性模型真实exit0于04:20:38 UTC / 12:20:38 UTC+8：36007细胞、RNA95038基因、17136个候选模型基因、2000HVG、50PC/30图PC；同池raw res0.4为24群，GSE-Harmony res0.4为15群。加权最大来源比例0.711→0.278，说明来源结构改变，不证明生物学准确率提高。RNA marker另发现B/Ig相关程序，提示原竞争panel范围不足；I已启动扩展B/Plasma/内皮/血小板及更宽T/髓系证据，不能把先前41191理解为完整排除了所有混杂。

Tcell_6的全156631细胞×95038原生基因输入已于04:28:00 UTC / 12:28:00 UTC+8完成，248354875 nnz；K已核验其精确V7 cluster_uid/major集合，52来源、17252模型候选基因，未按Treg预测筛选。通用exporter的uns scope措辞不适用于该子集，已以K精确集合核验侧车补正，矩阵不改。K的24线程全子集诊断于04:36:04 UTC / 12:36:04 UTC+8真实启动；细分未完成，不能宣称已将156631全部注释为Treg。
