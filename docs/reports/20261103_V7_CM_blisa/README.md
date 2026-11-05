# V7注释下的CoVarNet重算、稳健性与Xenium blisa

取证：2026-11-03 01:00 UTC至2026-11-05（UTC+8为+8 h）。运行目录为`${PC_ROOT}/Project_v5/runs/20261103_V7_CM_blisa_v1`。本页只收录节点、模块、来源和切片层面的汇总表，不含样本级H、频率矩阵或细胞行。运行细节见[运行记录](../../../runs/20261103_V7_CM_blisa_v1.md)；秩的取舍见[ADR-0004](../../../decisions/ADR-0004-v7-covarnet-rank-k9-comparison.md)。

## 结论

1. 此前没有可用的V7 CoVarNet，本次是第一次在冻结的V7注释上重算。
2. V7 K9中有5个模块在两个seed、来源等权和53次留一来源下都稳定。
3. **模块结构很大程度上是样本类型。** 与V5 CM01对应、也是最稳健的V7 CM05，在血液和骨髓样本中占活性的0.69和0.85，在肿瘤样本中只占0.07。V5 CM08对应的CM02富集于淋巴结；CM07富集于正常/癌旁组织。
4. 唯一主要出现在肿瘤组织中的是CM04（myCAF＋恶性上皮＋血管内皮＋IL1B巨噬，对应V5 CM07）。但它在留一来源下呈双峰（32/53复现、18/53丢失），在K4中也没有干净的父模块。
5. Xenium blisa只给CM04提供弱支持，且依赖成员定义。myCAF在所有panel、甚至用5,357个基因时都分不清，说明问题首先在V7的myCAF标签本身。
6. 因此，在只含肿瘤组织样本的CoVarNet出结果之前，不应对任何V7 CM做生物学解释。该运行已提交（job 30662894）。

## 1. 输入

- 冻结checkpoint `final_annotation_v7.pkl.gz`（SHA-256 `e4fda365…`，与V7 manifest一致）。在旧服务器导出逐细胞编码（[脚本](../../../scripts/v7_cm_blisa/export_v7_codes_oldserver.py)），与WEHI canonical zarr逐细胞比对：索引哈希一致，V7大类、恶性判定、strict资格均100%一致。
- 主分支：`strict_major_eligible_v7 & !exclude_from_primary_major_v2 & strict_downstream_subtype_eligible_v7`。节点用`final_annotation_with_malignancy_v7`，父类内频率后做minmax，样本至少100个细胞，节点至少出现在10个样本中。得到1,207样本、53来源、47节点（[摘要](tables/input_summary_primary_strict.json)）。
- legacy全显示分支：1,303样本、181节点。
- 程序：作者CoVarNet驱动（`scrna_r_covarnet2`），nsNMF；秩调查2–20，nrun 30；最终用两个seed（20261103/20261104），各nrun 30。

## 2. 秩

[主分支](tables/rank_survey_primary_strict.tsv)和[legacy分支](tables/rank_survey_legacy_all_display.tsv)按“第一次cophenetic下降超过0.02”的规则都会选K4。本轮仍用K9，原因见ADR-0004：

- 这是与V5（K9）的受控比较，同时改K会混淆注释效应与分辨率效应；
- cophenetic在2–20范围内差别很小且偏好小K，RSS没有拐点；
- K4太粗。

K4的4个模块跨seed相关0.98–1.00（[K9→K4对应](tables/v7_K9_vs_K4_W_corr.tsv)）。完整调查跑出的K9与单独的K9-only快速运行W逐位相同。

## 3. 稳健性（主分支K9）

[汇总表](tables/v7_primary_K9_robustness.tsv)，[逐来源表](tables/v7_primary_K9_LOSO_r.tsv)。r为按共享节点W列匹配后的相关。

| V7模块 | V5对应 | seed | 来源等权 | 留一来源 r≥0.8 | K4父模块（W r） | 主要样本类型 |
|---|---|---|---|---|---|---|
| CM05 | CM01 | 0.98 | 0.97 | 53/53 | 0.88 | 血液0.69、骨髓0.85 |
| CM06 | — | 0.99 | 0.99 | 52/53 | 0.72 | 无明显偏向 |
| CM07 | — | 0.97 | 0.80 | 52/53 | 0.72 | 正常/癌旁0.29 |
| CM03 | — | 0.96 | 0.69 | 52/53 | 0.55 | unknown类0.28 |
| CM02 | CM08 | 0.90 | 0.76 | 45/53 | 0.58 | 淋巴结0.33 |
| CM04 | CM07 | 0.92 | 0.91 | 32/53（18个<0.4） | 0.56 | 转移瘤0.19、原发0.14，血液/骨髓≈0 |
| CM01 | — | 0.87 | 0.81 | 20/53（13个<0.4） | 0.49 | 原发0.12 |
| CM08 | — | 0.23 | 0.58 | 1/53 | 0.40 | 不稳定 |
| CM09 | — | 0.66 | 0.24 | 3/53 | 0.57 | 不稳定 |

样本类型来自`20261010`运行的`sample_class`（[按类汇总](tables/CM_activity_share_by_sample_class.tsv)；1,197/1,207个样本可匹配，268个为unknown）。CM04的丢失并不集中在大来源：去掉只有3–4个样本的来源也可能让它丢失，更像NMF存在两个相互竞争的解，而不是依赖某个来源（解释）。

与V5对齐（[主分支](tables/primary_v7_v5_cm_alignment_seed_20261103.tsv)、[legacy](tables/legacy_v7_v5_cm_alignment_seed_20261103.tsv)）：V5 CM01在两个分支、两个seed中都复现；V5 CM07和CM08只是部分复现，legacy分支中CM07的对应随seed变化。

## 4. Xenium blisa

- 29个肿瘤切片，18个[单位](tables/xenium_units.tsv)。其中S16是TMA：uni/multimodal成对合并，每个单位包含多名患者。
- 参考集为V7对象中的250,898个细胞，训练与留出按样本分开。每个panel训练logistic回归，留出F1≥0.5的标签才算可识别（[逐panel F1](tables/)），每个panel 27–38/55个。
- blisa 1.0.3，CellChatDB human固定为本地副本，bin 50/25/100 µm，999次置换。每个切片可测的LR对只有7–57个（5K panel为324个，[运行汇总](tables/blisa_run_summary.tsv)）。BH显著的热点只出现在4个单位里（[LR复现表](tables/blisa_LR_recurrence_bin50.tsv)）。
- 检验方法：模块内有序发送→接收对的CCI均值，与谱系匹配的随机成员集合比较，再在各单位间做符号检验。用W top10定义成员时，CM04在两个seed下为15/18和16/18个单位为正，且在25、50、100 µm下一致；用作者网络节点定义成员时没有信号。效应小（中位z约0.55），主要来自通用的免疫/内皮LR。
- [CAF可分性检查](tables/caf_separability_check.txt)：只用成纤维细胞、5,357个基因时，myCAF_ACTA2的F1为0.54，误分到iCAF/ECM-CAF和Other。该群RGS5/MYH11最高，可能混有周细胞/血管平滑肌细胞（解释，未验证）。

## 5. 最值得做的下游

按价值和代价排序（建议，未执行的部分已注明）：

1. **只含肿瘤组织样本的CoVarNet（已提交）。** 698个原发/转移肿瘤样本、42个来源、46个节点，K9两个seed，加2–20调查。完成后重复留一来源与V5对齐。若CM04在这里稳定，才值得继续研究；若消失，说明CM04只是肿瘤与非肿瘤组织的对比。
2. **修正CAF状态定义后再判断CM04。** 用连续程序分数（收缩型ACTA2/MYH11/RGS5 vs ECM型POSTN/FAP/COL）替代离散的myCAF标签，并单独标出周细胞。看CM04依赖的是哪一种程序。
3. **只针对CM04的scRNA通讯。** 在每个样本中，对CM04成员对与丰度匹配的非成员对运行LIANA rank_aggregate，再在恶性上皮中计算接收端响应（独立靶基因，剔除LR与身份基因）。这比在全部1,300个样本上做全对LIANA更省，也更能回答具体问题。需要先在scanpy环境中补装decoupler。
4. **暂缓。** 在靶向panel上继续做blisa（每切片LR太少，myCAF不可识别）；在补齐治疗阶段和患者映射之前，暂缓临床关联。

## 6. 未完成

- 只含肿瘤样本的CoVarNet、来源等权的第二个seed仍在运行。
- 逐细胞的V5→V7细分流向、全量LIANA、接收端响应尚未做。
- GSE169609和GSE299623在V7 atlas中仍然存在（TME线已永久排除）。单独去掉它们时，CM05/CM02均保持r≥0.8；但去掉GSE169609后CM04降到0.24（双峰中的丢失一侧）。没有做同时去掉两者的版本；第1项肿瘤样本运行之后应补做。
