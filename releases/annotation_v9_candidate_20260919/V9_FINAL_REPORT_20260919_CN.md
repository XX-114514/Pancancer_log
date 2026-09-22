# V9 泛癌注释：全流程结果报告（第二版，2026-09-19）

> 本版取代 `V9_FINAL_REPORT_20260918_CN.md`。与第一版相比，L2 层修正了两处问题（见第 3 节 IMPL-19、IMPL-20），并在修正后的 L2 基础上重跑了 S4–S7。
> 对象：4,676,787 个细胞（54 个 GSE、1,322 个样本、43 个癌种）。未通过的门槛按原值报告，并给出原因，不做事后调参。
> 技术细节与逐次运行记录见 `V9_IMPLEMENTATION_AND_RUN_LOG_20260917_CN.md`；运行谱系见 `runs/.../13_v9_annotation/RUN_HISTORY.tsv`。

---

## 0. 与第一版相比的变化

| 项目 | 第一版（int_v3） | 第二版（int_v4） |
|---|---|---|
| L2 作者标签宏平均 F1（预注册口径） | 0.740 | 0.753 |
| L2 作者 F1（剔除 AML 骨髓 GSE116256） | 0.893 | **0.933** |
| T_NK_ILC / Mast 逐类 F1 | 0.884 / 0.753 | **0.981 / 0.904** |
| 分选纯 T 数据集 GSE179994 中 T_NK_ILC 比例 | 约 52%（45% 未解析，1.46 万被误标为其他谱系） | **97.6%** |
| 留一 GSE：宏平均 F1 < 0.6 的 GSE | 4 个 | **0 个**（加权 F1 中位数 0.955） |
| SingleR 对照总体一致率 | 0.812 | **0.836** |
| HSPC / 红系 L2 细胞数 | 5,233 / 4,501（**绝大多数为假阳性**） | 45 / 130 |

第一版 L2 中的 HSPC、红系、粒细胞以及 `signature_reassigned` 细胞**不应再使用**。

---

## 1. 最终注释结果

### 1.1 L1 大区（`S2_L1/20260918_l1_v7`，与第一版相同）

| L1 | 细胞数 | 占比 |
|---|---|---|
| Immune | 2,907,799 | 62.2% |
| Epithelial | 1,004,360 | 21.5% |
| Stromal | 261,433 | 5.6% |
| Endothelial | 94,438 | 2.0% |
| Neural_crest_glial | 81,574 | 1.7% |
| unresolved | 306,138 | 6.5% |
| 排除样本 | 21,045 | 0.4% |

作者标签宏平均 F1 0.944（非盲）；双种子 ARI 0.972；未解析 + doublet 6.9%。

### 1.2 L2 大类（`S3_L2/*/20260919_l2_v13` → `S3_L2_integrate/20260919_int_v4`）

| L2 | 细胞数 | | L2 | 细胞数 |
|---|---|---|---|---|
| T_NK_ILC | 1,636,131 | | Glial | 44,225 |
| Epithelial | 999,530 | | Mast | 38,939 |
| Myeloid | 634,108 | | Granulocyte | 11,763 |
| B_Plasma | 491,191 | | Neuronal | 8,311 |
| Fibroblast | 210,005 | | Melanocytic | 4,424 |
| Blood_EC | 91,480 | | Lymphatic_EC | 4,095 |
| Pericyte_SMC | 50,985 | | Erythroid_Megakaryocyte / HSPC / Mesothelial | 130 / 45 / 10 |

未定型：`Immune_unresolved` 82,058、`Neural_crest_glial_unresolved` 22,501、`Stromal_unresolved` 2,201、`Endothelial_unresolved` 1,562、`not_annotated_at_L2` 343,093（L1 未解析或样本被排除）。
置信度：high 2,885,906、medium 914,855、low 325,126、none 207,807。
跨大区回流：355,019 个细胞进入回流，其中 246,697 个被重判成功。

各大区稳定性（双种子 ARI）：免疫 0.978、基质 0.793、内皮 0.703、神经嵴 0.596。

### 1.3 L3 细胞类型（试点，`S4_L3/*/20260919_l3_v3`）

**B_Plasma（370,290 个细胞，分辨率 0.2，5 个 cluster；bootstrap ARI 0.904，双种子 ARI 0.693）**

| 名称 | 细胞数 | 置信度 |
|---|---|---|
| Plasma | 108,277 | high（AUROC 0.992） |
| B_naive | 40,395 | high |
| Plasmablast | 1,895 | high（Plasma + cycling z ≥ 2） |
| B_Plasma_novel_c0_TNFRSF13B（记忆 B：TNFRSF13B/SCIMP/AIM2） | 122,065 | 待复核（B_memory AUROC 0.690 < 0.80） |
| B_Plasma_novel_c3_MYBL1（增殖/生发中心样 B） | 22,548 | 待复核（与 B_naive 差值不足 0.05） |
| B_Plasma_novel_c4_IL7R（**T 细胞污染**：IL7R、CD3D/E/G、CD2） | 3,412 | 待复核 |
| B_Plasma（父类保留） | 71,698 | low |

**Myeloid（536,799 个细胞；所有分辨率的 bootstrap ARI 均为 0.56–0.61，未达 0.80 门槛，按预注册规则取 ARI 最高的分辨率 2.0，21 个 cluster；双种子 ARI 0.501）**

| 名称 | 细胞数 | | 名称 | 细胞数 |
|---|---|---|---|---|
| Mono_classical | 32,613 | | Mono_nonclassical | 7,137 |
| Macro_C1Q | 31,338 | | mregDC_LAMP3 | 5,708 |
| cDC2 | 17,667 | | pDC / cDC1 | 4,152 / 2,958 |

另有 12 个 cluster 待复核（CCL18⁺ / APOC1⁺ 巨噬、S100A12⁺ 单核、增殖髓系 SPC25/UBE2C、IFN 应答 RSAD2、SELENOP⁺ 组织驻留巨噬等），共 287,833 个细胞；父类 `Myeloid` 135,807 个（low）。

**解读**：髓系在任何分辨率下都没有稳定的离散划分（bootstrap ARI ≤ 0.61），这与髓系细胞在肿瘤中呈连续状态谱的文献认识一致；因此 Myeloid 的 L3 结果只应作为候选，**不宜在复核前用于下游定量**。

复核队列（共 15 个 cluster，证据已填好）：`preregistration/manual_review_decisions_L3_pilot_TEMPLATE_v2.tsv`。

### 1.4 L4 状态（`S5_states/20260919_s5_v3`）

- 15 个状态程序中 13 个至少在一个谱系通过共表达一致性：IFN 应答、MHC-II 抗原提呈、S 期、G2M 期、解离应激在 12 个谱系通过，缺氧在 9 个；耗竭、细胞毒、组织驻留、SPP1-TAM、myCAF、iCAF、tip 各在其限定谱系内通过；IgA/IgG 可用基因不足 3 个，无法打分。
- 连续 z 分数与 `state_high`（z ≥ 2 且该谱系通过一致性）：`objects/v9_state_z.pkl.gz`、`v9_state_high.pkl.gz`。
- NMF 元程序：B_Plasma 44 个（327 个样本）、Myeloid 54 个（544 个样本）。**数量偏多，含较多只覆盖 3 个 GSE 的小程序，聚类阈值需复核后再用于下游。**

### 1.5 恶性轴（`S6_malignancy/20260919_s6_v3`，未重跑 CNV）

| 判定 | 细胞数 |
|---|---|
| not_candidate | 3,380,652 |
| malignant | 437,320 |
| uncertain | 396,625 |
| excluded_qc | 192,352 |
| not_evaluated_role_changed | 149,197 |
| cnv_no_malignancy_detected | 81,922 |
| not_evaluable | 38,719 |

3,353,123 个细胞直接继承 V7 CNV 结论；113,599 个细胞在 V9 中变为 candidate 但 V7 未评估；49,686 个细胞存在 reference/candidate 角色互换。**632 个样本列入 CNV 重跑建议（D-C1），未执行**，清单：`tables/rerun_suggestions_D-C1_not_executed.tsv`。

---

## 2. 验证结果（S7）

| 检查 | 结果 | 门槛 | 结论 |
|---|---|---|---|
| 作者标签 L1 宏平均 F1 | 0.944 | ≥0.90 | ✅（非盲） |
| 作者标签 L2 宏平均 F1 | **0.753**；剔除 GSE116256 后 **0.933**；剔除 HSPC/红系两类后 0.932 | ≥0.85 | ❌（预注册口径） |
| L2 逐类 F1 | Epithelial 0.982、T_NK_ILC 0.981、B_Plasma 0.964、Myeloid 0.963、Fibroblast 0.957、Glial 0.926、Mast 0.904、Pericyte 0.748；**HSPC 0.002、红系 0.101** | ≥0.70 | 两类不达标 |
| doublet 比例（各大区） | ≤0.9% | ≤5% | ✅ |
| L1 未解析 + doublet | 6.9% | ≤5% | ❌（逐 GSE 解释） |
| 双种子 ARI | 免疫 L2 0.978；基质 0.793、内皮 0.703、神经嵴 0.596；B_Plasma L3 0.693；Myeloid L3 0.501 | ≥0.85 | 仅免疫 L2 通过 |
| 留一 GSE CellTypist（L2 中位 F1） | B_Plasma 0.996、Mast 0.994、T_NK_ILC 0.994、Epithelial 0.991、Myeloid 0.988、Blood_EC 0.974、Lymphatic_EC 0.957、Fibroblast 0.929、Pericyte 0.898；Glial 0.649；红系、粒细胞 0.070；Melanocytic 0.005、Neuronal 0.0 | ≥0.80 | 9 个主要类别通过；**没有任何 GSE 低于 0.6**（加权 F1 中位数 0.955，最低 GSE138665 0.466） |
| SingleR benchmark（HPCA，20 万细胞） | 总体 0.836；T_NK_ILC 0.996、Blood_EC 0.991、B_Plasma 0.972、Granulocyte 0.929、Myeloid 0.885、Epithelial 0.882、红系 0.800 | 仅对照 | HPCA 主标签不含肥大、淋巴管内皮、黑色素、间皮 |
| 环境 RNA 敏感性（原始 vs decontX） | 免疫 0.996、基质 0.981、神经嵴 0.980、内皮 0.975 | ≥0.90 | ✅ |

**L2 门槛为什么仍未通过**：作者标签中的 HSPC 与红系只来自 GSE116256（AML 骨髓）。修正后 V9 不再把其他数据集里的零散细胞误标为 HSPC/红系，但在 AML 骨髓中，免疫大区本身以白血病祖细胞为主，V9 无法把它们与髓系可靠地区分，所以这两类的 F1 接近 0。剔除该 GSE 后，其余 8 类的宏平均 F1 为 0.933。

---

## 3. 过程中修正的方法学缺陷（共 11 项）

前 9 项见第一版报告与实施记录第 3–5 节。本版新增：

10. **IMPL-18 的假阳性（由我此前加入的规则造成）**：IMPL-18 允许"没过自己 cluster 的标签、但只超过另一个标签"的细胞直接改判。q99 阈值按定义会放过约 1% 的阴性细胞，IMPL-18 把这 1% 尾部变成了标签。第一版中全库**所有** HSPC（5,233）与红系（4,501）细胞都来自这一路径；在分选纯 T 数据 GSE179994 中，14,577 个 T 细胞被标成了 B/髓系/红系/肥大/HSPC/粒细胞。修正（**IMPL-19**）：改判只有在 ≥50% 的 kNN 邻居也超过同一标签时才保留。免疫大区改判数从 42,762 降到 4,248。
11. **逐 GSE 阈值不适用于分选数据（IMPL-20）**：逐 GSE 校准假设每个 GSE 内都有本地阴性细胞。分选数据（GSE179994 纯 T、GSE139249 纯 NK 等）不满足。修正：某 GSE 在该大区内 ≥80% 细胞属于同一个第一轮标签时，改用全库合并阈值。免疫大区中有 5 个 GSE 触发：GSE139249、GSE152048、GSE179994、GSE189889、GSE299340。

另有两个工程问题：4 个 S3 任务并行时，默认线程池超过了用户进程数上限（4,096），已统一限制为每库 8 线程；S7 的细胞筛选只接受 `confirmed|reassigned`，把上皮和救援细胞排除在验证之外，已放宽。

---

## 4. 已知限制

1. **L1 的作者 F1 非盲**（修正案 01 在观测后作出）；L2 的 IMPL-17–20 同样是在看到结果后加入的。独立证据是留一 GSE 与 SingleR 对照，二者都显示主要类别表现良好。
2. **AML 骨髓（GSE116256）**：HSPC、红系在 L2 基本无法注释（全库仅 45 和 130 个细胞）。
3. **稳定性**：基质、内皮、神经嵴的 L2 以及两个 L3 试点的双种子 ARI 低于 0.85；髓系 L3 没有任何分辨率通过 bootstrap ARI 门槛。
4. **L3 仅完成 B_Plasma 与 Myeloid 两个试点**，且 15 个 cluster 需人工书面结论后才能冻结。
5. **神经嵴大区**以 GBM 与葡萄膜黑色素瘤的恶性细胞为主，L2 大多为 mixed/回流，由 S6 处理。
6. **独立投票在 2M 老队列上较弱**（scATOMIC 只覆盖新队列）；报告中一律写"与独立分类器的一致性"，不写"准确率"。

---

## 5. 需要你决定的事项

| # | 事项 | 说明 |
|---|---|---|
| 1 | **L3 复核队列（15 个 cluster）** | `manual_review_decisions_L3_pilot_TEMPLATE_v2.tsv`。我的判断：B c0→B_memory、c3→增殖/生发中心 B、c4→T 细胞污染；髓系 c0/c1→CCL18⁺ 巨噬（Macro_C1Q/FOLR2 连续谱）、c2→APOC1⁺ 脂质相关巨噬、c4→Mono_classical、c9/c12→增殖髓系、c13→IFN 应答（状态）、c17→Macro_FOLR2_LYVE1 |
| 2 | **632 个样本的 CNV 重跑（D-C1）** | 批准后约 15 小时 |
| 3 | **低稳定性子集是否改用 scVI 表征** | 内皮、神经嵴 L2，髓系 L3 |
| 4 | **L3 匹配门槛是否在 V10 中调整** | 不在本轮事后放宽 |
| 5 | **其余未批准重算** | decontX 全量、scATOMIC 补 2M 队列、文献参考图谱 CellTypist、外部验证数据集 |

---

## 6. 产物位置

| 内容 | 路径（相对项目根） |
|---|---|
| L1 + L2 逐细胞注释 | `runs/.../13_v9_annotation/S3_L2_integrate/20260919_int_v4/objects/v9_L1_L2_annotation.pkl.gz` |
| L3 试点 | `.../S4_L3/{B_Plasma,Myeloid}/20260919_l3_v3/objects/v9_L3_annotation.pkl.gz` |
| 状态 | `.../S5_states/20260919_s5_v3/objects/v9_state_{z,high}.pkl.gz` |
| 恶性轴 | `.../S6_malignancy/20260919_s6_v3/objects/v9_malignancy.pkl.gz` |
| 验证 | `.../S7_validation/20260919_s7_logo_v4/`、`.../20260919_s7_singler_v3/` |
| 冻结清单 | `v9_annotation/preregistration/V9_FREEZE_MANIFEST_v2.json` |
| 复核模板 | `v9_annotation/preregistration/manual_review_decisions_L3_pilot_TEMPLATE_v2.tsv` |
