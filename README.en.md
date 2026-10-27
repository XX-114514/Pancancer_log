# PanCancer

[简体中文](README.md) | [**English**](README.en.md)

## A pan-cancer atlas connecting cellular composition with T-cell states and spatial organization

**Use traceable cell identities to identify co-occurring cellular combinations and interpret them in defined tissue and clinical contexts.**

The tumor microenvironment comprises multiple cell populations. Understanding this ecosystem requires answering three questions: which cells vary together, which states they occupy, and how they meet within tissue. PanCancer brings a pan-cancer composition atlas, an immune single-cell atlas, clinical cohorts and spatial protein imaging into a stratified analytical framework. Its focus is a **CTLA4–GZMB-associated T-cell co-variation core**, connecting compositional structure, state organization and clinical associations.

The current immune atlas covers **2,670,057 cells, nine sources and 400 physical libraries**, with 18 major classes, 270 descriptive clusters and five continuous T-cell state axes. Within the primary-tumor T-cell analysis scope, all seven CM02 core members remain among the top 10 members after fitting with equal source weights. A fixed functional RNA proxy is associated with treatment response and overall survival in two immunotherapy cohorts, respectively. In regional lymph nodes from non-small-cell lung cancer (NSCLC), a regulatory-associated CD4-like protein state shows local spatial enrichment with three CD8-like states.

[Atlas and versions](#atlas-resources-and-versions) · [CM02 core](#a-t-cell-co-variation-core-retained-after-source-balancing) · [Clinical associations](#a-fixed-functional-rna-proxy-connects-with-immunotherapy-outcomes) · [Spatial organization](#nsclc-lymph-nodes-reveal-defined-t-cell-state-organization) · [Methods and evidence](#methods-reproducibility-and-evidence)

### Key findings

| Scientific claim | Strongest quantitative evidence | Question addressed |
|---|---|---|
| An immune atlas connects identity with state | 2.67 million cells; 18 major classes; 270 descriptive clusters; 16 T-cell clusters and five state axes | Which cells can be compared within a shared coordinate system? |
| The CM02 composition core is retained after source balancing | Effective-basis cosine similarity approximately 0.979; all seven core members retained among the top 10 | Does the cellular combination depend on source weighting? |
| A fixed regulatory-associated–cytotoxic RNA co-enrichment proxy has clinical associations | IMvigor210 response OR 3.57; SU2C_MARK overall-survival HR 0.51 | Does this functional theme connect with patient outcomes? |
| NSCLC lymph nodes contain spatially localized state relationships | Local enrichment factor 1.38 for FOXP3/CTLA4-associated CD4-like and PD1/TIM3-associated CD8-like states | Do the associated states form neighborhood structures within tissue? |

These layers retain their respective measurement objects: **cell-composition CMs, an RNA functional proxy, relative protein states and clinical endpoints**. Each is presented below with its own denominator and statistical definition.

## An atlas retaining tissue context and cell states

![Immune atlas overview showing sources, major cell classes and a shared embedding](assets/publication/20261026/atlas_overview.png)

**Figure 1 | Tissue context of the immune single-cell atlas.** The original analysis figure shows the atlas and its composition across sources. The atlas contains 2,670,057 cells and 400 physical libraries; cells, libraries and patients are distinct counting units. Sources include tissue, immune-enriched material and CAR products, with module discovery organized by the corresponding material scope. Panel d describes native endpoint eligibility at the immune-v1 freeze; subsequent clinical metadata bridges are recorded under their respective new versions. [Atlas delivery and complete source table](inventories/runtime/2026-10-14/immune_atlas_v1/RESULTS.md).

The atlas combines the following features:

- **Preserve the full tissue context.** Immune, stromal, vascular, epithelial and other tissue components provide the background for interpreting cellular ecosystems.
- **Separate identity from state.** A population can simultaneously express memory, cytotoxic, interferon or proliferation-associated features. Continuous state axes complement cluster labels in describing these changes.
- **Retain sources and analysis units.** Primary tissue, metastatic tissue, blood, immune-enriched material and CAR products define distinct analytical contexts. Repeated material from a patient is handled according to the relevant design.

![Five continuous T-cell RNA state axes and cluster-level state profiles](assets/publication/20261026/T_continuous_states.png)

**Figure 2 | Continuous states beyond T-cell identity.** Five fixed RNA scoring axes describe naive/memory, cytotoxic, dysfunction-associated, interferon and cycling-related expression; the lower-right panel summarizes T-cell cluster state profiles. The scores and the figure label “dysfunction” denote associated RNA programs; functional activity is assessed through the relevant functional study design. This figure provides descriptive identity–state coordinates. The clinical and spatial analyses below use their separately frozen RNA-proxy and protein-state definitions.

## A T-cell co-variation core retained after source balancing

### CM02 places regulatory-associated, effector and memory-related T-cell components in a shared composition context

The full version identifier is **2026-10-14 immune-atlas CM02, v1, `T_whole_TME_primary_tumor`**. Its display name is **CTLA4–GZMB T-cell co-variation**.

**Scope definition:** `whole_TME` retains the capture-scope label of this historical version. Lineage enrichment must be reviewed separately using the original library-preparation evidence. This T-cell CM uses proportions within the T-cell major class and describes T-cell composition under the corresponding capture conditions.

| Frozen core member | Display marker | Biological interpretation |
|---|---|---|
| `T_CTLA4_235` | CTLA4 | Regulatory/activation-associated T-cell component |
| `T_GZMB_231` | GZMB | Cytotoxic effector-associated component |
| `T_GZMK_232` | GZMK | GZMK-associated effector/memory component |
| `T_LTB_237` | LTB | LTB-associated lymphocyte component |
| `T_IL7R_239` | IL7R | IL7R-associated memory/persistence component |
| `T_CTSW_240` | CTSW | Effector-associated component |
| `T_LEF1_243` | LEF1 | LEF1-associated naive/memory component |

Markers and these interpretations help describe the compositional theme. Formal fine-grained identities combine multiple markers with source evidence. **The biological focus of CM02 is the joint organization of its components.**

### All seven core members are retained with equal source weights

![CM02 core members and source-balanced stability of all four T-cell modules](assets/publication/20261026/cm02_source_balance.png)

**Figure 3 | Testing the composition core by changing source weights.** With inputs, preprocessing and factorization rank frozen, the original unit-weighted objective is compared with an equal-source-weight objective using four fixed seeds. CM02 effective-basis similarity is 0.979126–0.979140, raw W similarity is 0.949640–0.949663, nine of the top 10 members overlap, and all seven original core members are retained. The right panel shows all four components in this scope. The arrangement on the left indicates membership and does not encode network edges or loading magnitude.

**This experiment tests the stability of the composition core to source weighting.** The T-cell combination retains a similar structure when the relative contribution of sources changes, providing a stable anchor for context-specific ecological analysis. The four seeds share the same data; transfer across sources is evaluated using separate cohort designs.

[Source-balancing methods](methods/CM-source-balanced-sensitivity-20261024.md) · [Statistical definitions and comparison scope](docs/PROJECT_PUBLICATION_EVIDENCE_20261026.md#来源平衡的完整比较)

## A fixed functional RNA proxy connects with immunotherapy outcomes

The composition finding motivates a testable functional question: **is joint elevation of regulatory-associated and cytotoxic expression associated with clinical outcomes?**

To examine this theme across cohorts, predefined multigene panels separately score regulatory-associated and cytotoxic expression. The smaller of the two scores defines the co-enrichment proxy, and the principal model includes a pan-T background measure. Panels, directions and readouts were frozen before calculation of the new effects; outcomes were not used to select genes or optimize cutpoints.

![Clinical associations of the fixed CM02 functional proxy in IMvigor210 and SU2C_MARK](assets/publication/20261026/cm02_clinical_associations.png)

**Figure 4 | Associations of a fixed RNA proxy with two clinical endpoints.** The left panel reports the treatment-response odds ratio and the right panel the overall-survival hazard ratio, retaining BH correction within the original respective testing families. Each effect corresponds to a one-standard-deviation increase in the fixed co-enrichment score.

| Cohort | Endpoint and model | Analysis units / events | Principal effect and 95% interval | BH q |
|---|---|---|---|---|
| IMvigor210 | Treatment response; logistic regression | 298 / 68 | **OR 3.57 [1.80–7.06]** | **0.00144** |
| SU2C_MARK | Overall survival; Cox model | 127 / 60 | **HR 0.51 [0.32–0.82]** | **0.0202** |

**Clinical takeaway: the fixed regulatory-associated–cytotoxic co-enrichment theme is associated with higher odds of response and a lower mortality hazard in the two respective cohorts.** The cohorts address different endpoints. The RNA proxy measures a functional theme, while the composition CM retains its original cell-proportion definition.

This section evaluates **clinical association**. Predictive performance relative to either single arm, pan-T background and joint models is a separate incremental-value question. Complete predefined comparisons are retained in the [companion evidence notes](docs/PROJECT_PUBLICATION_EVIDENCE_20261026.md); association effects are reported using their defined statistical interpretation.

[Frozen clinical plan and panels](inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/frozen_plan_before_effects.json) · [Testing families and comparison definitions](docs/PROJECT_PUBLICATION_EVIDENCE_20261026.md#临床代理的固定定义)

## NSCLC lymph nodes reveal defined T-cell state organization

The clinical proxy describes a patient-level functional theme. Spatial analysis asks a further question: **how are the associated states arranged in tissue?**

In surgically sampled regional lymph nodes from NSCLC, the analysis anchors on a candidate CD4-like state with relatively high FOXP3/CTLA4 scores and compares three CD8-like protein states. Fixed checks account for MKI67 proliferation strata and local composition, neighborhood degree and spatial blocks. Effects are aggregated by registered patient unit.

![Spatial associations and local enrichment for three T-cell state relationships in NSCLC lymph nodes](assets/publication/20261026/nsclc_state_organization.png)

**Figure 5 | Three state relationships within one tissue context.** The analysis includes 11 registered patient units and 97 images, with all three predefined receiver endpoints displayed. The left panel shows spatial odds ratios stratified by MKI67; the right panel shows centered enrichment factors relative to the local conditional null model. These are distinct estimands.

| CD8-like receiver state | MKI67-stratified spatial OR [95% interval] | Local conditional enrichment factor [95% interval] | Patient units with positive local effects |
|---|---|---|---|
| PD1/TIM3-associated | **2.54 [2.05–3.39]** | **1.38 [1.27–1.51]** | **11/11** |
| MKI67-associated | **1.46 [1.24–1.68]** | **1.18 [1.07–1.21]** | **11/11** |
| CCR7/CD127-associated | **2.00 [1.33–2.33]** | **1.17 [1.03–1.31]** | **9/11** |

BH q is approximately 0.00209 for the three principal associations and 0.003 for the local conditional tests; corresponding BY q values are all below 0.01. Overall effects remain positive in leave-one-patient-out analyses. The shared permutation-resolution floor is not used to rank the three receiver states by strength.

**Spatial takeaway: regulatory-associated CD4-like states form tissue neighborhoods with multiple CD8-like states, and the enrichment remains supported after proliferation and local-background checks.** States are defined by relative protein measurements. The panel lacks a CD3 channel, so candidate populations are described as CD4-like/CD8-like. The same-source secondary family follows the initial spatial observations and provides state-organization evidence within this NSCLC context.

### Interpreting relationships in their tissue context

An external melanoma tumor-draining lymph-node (TDLN) source has completed all 15 primary testing slots under the frozen relative-protein-state and spatial framework, with independent numerical auditing. This source provides a **cross-cancer tissue-context comparison**. Primary slots, measurable denominators and structural missingness are retained in the companion evidence notes. Effects are estimated separately for each tissue context rather than pooled into a universal Treg–CD8 mechanism.

This perspective makes the next question specific: **which state combination, in which cancer, treatment stage and tissue region, has reproducible ecological and clinical relevance?**

[Spatial methods](methods/CM02-LN-regulatory-specificity-20261022.md) · [Complete spatial families and external context comparison](docs/PROJECT_PUBLICATION_EVIDENCE_20261026.md#空间组织的完整情境对照)

## Atlas resources and versions

Versioned resource layers connect composition discovery, immune themes and identity expansion across the full tissue background.

| Resource | Core scale | Scientific use | Authoritative entry point |
|---|---|---|---|
| **October 2026 pan-cancer CMs, v20261010** | 172 subclusters, 1290 libraries, nine CMs | Pan-cancer cell composition and source/cancer-context relationships | [Within-source and communication analysis](runs/20261012_CM_edge_coverage_v1.md) |
| **2026-10-14 immune atlas, v1** | 2,670,057 cells; 400 physical libraries; nine sources | Immune themes, T-cell states and material-stratified CM discovery | [Formal atlas delivery](inventories/runtime/2026-10-14/immune_atlas_v1/RESULTS.md) |
| **2026-10-14 immune-atlas CM02, v1** | Seven-member core in the primary-tumor T-cell scope | CTLA4–GZMB-associated T-cell ecology focus | [Scope-specific module names](inventories/analysis/2026-10-18/immune_CM_naming_v1/CM_NAMES_CN.md) |
| **Full-context native archive expansion** | 9,439,479 observations; 65 source keys | Unified identities, within-source cluster markers and cross-version correspondence | [Hierarchical annotation methods](methods/full-context-atlas-hierarchical-annotation-plan-20261023.md) |

Source keys, native observations, independent studies and patients are recorded separately in the full-context expansion. Server checks on 2026-10-26 confirmed completion of the full projection and within-subset marker production for 13,279 subsets defined by fixed rules. Coordinate precision, independent audits, unified identities, crosswalks and use-specific views undergo separate acceptance checks; the complete acceptance register defines the final unified-identity release. The expanded archive and the formal immune-v1 results retain their own version scopes.

### Other composition themes in the same T-cell scope

| 2026-10-14 immune atlas v1; `T_whole_TME_primary_tumor` | Descriptive name | Ecological interpretation |
|---|---|---|
| CM01 | DNAJB1–CCL4 stress/chemotaxis-associated T-cell co-variation | Stress- and chemotaxis-associated composition |
| **CM02** | **CTLA4–GZMB T-cell co-variation** | **Joint organization of regulatory-associated, effector and memory components** |
| CM03 | XIST–LTB/LEF1 mixed T-cell co-variation | A mixed theme interpreted with identity and source context |
| CM04 | GZMK–IL7R–CTSW T-cell co-variation | Effector/memory-associated composition |

The [complete scope-specific descriptive names](inventories/analysis/2026-10-18/immune_CM_naming_v1/CM_NAMES_CN.md) provide the full catalog. Factor instances in overlapping material scopes are interpreted within their respective scopes.

## From composition themes to state ecology

The next scientific steps address one question: **resolve co-infiltration into specific states, tissue relationships and incremental clinical value.**

1. **Cell identities provide stable coordinates.** Combine within-source clustering, expression evidence from each cluster, references and author-provided lineage information to release traceable identities and cross-version mappings. Record cluster-derived markers separately from evidence borrowed from a parent population.
2. **Distributional readouts characterize active states.** Compare the measurement stability of means, active-cell fractions and quantiles in existing native expression data, and examine whether small populations of highly active cells carry patient-level coordination signals.
3. **Tissue and clinical contexts guide ecological interpretation.** Within a defined cancer, material and treatment stage, use patient units to compare state combinations, spatial positions and matched clinical endpoints, and assess incremental value over single-arm readouts and general T-cell background.
4. **Communication analysis serves localized state relationships.** For candidates supported by state and tissue evidence, organize ligand–receptor evidence using sender expression, receiver-cell pathway/transcription-factor responses and spatial compatibility.

These questions form successive analytical layers: composition discovery nominates themes, state measurements define the objects, spatial validation locates them, and clinical designs evaluate their use.

## Methods reproducibility and evidence

### What each main figure establishes

| Main figure | Argument | Supporting evidence |
|---|---|---|
| Figure 1: Atlas overview | A unified cell atlas retains source and tissue context | Formal immune-v1 counts, coordinates and source acceptance |
| Figure 2: Continuous states | RNA state dimensions complement identity labels | Fixed T-cell state scores and cluster profiles |
| Figure 3: Core stability | The CM02 composition core is retained when source contributions change | Original/source-balanced weights; all four components and four seeds |
| Figure 4: Clinical associations | A fixed functional theme connects with cohort-specific clinical endpoints | Original clinical-family correction and independent numerical audits |
| Figure 5: Spatial organization | Defined protein states form neighborhoods in NSCLC lymph nodes | Proliferation strata, local conditional background and patient-unit resampling |

### Reproducibility resources

- [Homepage evidence, complete comparisons and figure notes](docs/PROJECT_PUBLICATION_EVIDENCE_20261026.md)
- [Measurement definitions and full comparison scopes for source-level claims](docs/PROJECT_PUBLICATION_EVIDENCE_20261026.md#各层证据回答不同问题)
- [Figure provenance and SHA manifest](assets/publication/20261026/figure_manifest.json)
- [Figure-rendering script using aggregate statistics](scripts/render_project_homepage_figures_v1.py)
- [Atlas analysis interface](docs/IMMUNE_ATLAS_V1_HANDOFF.md) and [all 270 descriptive clusters](inventories/runtime/2026-10-14/immune_atlas_v1/subcluster_catalog_aggregate_public.tsv)
- [Current status](STATUS.md), [scientific roadmap](ROADMAP.md), [task board](TODO.md), [formal run index](inventories/runs.tsv) and [complete run history](runs/README.md)
- [Historical annotation releases](releases/README.md), [data and publication policy](docs/DATA_POLICY.md) and [manuscript narrative principles](docs/MANUSCRIPT_NARRATIVE_PRINCIPLES.md)

The homepage presents reviewed scientific figures and source-level statistics. Native expression, patient and sample records, coordinate tables and learned weights remain in the corresponding analysis environments. Historical runs and complete comparisons are traceable through the links above.

### Main message

**PanCancer identifies cellular combinations in multi-source atlases, interprets them through continuous states, and gives them biological meaning through tissue location and clinical context. The CTLA4–GZMB-associated T-cell core is a focal anchor for this layered ecological analysis.**

---

Scientific evidence checked: **2026-10-26**. English translation and language navigation added: **2026-10-27**. Main figures follow the order of the scientific argument; see [STATUS.md](STATUS.md) for current execution status.
