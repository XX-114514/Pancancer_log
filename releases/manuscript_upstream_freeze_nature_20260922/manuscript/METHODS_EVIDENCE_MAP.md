# Methods—evidence map

本表用于写 Methods 时逐项追溯，不替代正式文稿。

| Methods topic | Authoritative source | Version boundary | Missing before submission |
| --- | --- | --- | --- |
| Cohort inventory and discovery workflow | `methods/project-v3-unified-workflow.md`; `inventories/datasets.tsv`; run `20260725_120105_unified_pancancer_covarnet` | Earlier 31-GSE discovery object differs from 54-GSE V7 universe | Reconcile study-flow denominators; resolve GSE278694 exclusion |
| V7 identity annotation | `releases/annotation_v7_20260905/methods/V7_REFERENCE_AND_ANNOTATION_LOGIC_CN.md` | Global frozen reference | Convert to concise English Methods; add citations |
| V7 CNV/malignancy | same V7 method, config and aggregate tables | Per-sample raw counts; skips are not negatives | Verify `genome='hg20'`; prepare exact software/version table |
| V8 lineage refinements | five run records dated 2026-09-13/14; external sentinels in `EXTERNAL_ARTIFACTS.tsv` | Taxonomy only; no global writeback | Complete held-out benchmark and promotion decision |
| V9 annotation candidate | external `V9_FINAL_REPORT_20260919_CN.md` and freeze manifest v2 | `FROZEN_CANDIDATE`; S6 mostly inherits V7 CNV | 15 L3 reviews, stability decisions, 632-sample CNV decision and promotion audit |
| Integration comparison | `runs/20260725_120105_unified_pancancer_covarnet.md` | Harmony, BBKNN, scVI have different metric strengths | Select/report comparison metrics and exact subsets |
| CoVarNet K=9 | unified run record and existing evidence-chain release | Current manuscript figures are V5-derived | Rerun on final identity or explicitly retain historical branch |
| LIANA/communication | `releases/communication_evidence_chain_20260906` | V5-derived | Version-consistent rerun/impact audit |
| Spatial validation | communication evidence reports and tables | Separate validation datasets | Panel-level n, unit, test and source closure |
| Clinical/immunotherapy validation | unified and immunotherapy run records | Fixed projection; some exploratory endpoints | Verify original metadata/endpoints and final audit status |
| Data/code governance | `docs/DATA_POLICY.md`; `docs/REPRODUCIBILITY.md` | Git release is metadata-only | Select public/controlled repositories and DOI/accessions |

Methods 不能只引用当前 Markdown 摘要。最终写作必须回到对应 run、配置、统计表和
外部对象哈希；任何找不到原始命令、环境或统计分母的内容应标记为未验证，而不是补写。
