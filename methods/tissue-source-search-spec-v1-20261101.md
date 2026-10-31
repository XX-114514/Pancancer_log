# 组织单细胞来源补充：检索规格 v1（2026-11-01冻结）

## 目的与顺序

用户确定的顺序：
1. 为当前分析集中只有单一来源的14个癌种补充独立来源；
2. 按同样的必需字段复核现有43个来源，把不合格的列为候选排除；
3. 再补验证来源。

本规格在看到检索结果之前冻结。后续若要修改，只能另开v2，并写明原因，不能按结果调整。

## 目标癌种（`analysis_main`，见[癌种统一v1](cancer-type-harmonization-v1-20261031.md)）

Diffuse glioma、Papillary thyroid carcinoma、Osteosarcoma、Nasopharyngeal carcinoma、Pancreatic ductal adenocarcinoma、Pancreatic neuroendocrine tumor、Hepatocellular carcinoma、Intrahepatic cholangiocarcinoma、Bladder cancer、Esophageal cancer、Basal cell carcinoma、Acral melanoma、Seminoma、Renal medullary carcinoma。

## 纳入标准（全部满足）

| 项 | 要求 |
|---|---|
| 物种/材料 | 人肿瘤组织（可附癌旁）；整体组织**未分选**，或同时提供CD45+与CD45−两部分。排除只有免疫细胞分选、血液、细胞系、类器官、PDX |
| 技术 | 液滴或板式scRNA-seq为主分支；snRNA另列为单独分支，不混入 |
| 数据 | 可下载的**原始整数计数**（10x filtered/raw矩阵、h5、mtx，或作者计数矩阵）；只有标准化矩阵的排除 |
| 必需元数据 | 患者编号、取样部位（原发/转移）、治疗状态；新鲜/冻存尽量获取。GEO与论文/补充表都补不齐的，排除 |
| 规模 | ≥3名患者 |
| 访问 | 公开下载。受控访问（EGA/dbGaP/GSA-Human）单列，另行决定 |
| 独立性 | 与现有43个来源及已排除来源不是同一研究（PMID/DOI、GEO Super/SubSeries、CXG合集）；同团队同癌种的，人工复核 |
| 已排除来源 | GSE278694、GSE162498、GSE211602、GSE299623、GSE169609、TNBC/CTCL两个CXG，以及1,516个排除单位、83个待核单位所属的数据：不重新纳入 |

## 检索

- GEO：E-utilities `db=gds`，条件为`"Homo sapiens"[Organism] AND gse[Entry Type]`、单细胞关键词，再加各癌种的关键词。只保留"Expression profiling by high throughput sequencing"类型；记录查询串、日期和命中数。
- CELLxGENE：curation API列出全部公开合集，按`disease`（癌种）、`assay`、`suspension_type`、`tissue`过滤；记录合集DOI，用于和GEO比对是否为同一研究。
- TISCH2/3CA等整理库只用作线索，从中取原始编号，不用它们处理过的矩阵。

## 筛选流程

检索 → 候选表（自动标注：已知来源、项目记录中出现过、同一PMID）→ **用户筛选**（纳入/排除及原因）→ 下载清单（URL、大小、上游md5）→ WEHI下载并生成SHA256SUMS → 元数据侧表（GEO字段、论文补充表）；必需字段不全的退回筛选 → 现有流程：原生计数、QC、来源资格门、癌种映射 → 细胞指纹查重 → 加入分析集。
