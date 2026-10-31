# 2026-10-31 数据漏洞1–8作为持续目标：癌种统一与scrna_r环境重建

## 用户要求

把剩余数据的8项漏洞（治疗阶段、单位≠患者/伪重复、跨研究别名、癌种名称、材料信息、独立验证、非RNA特征、计算环境）作为持续目标推进，**优先**癌种名称统一和环境构建。用户同时定下规则：排除的数据不再整合；物理删除需另行确认清单。

## 开始状态

- 仓库`f7fd872`之后；分析范围为770单位/43来源/2,968,368行（`20261027_TME_module_scope_registry_v2`），单位材料表`20261030_TME_CM_technical_material_design_v1`。
- WEHI上R环境只有CoVarNet子集`scrna_r_covarnet2`；完整scrna_r未重建（Anaconda默认源不可达）。

## 癌种统一（完成v1）

方法和结果见[方法](../../../methods/cancer-type-harmonization-v1-20261031.md)，表[cancer_type_map_v1.tsv](../../../inventories/cancer_type/cancer_type_map_v1.tsv)。取证在2026-10-31 11:35—11:55 UTC：
- GEO系列文本42个、GSM特征17个系列；OncoTree官方API（stable 2025-10-03）。
- 与770单位表关联：770/770匹配，0个原标签冲突；35个写法归为23类；有≥2个来源的分层9个（595个单位）。
- 纠正GSE182109标签（混合级别胶质瘤）；GSE297041的`CESC`没有依据。

## 环境重建（进行中）

- 新增`scripts/wehi_migration/jobs/build_R_full_env.sbatch`与`wehi_install_R_github.R`：conda-forge r-base 4.3.3（与旧环境同一build `h5074ccb_16`）加系统库，一次求解；然后按旧版本安装全部523个CRAN/Bioconductor包（复用`wehi_install_R_exact.R`），9个GitHub/其他包按记录的commit安装，失败时复制旧编译件并标注。
- 失败1（作业30614865）：`imagemagick`需要ICU ≥78，与r-base build的ICU 75.1冲突，求解失败。逐个干跑定位：`libgdal-core`与`libxml2-devel`（仅libxml2 ≥2.14有）冲突。修正：去掉imagemagick和libxml2-devel，固定`libxml2=2.13`（旧环境2.13.9）、`hdf5=1.14`（求解器原本选了2.2.0；旧1.14.5，hdf5r 1.3.12需要）、`jags=4.3`（原本选了5.0.1；rjags 4-17需要JAGS 4）。干跑通过：icu 75.1、hdf5 1.14.6、jags 4.3.2、libgdal-core 3.10.3。
- 已提交作业30614882（环境`${PC_ROOT}/envs/conda/scrna_r_full`），待完成后检查`rebuild_report/`和核心包加载结果。

## 其余目标（未开始）

3 跨研究别名、1 治疗阶段、2 患者/伪重复、5 材料、6 独立验证、7 非RNA特征排除表——按优先顺序在环境构建完成后推进。

## 建议commit message

Harmonise cancer-type labels to OncoTree for the 43 qualified sources; rebuild the full R env on WEHI
