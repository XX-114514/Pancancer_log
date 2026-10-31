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

## 进展更新（12:20 UTC）

- **漏洞7 非RNA特征排除表**（完成v1）：[non_RNA_features_v1.tsv](../../../inventories/feature_exclusions/non_RNA_features_v1.tsv)。core中有GSE274229的229个（211个抗体、14个hashtag、4个同型对照，按GEO `feature_types`判定，在库中计数全为0），以及GSE169609的5个`TOTALSEQC09xx`（类型未核，计数非零；该来源已排除）。immune和public union两个库的特征名中没有抗体/hashtag类特征（已检查）。跨来源分析只使用RNA特征。
- **漏洞3 跨研究别名**，第一层（GEO元数据）：42个GSE的PubMed ID互不重复，彼此之间没有SuperSeries/SubSeries关系；只有GSE273127（食管）和GSE274934（肺）的第一贡献者相同，是同一团队、不同癌种。没有发现同一研究重复发布。第二层（表达指纹）：新增`scripts/data_gaps/cell_fingerprint_duplicates.py`，对43个来源的每个细胞计算与基因顺序无关的计数指纹（完整非零计数多重集，以及最大50个计数），查找跨来源相同的细胞；作业30614885运行中。
- 环境作业30614882运行中（Bioconductor基础包已开始编译）。

## 进展更新（2026-11-01）

### 环境（漏洞8）
- 作业30614882跑了2.5 h后失败：conda-forge r-base不带recommended包（这一点在后续轮次已自行补上）；udunits2头文件缺失（units/sf失败）；求解器选了gcc 15，readxl等旧C++代码编译失败；Cairo缺X11头文件；systemfonts找不到`ft2build.h`。Bioconductor包装成了3.18的最新补丁版本，而不是旧版本（如zlibbioc 1.48.2，旧为1.48.0）。
- 修正：编译器固定gcc 11（旧环境为11.2）；加`udunits2`和xorg头文件；`CPATH`加入freetype2；`wehi_install_R_exact.R`的Bioconductor包先按精确版本取3.18的源码包或Archive，取不到再用BiocManager。上次失败的环境已删除（失败报告保存在`envs/build_logs/attempt_30614882_rebuild_report`），重建。
- 作业30615342在激活环境时失败（编译器激活脚本读取未定义变量，与`set -u`冲突）→ 激活前后放宽`-u`，以RESUME方式继续（作业30615350，运行中）。

### 漏洞3 表达指纹
- 首次运行（30614885）只用10分钟就算完全部指纹，但写parquet时失败（scanpy环境没有pyarrow）→ 改为写tsv.gz，重跑（30615343）。

### 漏洞1 治疗阶段
- 单位↔GSM桥：按GSM编号、标题完全相同、标题方括号内样本名、唯一子串的顺序匹配，731个GEO来源单位中549个唯一匹配（脚本`scripts/data_gaps/unit_geo_characteristics_v1.py`；单位级结果只留服务器）。
- 证据：GEO系列声明3个来源、GSM字段6个来源、样本名标签（GSE182109的ndGBM/rGBM/LGG）、论文全文（PMC，经E-utilities下载到WEHI，只抽取治疗相关句子）。证据表[treatment_evidence_v1.tsv](../../../inventories/data_gaps/treatment_evidence_v1.tsv)。
- 结果（[按来源汇总](../../../inventories/data_gaps/gap1_treatment_resolution_by_source_v1.tsv)）：466个待核实单位中206个有了证据。其中未治疗141个；有既往治疗或其他情况43个（GSE294109的19个是ICB/TKI治疗后切除，GSE182109的18个复发GBM，GSE184362的5个是复发前已手术或碘治疗）；GSE182109的4个LGG论文未说明。剩余260个单位分布在17个来源。其中GSE148071、GSE274229、GSE152048论文写明来源内有混合治疗状态，需要补充表；GSE132465、GSE222315没有PMC全文。
- 同时得到细亚型的单位级信息：GSE182109分级（IV 35 / II 8 / III 1）、GSE161529受体亚型、GSE159115、GSE197461、GSE273127、GSE274934的组织学（见`unit_geo_characteristics_by_source_v1.tsv`），用于癌种v2的单位级细分。

### 漏洞3 结论（表达指纹，作业30615343，9.6 min）
- 43个来源4,906,058个细胞（总计数≥500、检出基因≥200）。完整指纹跨来源相同0个；top50指纹跨来源18组、9个来源对（每对1–9个），视为低复杂度巧合。
- 来源内"重复"305,390个，全部来自GSE294109和GSE256326，每个恰好2份——这两个来源在immune和public两个库各有一份，本筛查把两库都读了。这同时证实两库副本的计数向量逐细胞完全一致（此前只核对过细胞ID）。
- 结论：43个来源之间，从GEO元数据和细胞级计数都看不到同一份数据被重复收录。无法排除的是：同一患者在两个研究中各有不同样本，计数矩阵检测不到这种情况，需要原始reads的基因型。

### 漏洞2 患者键（[按来源汇总](../../../inventories/data_gaps/gap2_patient_keys_by_source_v1.tsv)）
- 381个无患者依据的单位中，123个按来源命名规则（对照GEO标题与论文）归并为75名患者：GSE274229 44→25（阶段加文件名第二个编号；6 mCRPC + 14 mHSPC + 5 localized，与论文一致；该编号只在同一阶段内唯一）、GSE127465 18→7、GSE294109 25→15（多区域）、GSE178481 15→9（多区域）、GSE155698 16→15、GSE138709 5→4。即这些来源原本多算了48个"独立单位"。
- 245个单位按标题判断为一个单位一名患者；GSE163558（原发/转移样本无法配对）和GSE228501共13个单位仍未知。
- 有记录的389个单位保留原patient_candidate_block。
