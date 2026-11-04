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

### 漏洞5 材料信息
- 取样部位（[按来源汇总](../../../inventories/data_gaps/gap5_sampling_site_by_source_v1.tsv)，脚本`gap5_sampling_site_v1.py`）：521个"肿瘤阶段未注明"的单位中，449个按GEO source name / tissue字段或来源设计判定为原发器官部位（GSE182109、GSE138665、GSE138709、GSE222315按来源规则），3个为转移部位（腹膜），69个未解决（GSE179994 46、GSE152048 11、GSE139829 8、CXG SPECTRUM 4）。取样部位与疾病分期分开记录：GSE274229前列腺活检来自mHSPC/mCRPC患者，GSE305078含advanced TNBC。
- 新鲜/冻存（[证据](../../../inventories/data_gaps/fresh_frozen_evidence_v1.tsv)）：GEO简要视图不含protocol文本，改为从论文方法部分抽取。429个待核实单位中，188个（9个来源）论文写明为新鲜处理；GSE299340原文有歧义（"fresh…stored at cryogenic temperature"）；GSE161529只写明正常组织为新鲜；其余195个单位（11个来源）方法中没有写明。

## 2026-11-01 用户新指令：先补14个单来源癌种，再复核现有来源，最后补验证来源

### 环境（漏洞8）基本完成
- 作业30615350：510个包为精确版本，余17个不一致。根因：ragg 1.2.6与当前freetype头文件不兼容；lme4 1.1-34不能用Matrix 1.6-5的CHOLMOD头文件编译；旧环境中部分Bioconductor包是3.17版本。
- 修正作业30615876（`fix_R_full_env.sbatch`；第一次提交时`install_version`没有设置CRAN源而失败）：ragg用`-fpermissive`；lme4改用1.1-35.5，登记为版本偏差（`rebuild_report/R_version_deviations.tsv`）；安装脚本同时查3.17和3.18的源码包。结果：Seurat、CellChat、copykat、scATOMIC、infercnv、NMF、SingleR、GSVA全部LOAD_OK；GitHub包按commit安装（COSG、SingleR、sceasy为复制旧编译件）。presto已改回1.0.0。monocle3需要旧环境没有的`terra`，正在补装（作业30615917）。

### 检索规格与候选（第1步）
- [检索规格v1](../../../methods/tissue-source-search-spec-v1-20261101.md)在看结果之前冻结。
- `search_tissue_sources_v1.py`：GEO E-utilities检索14个癌种，加CELLxGENE curation API。共969行（GEO 953、CXG 16），检索日志见`inventories/source_search/search_log_v1.json`。
- `screen_tissue_candidates_v1.py`：第一层按系列级规则过滤（人、3–300个样本、非已知来源、非功能实验或bulk、有发表、有计数文件、标题为肿瘤组织）；第二层对每个癌种前20个逐个读GSM记录，看物种、患者/组织/治疗字段、分选/血液/非scRNA用词。共152个进入第二层，结果见`shortlist_auto_v1.tsv`。
- 人工预判（依据标题与GSM字段，**尚未核原始矩阵与论文**）：[candidates_curated_v1.tsv](../../../inventories/source_search/candidates_curated_v1.tsv)，A/B/C三级，共63条。胶质瘤没有A级（自动检索偏向plate-based研究，需要定向补检）；PanNET、RMC只有B级。

### 现有来源复核（第2步）
[existing_source_fitness_v1.tsv](../../../inventories/data_gaps/existing_source_fitness_v1.tsv)：按冻结的必需字段（患者编号、取样部位、治疗状态、≥3名患者）复核。保留23个来源（449个单位/357名患者）；待补充表17个（302/255）——补不齐的排除；建议排除3个：GSE163558（患者无法配对）、CXG SPECTRUM卵巢（6个单位来自同1名患者）、GSE228501（患者与治疗都不明）。**尚未执行，等用户确认。**
- monocle3 1.3.1（GitHub tag）与presto 1.0.0均LOAD_OK（作业30615914/30615917，`jobs/fix2_monocle_presto.sbatch`）；为满足monocle3，新增旧环境没有的terra 1.9.50。偏差全部记录在`scrna_r_full/rebuild_report/R_version_deviations.tsv`：lme4、terra、monocle3来源。**漏洞8的R环境部分完成。**

### 用户确认排除3个来源（2026-11-01）
- 用户确认排除GSE163558、CXG SPECTRUM卵巢（`CXG_b252b015…`）、GSE228501，**不再整合**。登记于[source_exclusions_v1.tsv](../../../inventories/data_gaps/source_exclusions_v1.tsv)，`inventories/datasets.tsv`新增3行。
- 分析集v2：751个单位、40个来源、2,893,337行（服务器`analysis/data_gaps/analysis_units_v2_private.tsv`，SHA256 `9cb95a79…c8eb`）。排除只体现在单位清单；冻结存储与原生数据不改动，也不物理删除（被排除的细胞与其他细胞在同一个zarr中）。
- 影响：精原细胞瘤在分析集中没有来源了（只能靠新来源补）；卵巢癌、胃癌各剩1个、2个来源。
- GSE163558的排除理由只依据GEO记录，没有查论文。

### A级候选详细核查（进行中）
`detailed_check_candidates_v1.py`（Slurm `detailed_check_candidates_v1.sbatch`，作业30618920）：对22个A级编号逐一读取GEO系列与GSM完整记录（物种、技术、分选用词、样本来源、患者字段），列出补充文件（原始计数与标准化文件），下载1个<60 MB的计数文件检查是否为整数，并从PMC论文方法部分抽取治疗、新鲜/冻存、分选、患者数相关句子。第一次用nohup在会话shell中启动，进程无输出即消失，改为Slurm。

### A级候选详细核查结果（2026-11-01）
- 自动部分：作业30618920（4 min）。整数探针第一版有误（会选中`filelist.txt`），改为`count_file_probe_v1.py`：从GEO记录字段取文件名，只流式读取压缩文件的前几MB并解压检验。Slurm作业30619015里有9个系列没有解析到数值（可能是下载时的瞬时问题），直接重跑后全部有结果。
- 人工汇总：[candidates_detailed_check_v1.tsv](../../../inventories/source_search/candidates_detailed_check_v1.tsv)。22个编号中，通过13个，有条件7个，不通过2个（GSE243238是bulk RNA-seq；GSE197778只有1名患者）。
- 主要条件与风险：
  - GSE150430：GEO上只有标准化矩阵，原始数据在CNGB，需要从FASTQ重新比对。
  - GSE151530：很可能包含GSE125449的样本（同一NIH团队，论文引用了后者），GSE189903也来自同一团队，需要用细胞指纹去重。
  - GSE206245：GEO只存放了论文队列的一部分；它与GSE150430、现有的GSE162025来自同一机构，需要查患者重叠。
  - GSE301651：属于新辅助免疫化疗队列，scRNA患者的治疗状态未明。
  - GSE254918：每个癌种只有3名患者，刚好达到下限。
  - GSE160269：CD45+与CD45−两部分需要按患者合并。
- 癌种覆盖（只算通过的）：PTC 2、骨肉瘤 2、PDAC 2、肝细胞癌 3（GSE149614、GSE125449、GSE189903）、ICC 2（GSE125449、GSE189903中的ICC部分）、膀胱癌 1（GSE277524）、食管癌 2、肢端黑色素瘤 1（GSE215120）。NPC、BCC只有有条件的候选；精原细胞瘤没有可用候选；胶质瘤、PanNET、RMC还没做A级核查（没有A级候选）。

## 2026-11-01 下载、处理并构建分析就绪zarr（用户要求）

- 13个"通过"的补充来源已下载、读取并完成QC，与现有分析集v2合并为`pancancer_analysis_ready_v1.zarr`：3,818,042个细胞、53个来源、1,020个单位、793名患者、22个癌种、12,498个共同基因。逐步记录见[run](../../../runs/20261101_analysis_ready_zarr_v1.md)，数据说明见[docs/ANALYSIS_READY_ZARR_V1.md](../../../docs/ANALYSIS_READY_ZARR_V1.md)。
- 独立回读验证全部通过；跨来源完整指纹相同的细胞为0。
- 构建后发现GSE125449的P编号在两个Set之间重复且对应不同癌种，不能作为患者ID，改为按样本计患者（只改obs）。
- 未完成/限制：现有部分没有扣除环境RNA，补充部分没有估计环境RNA；GSE163203低深度；GSE277524只用QC下限识别细胞；GSE160269为分选后合并；精原细胞瘤没有来源；数据在scratch上，持久空间尚未落实。

## 2026-11-01/02 用户指令：新来源去富集/低深度，统一QC+双细胞+环境RNA，整合降维出UMAP；扩充数据（覆盖TCGA大多数癌种，约500万细胞）

### 新来源的排除规则（用户确定，适用于所有新来源）
- 不使用富集数据：谱系分选（FACS/MACS/CD45±等）一律排除；只去死细胞的活细胞FACS不算富集，保留。
- 深度与其他数据差异过大的排除：来源中位UMI低于各来源中位数之中位数（4,127）的1/3、或高于其3倍即排除（现有40个来源的范围为0.43–2.95倍）；测到的基因少于15,000的排除（作者已做过基因过滤）；少于3名患者的排除（检索规格既有条款）。
- 据此排除：GSE160269（CD45±分选）、GSE163203（Microwell-seq，0.16倍）、GSE244983（11,622个基因）、GSE225857（免疫/非免疫分开计数，无法排除分选）、GSE117570（10,147个基因）、GSE190888（2名患者）。

### 环境RNA统一处理
- 现有来源原本不一致：core只有估计值没有校正矩阵；部分immune/public来源用了其他方法。现对**全部来源**统一做decontX：按单位快速聚类后调用项目既有worker（`scripts/decontx_library_worker.R`算法与参数不变；独立的decontX包需要rstan，编译失败，改用celda 1.18.2中的`decontX()`，复制为`decontx_library_worker_celda.R`）。
- 原51个来源全部完成，没有单位失败；各来源中位污染率0.2%–32%（中位3.8%；最高为GSE138665 32%、GSE141526 17%、GSE179994 16%）。整合时去掉污染率≥0.5的细胞。

### 数据扩充（检索v2/v3）
- v2（TCGA缺失或不足的癌种）：31个候选，判定见`inventories/source_search/v2/candidates_detailed_check_v2.tsv`；读取14个（HNSCC 5、皮肤黑色素瘤 3、卵巢 2、子宫内膜/卵巢 1、软组织肉瘤 1、DLBCL 1、间皮瘤 1），全部通过深度和基因数检查，decontX运行中。
- v3（常见癌种的独立来源）：27个候选，见`v3/candidates_detailed_check_v3.tsv`；读取11个，其中2个按上述规则排除，余9个在做QC。与现有来源属同一研究的排除（GSE144735与GSE132465同一论文；GSE208653与GSE197461同一论文）。
- 读取中遇到的问题：GSE212527矩阵未过滤（改按作者细胞列表）；GSE181919材料字段拼接错误；GSE173278文件名不同；下载脚本第二次运行时重写只读的SHA256SUMS失败（已改为只给新目录写）。

### 整合流程（integration v2，脚本`scripts/integration_v2/integrate_v2.py`）
沿用此前TME标准整合配方：seurat_v3 HVG（按来源分别计算后合并排名，6500个候选→4000个，去掉技术基因），decontX计数normalize_total 1e4 + log1p，中心化ARPACK PCA 50个主成分，按来源做Harmony（theta 2，lambda 1，20轮），15近邻，UMAP。等全部来源完成decontX后在bigmem节点运行。

## 2026-11-03 恶性细胞判定 v1 全量运行检查（UTC 11:35 / UTC+8 19:35）

- 已观察：array 30645933 完成 68/72 个来源（1,162个单位）；无 Traceback；inferCNV 未运行的66个单位都没有候选细胞（全为免疫/基质），不属失败。160个单位借用了同来源的参考细胞。调用合计：non_malignant 3,003,312；malignant 1,025,620；uncertain 361,940（全部为 methods_disagree）。
- 已观察：4个作业（GSE125449、GSE183904、GSE253013、GSE269936）在 scATOMIC 中停止进展超过20 h（最后写入在2026-11-02 08:54–12:51 UTC）。无法登录计算节点查看进程。已取消，R worker 增加超时（scATOMIC 2 h，单核重试 4 h；inferCNV 6 h）并按进程组终止，从 pass1 缓存续跑（job 30653005）。
- 已观察：GSE299023（骨肉瘤）0% 恶性；GSE152048 3.3%；GSE182434（DLBCL）0.1% 恶性、22.9% uncertain，scATOMIC 把 CNV 阳性 B 细胞判为 "B Cell"。解释：肉瘤/间皮瘤的恶性细胞表达 COL1A1/DCN，被当成基质参考；scATOMIC 没有淋巴瘤类别。修正：DLBCL、骨肉瘤、软组织肉瘤、间皮瘤采用 CNV 主导模式。成纤维/周细胞样聚类改为候选，参考只用免疫和内皮细胞；CNV 阳性即判恶性。6个来源重跑（job 30653009），旧结果移到 `malignancy_v1/superseded_rule_v0/`。
- 已观察：GSE179994 全部46个单位100%为免疫谱系，没有基质或上皮细胞；采样部位原本未解析。解释（假设：T 细胞分选设计）：属谱系富集数据，按用户"不用富集数据"的规则不进入 TME v1。保留在 integration v2，未删除。
- 待办：剩余10个作业完成后运行 integrate_tme_v1、report_dataset_detail_v1，并补全 docs/MALIGNANCY_TME_V1.md。

## 2026-11-04 四个高免疫比例来源核查与决定（UTC 04:40 / UTC+8 12:40）

- 用户要求（2026-11-04）：核实第一项；若确认未富集就先不排除；不确定细胞本次不进入 CoVarNet；所有决定标明清楚。用户提供了四个来源的原文方法证据。
- 已观察（本分析对象）：
  - GSE201347 只含12个肿瘤样本，GEO 的10个非肿瘤样本未纳入。
  - GSE162025 只含10个肿瘤样本，无 PBMC，共80,516个细胞，与原文肿瘤82,622个接近。
  - GSE189903 含 T1–T3（核心）和 B（边缘），N（癌旁）未纳入。
  - 这四个来源的 author_cell_type 都为空，无法做作者标签映射核对。
- 已观察（恶性比例偏低的原因）：共识规则要求 scATOMIC 和 CNV 同时为阳性。各来源"scATOMIC 正常 + CNV 阳性"的候选细胞数：GSE189357 3,808；GSE201347 5,030；GSE189903 1,296；GSE162025 1,197。
  - 在共识规则下，这些细胞一部分归入 uncertain，一部分归入 non_malignant。
  - GSE162025 的 malignant 加 uncertain 共3,176个，原文恶性上皮细胞为2,787个。
  - 解释：比例偏低是共识规则保守所致，不是标签映射错误。所有百分比都是本分析对象中的值，不是原文报告值。
- 决定：
  - 四个来源都保留，标为"未报告免疫特异富集"。
  - 处理信息分别记录：胃癌去死细胞；肝癌细胞悬液冻存复苏；NPC 活细胞（PI 阴性）分选，且有单独分离的配对 PBMC，PBMC 未纳入本对象。
- 纠正：
  - GSE189903 的 fresh_frozen 由 fresh 改为 cryopreserved_or_frozen。证据：PMC9729309 Methods 写明悬液在液氮中冻存，复苏后用于 10x 5′。旧记录有误。
  - GSE179994 的排除理由改正。GEO 方案没有分选；公开文件只有作者的 T 细胞子集（all.Tcell.rawCounts）。此前写的"T 细胞分选设计"是未经验证的假设，是错的。仍不进入 TME v1 和 CoVarNet v1，但理由改为作者事后计算筛选。
  - uncertain 细胞保留在 TME v1 zarr 中（covarnet_v1_include=False），不进入 CoVarNet v1。CoVarNet v1 共3,015,706个细胞。
- 全部决定记录在 `inventories/tme_v1/source_decisions_v1.tsv`。TME zarr 在整合完成后用 patch_tme_decisions_v1.py 写入。
