# 泛癌注释marker、panel与作者代码逻辑核查

取证节点：2026-11-04 23:14:30 UTC / 2026-11-05 07:14:30 UTC+8；本机日期2026-11-05 Australia/Melbourne。用户要求逐篇整理注释marker/panel，从作者代码解释注释过程，并提出本项目如何学习和借鉴。

## 起点、范围与计划

前一轮已交付8项研究的680条分层标签。本轮从标签回溯代码、正文明确panel、作者功能签名、DEG及蛋白验证，明确区分这些证据的用途。读取仓库AGENTS、README、STATUS、ROADMAP、TODO、RUNTIME_SYNC及前轮任务日志；主工作区已有大量并发修改，全部保留。未连接服务器、未启动注释/训练、不改V7或当前scVI工作图谱。

计划：六项可取得作者代码的研究做静态检查；Qian以公开作者PDF及补充表核查Methods，Yang如不能取得可核验代码则保留缺口；提取明确基因/蛋白列表并保留代码行/表位置，另给本项目建议，不将建议冒充作者实现。

## 实际读取的原始资料

- [Qian作者PDF](https://www.qian-lab.org/assets/pdf/single-cell_pan-cancer_blueprint_2020.pdf)、[S5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0355-0/MediaObjects/41422_2020_355_MOESM17_ESM.xlsx)、[S10抗体库](https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0355-0/MediaObjects/41422_2020_355_MOESM22_ESM.xlsx)。S10确为198行，保留蛋白原名、Clone等，不自动当成RNA基因。
- [Lodi作者仓库](https://github.com/lambrechtslab/Pan-cancer_Lodietal)：Cancer type analysis.R、T/B/DC/EC/MM五谱系脚本、Akoya.R、marker_genes.txt及前轮分层映射；检查MarkerToAdd/xx占位、具体分辨率与手工覆盖顺序。当前仓库树取证SHA为35f054bb81e1b321718495e096171580b3bc130c，旧下载文件另以本地快照为准，不将源文件发布历史等同于论文冻结版本。
- [Chu TCM](https://github.com/Coolgenome/TCM)：FindCluster.R、snn-marker.R、qc-by-cluster.R、bubble-plot.R、fig1/figure1_v2.r、fig2/figure2.r、TCellMap.R、TCellMap2.R及GSE169246查询实例split.R；[补充表S3–S10](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-023-02371-y/MediaObjects/41591_2023_2371_MOESM3_ESM.xlsx)。仓库树869e7c9db553ba04f1231a3bf2655e533157f06c；下载文件Git blob SHA在工具回执保留。
- [Yang GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233236)：沿用20个公开标签作为范围索引；Cell正文403，本轮未取得可核验完整作者注释代码/marker表，不将标签后缀拼成完整panel。
- [Kang Zenodo](https://doi.org/10.5281/zenodo.10651059)：已有公开analysis_codes.tar中Step1.NMF_processing.py、Step2.NMF_module_visualization.ipynb、Step4.Reference_component_analysis.ipynb，静态读取代码单元；[作者Supplementary Data2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-48310-4/MediaObjects/41467_2024_48310_MOESM5_ESM.xlsx)五组99个状态列表。
- [Luo PanCAF](https://github.com/Xiaxy-XuLab/PanCAF)：Main_Part1.R、Main_Part2.R与Supplementary Data2/4，区分数字群、CAF解释、局部比较panel及缺实际标签的C24步骤。
- [Ma DC原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/)、[作者DC_altas代码](https://github.com/matyhoc/DC_altas)：Other_analysis.R、ROGUE.R、README、DCTypist/dc_predict.py，仓库树c486a586984e9c44290e257bc35edc24f775a224；核心三轮聚类从原文Methods核查，不冒称公开仓库含完整初始注释流程。
- [scATOMIC](https://github.com/abelson-lab/scATOMIC)：run_scATOMIC.R、classify_layer.R、add_unclassified_automatic_threshold.R、automatic_threshold.R、参考准备代码；检查父层路由、拒判与可停留父类，未展开全部模型分支特征或运行预测。

## 已观察的关键结果

1. Lodi先按癌种认大类，合并同谱系，初轮检查低质量与doublet，清理后重算谱系HVG/PCA/整合和聚类；DEG在RNA assay计算。T基线res0.7并用更高分辨率补稀有群；B因Harmony不足改CCA，用res0.8/1/1.9/2分别补标。这是人工多轮注释，不是一套统一分辨率或marker argmax分类器。
2. Lodi公开亚型绘图panel存在`MarkersToAdd`，大类RenameIdents有`xx`占位；marker_genes.txt实测仅21个唯一cluster。不能据其文件名声称70共享亚型的完整marker都已取得。代码还有objectNew、regression_conditions、Ident拼写及res3.6不在所列扫描网格等不完整之处，不能原样复制运行。未分配细胞默认Low quality的做法不作为本项目规则。
3. Qian原文按已知类型分离和每群至少10个两检验FDR<0.01、2倍DEG选择分辨率；从B和髓系中找回DC，再合并子聚类。S5 73个分层标签中72个可提取DEG；功能表C2_CMA1没有对应Mast cell DEG列表，单列缺口。PRDM1high/PRDM1-high等原名差异通过明确别名桥接，仅为找回来源，不覆盖作者原标签。
4. Chu补充S4/S6共有36个完整功能签名；明确绘图panel包含Fig1D 15基因、Fig2B CD8特征、p-Tex/Tex比较。身份、功能与状态证据各有用途。TCellMap使用FindTransferAnchors和MapQuery；代码有未定义hvgQ和外部路径依赖。一个查询实例的CD4与max(CD8A,CD8B)差值阈值0.1527只属于该实例，不提升为全图谱阈值。
5. Kang Step1输入已经有anno_category，不是初始身份注释代码；谱系×Patient_Organ_Tissue做K5–9 NMF。Step2检查MT/RP、单研究和跨谱系ambient；源代码10/150权重、150模块、top50中2个MT/RP、0.2分数差等参数为具体示例，不能直接用于我们所有谱系；原NMF未固定random_state。
6. Luo的cluster与解释标签数量不同，DC重聚类后将11个初始数字群合为6群。代码父编号语境及C24结果列表仍有限制；绘图基因可能含竞争谱系和ambient候选，不能将所有绘图基因合成一个阳性判定gene set。
7. Ma原文为髓系→DC→DC亚群三轮；明确正向与竞争身份panel。ROGUE是群内异质性指标，不是身份正确率。LAMP3来源用SingleR及参考top10 DEG、主分析max score>0.3，仍为计算推断。DCTypist代码调用CellTypist，公开类名与最终论文有差别；代码为标签多数投票且np.random.choice默认有放回，正文称模型平均。这里记录差异，不在未核清版本时擅自对齐。
8. scATOMIC低分子类可以回到Dendritic Cell或Macrophage or Monocyte等父类；参考准备明确排除ILC、Eryth、Platelet及部分非常规T。因此不能用其预测否定这些稀有类别，也不能将自动阈值当本项目通用阈值。

## 交付与建议的性质

生成`${OUTPUT_ROOT}/pancancer_annotation_code_panels_20261105/泛癌注释marker-panel与代码逻辑.xlsx`，三个工作表：逐篇代码与借鉴、作者marker与panel、本项目建议。作者资料表595行，包括：22代码panel、11正文panel、1正文蛋白判别组合、2蛋白hub panel、198个CITE-seq抗体条目、36个完整T功能签名、3个正文完整功能签名、99个完整NMF状态列表、201条DEG原表前10摘录，以及22条明确缺口（Qian1、Yang20、scATOMIC1）。条目数不是独立基因数、细胞类型数或全部可用panel数。

DEG摘录明确是原顺序前10，不等于完整DEG表或人工判定panel。完整签名保留原顺序和重复项；保留蛋白名称、OSA1等原写法；TCellMap跨行HSP90AB1仅拼接源token并备注，不静默大规模基因别名标准化。作者panel和15行本项目建议分表；建议panel没有经过本项目表达验证。

综合建议：父身份、cluster_uid、marker名称、可选标准亚型、状态和不确定原因分字段。父身份可靠但子类缺证据可标coarse_supported；父身份冲突、doublet/低质与参考缺失须分别检查。红系/HSPC/ILC可在证据允许的粒度停止，但不是这些类在生物学上没有亚型；尤其HSPC需排除内皮CD34共享，ILC不能仅以TCR dropout认定。身份构图如评估去除IFN/应激等影响，要保留未抹除状态的RNA层用于状态读出；本项目计划状态输入，因此不照搬作者全部状态回归。

## 执行、验证与失败

本地自编inspect_sources.py读取公开资料并静态抽取Kang notebook代码；build_panel_data.py提取有来源位置的列表；artifact-tool单一build_workbook.mjs生成XLSX；verify_workbook.py以ZIP/XML读回8×7、595×10、15×6字段，与源JSON逐项一致。三个原生筛选表、冻结首行/首列验证通过；三表预览目视检查。工作簿SHA256：`dd33d09a33520f1281ddcbb88eea2b05fbe7193235f0df285d7bda3ed945d476`。

公开Qian PDF及抗体表只读下载；抗体表SHA256为c1d7243ac859b43912cfd16539dd4c05c4b4384874833dd42a42cb696a6647e3。未下载新的逐细胞人类数据、未执行作者分析代码或加载DC pickle模型。原下载和脚本保留本地，公共仓库只提交本次脱敏摘要及daily引用。

失败及修正：环境没有bs4，改用Python标准库HTMLParser；Nature下载页部分只含摘要/引用，未冒称读到缺失Methods，Qian改读作者公开PDF、Ma读正常开放PMC；Chu PMC验证码、Yang Cell403与Zenodo页面超时均未绕过，使用已有正常公开补充表/作者代码/归档。起初API_QUICK_START路径多写references，定位到artifact_tool_docs后正确读取。作者脚本中发现的不完整变量与命名差异保留为研究边界，未修补后冒称复现成功。

仓库轻量validator实际exit0，0 errors / 1项既有warning，无失效本地相对Markdown链接；git diff --check通过。隔离检出基于远端7bf93091b1c6be1463f3d481f0c0dac9ab667875，仅准备本次两份脱敏日志，普通同步后独立ls-remote核验并另存本地收据。建议commit：`Review author marker panels and annotation code logic`。后续建议先以T/NK和B/ASC各选一组本项目cluster建立同一套证据审阅表，核对父身份、群DEG、竞争身份、状态、来源支持与细胞级参考预测；本次未启动该实验。

## 追问：QC顺序、父标签与计算方法澄清

2026-11-05 00:26:58 UTC / 08:26:58 UTC+8补充取证。用户追问Lodi先分大类后QC与本项目前置QC的区别，父标签为何需复查，以及marker计算差异。明确纠正前轮易误读的简写：Lodi的Cancer type analysis.R第162–194行在大类之前已经做表达量过滤、按orig.ident逐样本DoubletFinder及第二档表达过滤；谱系脚本后面的QC为二轮复查，主要检查已有pANN预测、质量分布及跨谱系表达，并非证明必须在合并谱系后重新运行doublet算法。前置QC与后置复查互补，本项目先QC本身没有错误，也无证据说明Lodi仅因该顺序就更准确。

父标签不是不可修正的边界：Qian Methods明确pDC可与B共聚类、其他DC在髓系，需要用marker找回DC后再合并细分。对一个群同时有T/DC信号，必须区分不同细胞子群混在同群与同一barcode的竞争程序；不能仅凭群均值同时有两类marker就宣布doublet，也不能仅凭一个阴性marker或RNA dropout排除身份。

读取本项目V7冻结方法、20260823注释总结、20260903 marker图口径与annotation_generation_followup，比较范围限定首页旧5M/V7，不把十月不同分支混为同一算法。V7有来源逐细胞标签、真实全库marker表达和CellTypist仲裁；已记录top score≥0.25、margin≥0.08、负marker与跨来源/样本门槛。旧谱系细分重建BBKNN/Leiden但复用共享PCA，不能说本项目没做谱系内聚类。V7面板展示359基因的群内比较不是全转录组DEG；V7后续有全基因描述性补名，二者不能混同。

计算对照：Lodi的RNA LogNormalize和本项目log1p(CP10K)属于同一基础变换，只有相同基因宇宙和文库分母时才逐值可比。Lodi FindAllMarkers在RNA assay找群特征，人工结合canonical/竞争身份后RenameIdents；公开亚型panel有占位，未给完整自动身份打分公式。AddModuleScore为基因集均值减表达匹配对照均值，代码中主要用于IFN/应激/缺氧等功能，不等于最终身份分类器。不能未核旧评分源代码就声称V7分数与AddModuleScore数值相同。

另外核查Seurat官方v4.3.0 DotPlot源码：avg.exp为mean(expm1(data))，显示颜色再log1p和按基因缩放；本项目旧v5图文档为平均logCP10K后逐基因min-max。mean(log1p(x))与log1p(mean(x))不相同，图色不能直接当注释正确性证据，也不冒称Lodi实际版本必为v4.3.0。官方NormalizeData/AddModuleScore/FindAllMarkers/DotPlot文档已查。

本项目后续refine_immune_native_markers.py第74–99行有mean CP10K、delta detection、log2((mean+0.01)/(rest+0.01))及跨来源支持，是明确标作无p值的描述性命名；measure_TME_full_cluster_markers_v2.py第73行已经实现Scanpy官方Wilcoxon/BH/common-gene比较，故不能将项目整体概括为只有固定panel或从未实现DEG。这里只核代码存在，不刷新后续运行是否完成/被V7使用的状态。没有修改科学代码、工作簿、细胞标签或服务器作业。

检索失败及边界：首次猜测marker审计文件目录失败，改由文件清单定位到07_fine_annotation；一次Seurat源码筛选正则转义错误在执行前失败，改正后只读成功。原V7评分全部聚合公式未在本轮本地冻结脚本中重新展开，因此只比较可核验的归一化、规则门槛和证据组织，不捏造精确权重。建议commit：`Clarify initial and lineage QC and marker annotation calculations`。

补核[Seurat FindMarkers官方定义](https://satijalab.org/seurat/reference/findmarkers)与[Scanpy rank_genes_groups官方定义](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html)：Seurat常规p_val_adj为按全基因数Bonferroni；本项目上述后续Scanpy代码显式使用BH，并非可直接逐值比较的同一输出。Scanpy文档也注明logfoldchanges为基于mean-log的近似；不得将其与所有Seurat版本或本项目CP10K+0.01描述性log-ratio强行等同。比较cluster背景（全TME其余 vs 同谱系其余）本身也影响marker排名；同谱系不显著的共享身份marker仍可为父身份提供表达证据。

本次追问轻量validator再次实际exit0，0 errors / 1项既有warning；仅为新增澄清日志检查，不是重跑或验证科学注释。仅追加本任务Codex/daily记录，原工作簿和源代码不改。
