# V7文献指导的逐类注释复查

## 请求、范围与初始状态

用户要求严格参考文献，从冻结V7逐类排查unresolved/coarse并推进完整注释。文献骨架为Lodi 2025谱系内复查和重新建模、Qian 2020 marker/竞争谱系核查，T细胞辅以Chu 2023；公开脚本中的占位符和未定义参数必须显式标记，不能宣称逐参数复现。

本节点取证UTC 2026-11-05T02:29:53.460713+00:00；UTC+8 2026-11-05T10:29:53.460713+08:00。原始V7不改，修订使用新目录。主工作区已有14个tracked改动和大量untracked；全部保留。公开发布只使用独立干净检出和显式路径。

## 已执行与观察

- 通过既有Tailscale跳板只读访问原项目档案，读取服务器祖先AGENTS；不删除、不覆盖旧对象。最初以近期项目根查找未命中，依据历史采集器纠正为V7原项目根；没有迁移证据。
- 核对冻结checkpoint SHA256：`e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e`，4676787细胞×128字段，index唯一。
- 冻结字段major/final/confidence与原始cluster、requires_recluster、strict eligibility等同时存在，可以区分父类改变后未重新分群与旧cluster未达亚型门槛。全量诊断及历史汇总表schema已只读取得；原始细胞/样本行未公开导出。
- 现成环境：Seurat5.3.1、Harmony1.2.4；另有同规模95038基因原生矩阵，尚须逐轴和来源核验，不能用其后加标签替代V7。
- `Tcell_4`旧扩展panel中GZMK/GZMA/CCL5通过，690/698样本支持，但CD8_GZMK与NK竞争得分margin0.00723，未通过旧规则；这只解释旧程序未决，不证明T/NK身份。`Tcell_6`旧扩展panel仅TIGIT通过、经典锚点0，不可凭BATF命名认定Tex。
- Platelet旧coarse记录的原因是`rare_lineage_no_forced_subtype`，本身并非注释失败。

## 执行记录与失败

命令采用 `${PROJECT_ROOT}` 逻辑路径：Python读取冻结pkl、历史TSV和只读HDF5/Zarr schema，SHA核对后聚合；软件版本由importlib.metadata与Rscript只读查询。大输出留本地私有取证目录。

首个JSON序列化遇到pd.NA字典键，修正递归类型转换后完整成功；一次SSH连接超时，独立验证跳板和目标后重试只读取证，没有重启科学作业。全量JSON传输耗时不能当成远程计算仍在运行；进程只读检查证实输出阶段无对应Python计算。

## 当前边界、下一步

尚未完成新注释，不把历史阈值命中或代码存在视为发表级验收。当前节点正在建立全谱系原因表和表达核查；随后先对T/NK/ILC与稀有谱系检查父身份、同细胞共表达和数据可测性，必要时在谱系内重新计算HVG/PCA/批次校正/图/分辨率及RNA marker。保留大类、marker描述群和状态是允许的科学终点；父类未确认须继续标记review。

验证包括计数守恒、checkpoint指纹、精确index对齐、缺失/未测基因分开、旧标签保留和新增输出读回。公开前检查隐私和diff，运行仓库轻量检查，普通push并独立核验。

建议commit：`Record frozen V7 literature-guided annotation audit checkpoint`。

## 阶段A/B实测完成与Mast初诊恢复

- A完成UTC 2026-11-05T02:32:35.577879Z（UTC+8 10:32:35）：106标签、26父类、143标签/旧cluster/方法交叉记录，4676787细胞守恒；359基因panel逐轴精确对齐，66056表达汇总及2576同细胞程序记录。原checkpoint SHA未变，新增标签数0。
- 实测481271细胞的requires_recluster为true，不能仅改名视为完成。旧strict_major字段只作历史判定证据。
- B完成UTC02:36:28.849535Z（UTC+8 10:36:28）：全部45036个V7 Mast候选、95038基因、57129234 nnz，精确cell index对齐；raw H5AD SHA f2f57e8842c010b81b8b8d59eadd17195831d149e1abd275b5bff2f544cc8124。矩阵、细胞metadata及PID留服务器私有。
- Mast_0/3两基因Mast身份组合阳性98.8%/99.2%，TPSAB1/CPA3/TPSB2表达强，支持先按Mast父身份复查。coarse_provisional与pending群QC/污染分布不同，不能自动一并转成成熟标签。
- 首次Seurat初诊真实exit1（UTC02:38:06.510378），CreateSeuratObject将下划线转连字符后出现重复行名，未进入聚类。失败目录和代码保留。v2新增显式唯一feature ID↔原始gene symbol映射，不丢基因、不合并表达，另开stage_c2；R parse通过、8线程真实启动。
- 本地系统Python缺pandas，改用已有bundled Python处理聚合结果；中文stdin编码失败改用UTF-8 patch写文件。均不涉及科学结果变化。

### Tcell_6算法原因进一步纠正

旧逐基因证据显示FOXP3/IKZF2/IL2RA均通过、Treg全6个marker通过。旧winner的canonical=0不能解释为整群没有Treg锚点。旧程序对每个panel的通过基因取平均，先选最高分panel，再检查该winner：TIGIT单基因的CD8_exhausted均分0.47666，高于6基因Treg均分0.32886，于是选择缺锚点的CD8候选并拒绝，没有继续评估完整Treg候选。

此前“缺经典marker”应限定为旧算法选中的候选，不能泛指所有候选。当前T6 CD3D93.5%、FOXP3 46.2%、IL2RA38.7%、CTLA4 64.9%，应优先检查CD4/Treg富集和内部混合，不整群强标CD8耗竭。此为候选排序限制；后续逐panel审查和谱系内分群。

E阶段于02:45:24.695267 UTC / 10:45:24.695267 UTC+8完成：842可评估样本中841支持Treg panel相对同样本其余T细胞均值差>0.10，沿用旧支持口径作描述性复查而非显著性检验。57042细胞至少两个Treg经典锚点阳性，53819同时至少两个T身份marker阳性。全基因候选CellTypist Regulatory T cells=142908/156631，scATOMIC T regulatory cells=85227；旁证范围与粒度不同，不能以模型票数或cluster富集整群赋Treg。没有新标签改写。

首轮公开节点：A/B/D/E已完成；C2真实HVG/PCA/邻接图/8分辨率完成，RNA DE仍运行。公共文件只包含汇总、方法和通用代码；新文件显式复制到独立干净检出，已有并发改动保留。轻量仓库验证0 errors/1既有warning；提交、普通push及独立SHA核验另记收据。

## 当前终端与首轮判读完成

C2真实exit0于02:47:11.432208 UTC / 10:47:11.432208 UTC+8；F review完成02:52:12.736680 / 10:52:12.736680。95038基因原始计数做完整文库归一化；2000谱系HVG、50PC、初诊30PC、8分辨率、RNA Wilcoxon，未因未知Technology使用虚构批次。27群只作初诊分区，来源效应和源基因名差异显著，不发布27亚型。49个实际Mast来源header重读并修正marker可测分母，低于80%覆盖图示×。未完成source稳定ID/别名统一，不能据源特异DEG定义新亚型。

27群人工工作性review：20群/38350细胞父身份支持候选；7群/6686细胞竞争程序或QC待复查。不是通过最终QC、排除名单或成熟亚型。cluster10富集T程序且doublet中位0.110，cluster15髓系程序及doublet/decontX升高，cluster13/25/26有间质/上皮竞争；先保留复查。cluster16共748细胞，TOP2A/BIRC5/RRM2等增殖程序与Mast身份共存，记录状态候选。私有cell review overlay包含旧标签原样保留和新review字段。

独立算术复核首次在最后来源计数断言失败：Pandas categorical value_counts保留未使用水平，与R实际非零来源数定义不一致；837条marker重算在此前已通过。修正为实际字符串来源计数后v2于02:55:52.539496 UTC / 10:55:52.539496 UTC+8 PASS：45036细胞、57129234 raw非零项全测量掩码、0个非零落在未测基因、837条marker记录、27群来源/计数、旧字段和原V7 SHA均核对。生物学路线仍属review，不能把算术PASS当最终身份验收。

公开发布前发现远端有并发新提交29c5f8d，保留原隔离检出的改动，另从该最新main建立新隔离检出。首次受限网络无法补取promisor对象，失败目录保留；授权网络下建立v2新检出成功。不覆盖并发首页、当前工作或其他项目记录。新图已查看，标签及marker可读，无患者/样本行。

早期未压缩aggregate回传最终SSH exit255、只收到部分stdout；它是传输失败，没有重启科学作业。独立压缩回传已经成功，解码JSON和所有计数/输出后才用于报告。此前pending段落是准备同步的本地快照，未在当时push；实际公共发布以最终32文件节点及独立远端SHA收据为准。公共TSV经聚合JSON重建，服务器原始artifact SHA作用于逻辑run路径；公共文件字节以本次Git blob为准。32文件人工diff/隐私扫描未命中内部路径、账户或凭据，约0.95MB；仓库完整轻量检查0 errors/1 warning，既有并发内容保留。


Publication verified 2026-11-05T03:02:32.127141+00:00 UTC; content commit `4fa722cdcb3c535997d2c3d267cc986b05ded957`, 32 reviewed files, ordinary push actual0 and independent ls-remote exact full-SHA match. No concurrent work replaced. Final repository validator: 0 errors / 1 pre-existing warning. Diff whitespace check required explicit NA for unmeasured trailing TSV fields; 837-row numerical readback unchanged. Whole-project final annotation remains incomplete; first-stage review and its numeric audits are complete.

## 继续：marker图对应类型标注与Mast逐细胞复查

用户要求marker图同时标出marker对应的亚型，并继续注释。2026-11-05 04:04 UTC / 12:04 UTC+8开始本节点；原15个tracked并发改动保留。新图将同时标列分组（父身份、候选亚型、状态和竞争谱系）与行的工作性身份/复查标签，原始图和F输出不覆盖。不能把marker的文献关联自动当作本数据已确认的亚型。

已重读实际run/method/task/daily、当前项目首页和任务状态。补充读取泛癌Mast原论文作者pan-MC Fig1/2代码、JCI2025 MCT/MCTC原研究及JEM2023跨器官Mast研究，以公开原panel标注而非凭空命名。代码SHA留本地来源收据，公开记录只保留相关解释与引用。

服务器输入合同明确：全基因X为raw source counts；obs中的decontX为污染指标，原生全基因层未做decontX，也不存在可直接复用的全基因校正层。不得把指标当成已校正表达。下一阶段用原生marker共表达、样本内doublet分位、深度/污染、来源复现逐细胞细分7个review群；保留全部细胞，既不整群自动删除，也不以一个marker强行命名。

### Marker目标图与逐细胞QC实际完成

G于04:12:11.277769 UTC / 12:12:11.277769 UTC+8完成：60基因/1620 marker-cluster记录、49来源；完整core的1322个原样本/方法背景组用于分位校准，Mast中44818细胞可评估。初始任意2/5 Mast marker规则得到565联合候选。H进一步要求一个tryptase（TPSAB1/TPSB2）加一个CPA3/KIT/MS4A2，均≥2 UMI，防止同一功能家族重复计为独立身份支持；42578细胞满足，41191无检测到的强竞争程序，559为联合doublet候选，828竞争程序复查，2365低marker复查，93为弱Mast但强竞争谱系复查。全部45036保留，没有最终排除或全项目最终标签发布。

MCTC相关CMA1/CTSG阳性候选在raw≥2 UMI口径为1312细胞/40来源，MCT相关panel为1097/33来源；是非互斥的RNA marker关联候选，不视作蛋白定义亚型。pan-MC的C9参考panel与细胞毒淋巴程序共享基因，因此额外并排T身份panel，不能仅凭CCL5/IL32/GZMB赋Mast亚型。

大图单包回传网络超时，科学producer均未重启；跳板独立登录成功。改为14KB块，逐块与整包SHA核验成功后本地重绘；聚合包SHA d9f78d346a791065e2ea99eef603a6f985fb95ba455860fd7685e5d495ee3fd8。3张v3图均已视觉检查：列有目标类型/候选亚型/状态，行有工作注释和n，marker仍保留。独立临时绘图环境Matplotlib3.10.7，不修改全局环境；初次本地缺库及py_compile用户缓存不可写已通过临时环境和AST检查解决。

H输入准备04:15:06.569299 UTC / 12:15:06.569299 UTC+8；36007细胞进入保守QC敏感性模型，17136个来源可测性模型候选基因，RNA完整95038基因保留。H真实R exit0于04:20:38.264302 UTC / 12:20:38.264302 UTC+8，marker已生成。通用supervisor的stage字面沿用initial，R代码SHA/专属目录/实际complete均为H，以H实际manifest为准。Tcell_6全156631细胞原生输入提取已实际启动；尚未据此声称Treg细分完成。
