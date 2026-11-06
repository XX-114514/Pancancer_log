# V7授权清理与谱系重建

用户于2026-11-06明确授权清理，要求沿已讨论的QC/注释逻辑继续排查。开始取证00:22:18 UTC / 08:22:18 UTC+8（Melbourne11:22:18）。本轮实际创建新的QC决定、排除/暂存名单与分析视图，随后从counts重新计算谱系HVG/PCA/图/多分辨率和RNA marker；不覆盖原V7，不要求再次确认。

初始状态：冻结4676787细胞V7及原核心对象保持，54来源历史过滤已实证，28999个not_run是真正无分数；T6有845个联合候选（K-09其中245），Mast有1323个联合候选。此前仅打review标签，未作为新增清理名单执行。主工作区15项tracked并发改动和其他未提交内容保留，公开同步使用独立检出及显式文件。

已读：仓库规则/README/STATUS/ROADMAP/TODO、原doublet核验、Mast/T6文献复查及panel修复方法；服务器原祖先规则只读核对。新科学run为`${PROJECT_ROOT}/runs/20261106_v7_QC_clean_reannotation_v1`，原annotation_rebuild输入只读。没有删除旧文件或修改生成的全局状态指针。

## 预注册处理与边界

- 本轮先做Mast/Tcell_6的实际清理重建，同时补查V7全部既往未评分记录单元。范围外谱系保留旧QC状态并注明尚未进行本轮完整复查。
- 多谱系程序+高原始背景内doublet分位作为入口，再审阅独立谱系锚点。锚点支持的联合候选进入新主分析QC排除；锚点不足的暂存复核，不在本轮主分析中硬赋亚型。此为操作性QC决定，不称生物学真值确认。
- 不按固定比例删除，不把未测marker或not_run当阴性。Mast父类支持与QC分别记录；不满足Mast精细分析入口者保留以待父类/环境RNA复查。IFN/增殖等状态不作为删除理由。
- 未评分单元按完整记录样本单元建模，不把Mast/T6子集直接当作原文库做Scrublet。低于100细胞保留不可评估；≥100的单元用Scrublet0.2.3、expected_rate=.04、seeds0/37/20261106、min_counts2/min_cells3/min_variability85、PC=min(30,max(5,n/5))、精确近邻。三种子成功、可检测模拟比例≥.1、各调用率≤.15且调用率极差≤.05才进入阈值/直方图审阅；3/3同细胞阳性才是补充排除候选。稳定性门槛为项目保守复查条件，不是论文通用概率或15%删除配额。
- 原方法建议核对模拟分数直方图和嵌入；依据[Scrublet原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/)、[作者实现](https://github.com/swolock/scrublet)与[Scanpy官方文档](https://scanpy.scverse.org/en/latest/api/generated/scanpy.pp.scrublet.html)。数值成功不自动代表可用，阈值审阅单独完成。

## 实际阶段

00:24:18 UTC / 08:24:18 UTC+8：原阶段终端、软件和资源预检完成；Scrublet0.2.3、DoubletDetection4.3.0.post1、scanpy1.11.5；Rscript存在，可用内存约840GiB。新run拒绝覆盖旧目录。

A于00:26:31 UTC / 08:26:31 UTC+8完成逐细胞锚点复查：T6的845联合候选中783有不同组件锚点支持、62不足；Mast1323中1150有支持、173不足。K-09的245中227有支持、18不足。基于全基因raw counts，不以marker总数/单一共享基因冒充谱系身份。该阶段输出证据，没有在A阶段自动剔除。

B补测已使用独立detached supervisor真实启动，最多4 worker、各1线程，写新私有输出与真实terminal，不创建定时推送。补测拟合、直方图审阅与最终名单尚待实际终端。下步固定已审阅规则，生成主分析排除/holdout/保留视图，再重新建模。

完整细胞行、样本标识、矩阵、PID和stdout留服务器/私有本地。公共产物仅审阅后的汇总、方法和通用代码。验证包括原SHA/index保留、全量计数守恒、排除/holdout/保留互斥、模型输入精确对应keep名单、原始计数拷贝一致、可测性分母、独立重算与真实退出。失败恢复和发布收据在本日志追加。

B实际terminal0，96单元×3种子均拟合完成；27单元稳定，69单元阈值/可检测性待审，19单元<100。17阳性单元直方图均检查；34采用、1数值警告暂存。C00:39:54 UTC实际完成1967排除/236暂存/4674584保留；独立重算名单、输入CSR及10x逐元素通过。发现初版聚合holdout字符串错误，新审计表纠正，不影响实际mask/input。D于00:42:18 UTC / 08:42:18 UTC+8真实启动两谱系，各16线程，R语法解析通过；终验尚待。

公开检查点只复制本任务新增文件，在origin/main独立worktree上更新索引，保留原工作区15项tracked并发改动及其他untracked。检查代码AST、公共汇总无cell/sample行、diff及repository validator，普通push后独立核SHA。建议commit：Record authorized V7 QC cleanup and clean lineage refit launch。

00:46 UTC公开push尝试被自动审批拒绝，理由为此前40+8文件授权不覆盖本轮25文件的私有分析衍生汇总及脚本；未执行网络push。已请求对已审阅本地commit84adfc6的明确授权，科学分析继续。普通推送仍待授权，不另行绕行。

00:48 UTC用户明确授权公开25文件。普通push成功，独立ls-remote核到84adfc683cb6727e784179e0af1746d41b7fed6d。此前自动审批拦截已由本批明确授权解决。后续新结果未包含在该25文件提交中。
Mast重拟合实际00:47:37 UTC terminal0，38959细胞；一组FindAllMarkers未过min.diff.pct，保持待marker解释。E首次汇总已写表但manifest因numpy.int64 JSON失败，不作完成；新v2目录转为Python int并重新汇总，保留失败历史。

Mast扩展NK三组件后36个校准有效的联合候选于01:02:17 UTC实际排除，v2全局2003／236／4674548，独立重算/原始计数一致。两个当前模型于01:09:01／01:09:07 UTC正常terminal0；不同版本编号不继承身份。先读完整RNA DE再制定工作标签；T6补Ig面板到118基因，Mast110，共3396条marker独立复核。当前无正向DE分区T6 0/7/10/11、Mast1/5/7；T6有10条Quick-TRANSfer和2细胞DE跳过，明确未接受最终细分。

大类标签按簇多marker和来源人工复核，不把逐细胞raw2阳性作为每个细胞的大类必需门槛。H于01:19:39 UTC / 09:19:39 UTC+8创建201667行两个原候选池工作侧表并读回完全一致；T6大类T154974、混杂662、低信息150；Mast模型38923，大类支持但其中2150上下文复查，另4754原候选父类/竞争RNA待复核。未改canonical V7，未自动替换CoVarNet输入。

图全部显式marker对应类型／候选状态＋行工作注释／n，统计使用来源可测分母；10张图SHA、SVG XML及全部行标签读回通过，人工查过代表性身份、状态与扩展上下文图。E首次numpy.int64序列化失败保留，v2/v3新输出恢复；初次公共检查点命令一次相对cwd错误已改绝对路径，未影响分析。Pandas分片PerformanceWarning不影响写出/读回，记录不隐藏。

当前公共结果包仅脱敏汇总、方法、通用代码与PNG；所有细胞/样本行、矩阵、PID和完整日志留私有。公共检查点25文件已授权并核远端84adfc683cb6727e784179e0af1746d41b7fed6d；其后结果属于追加内容，准备具体清单供授权。未把本地提交当云端同步成功。验证：全部脚本AST、逐表数值/守恒/范围、最终原V7 SHA、图读回、diff/privacy检查及repository validator。建议commit：Complete V7 QC v2 and marker-guided Mast T6 working annotation。

最终公共检查实际结果：repository validator 0 errors／1个既有历史路径warning（2026-10-13旧provenance）；本批88文件privacy/schema/AST检查通过、无细胞/样本行或矩阵、diff --check通过。GitHub原25文件提交已独立核验；本追加79新文件＋9更新共88文件完成本地审阅，尚待本批明确公开授权。
