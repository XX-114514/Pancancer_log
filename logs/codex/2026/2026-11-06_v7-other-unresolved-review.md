# V7其余unresolved逐类复查

用户于2026-11-06要求继续分析其他unresolved，并授权本次提交。开始取证03:43:05 UTC / 11:43:05 UTC+8（Melbourne14:43:05）。此前88文件已普通push，独立ls-remote核验2c19567e075b1251cd0a109c03d178eb9cf8921d；本轮科学工作另建版本，不覆盖前轮Mast/T6或冻结V7。

初始科学状态：全4676787细胞QC-v2视图2003排除／236暂存／4674548保留；Mast/T6工作侧表和marker已读回，细分未接受。主工作区既有tracked/untracked并发工作保留，本轮公开副本从2c19567独立检出。

初始范围：先复查Erythroid、Platelet、HSPC、ILC，联合造血分化边界、T_NK_ILC边界及原Ambiguous/Unknown作召回入口。旧父标签不是终判；先确认测量范围、大类组合及竞争身份，再处理有多组件证据的QC候选，从counts重建，最后输出大类、状态、命名marker与原因。HSPC_like不证明正常HSC；IL7R/GATA3与成熟T共享不能当ILC确证；红系不只凭HB，血小板不只凭PF4/PPBP。允许大类或大类+marker背景，不为每个计算分区强造亚型。

已读：仓库AGENTS、README、STATUS、TODO、ROADMAP、RUNTIME_SYNC、前轮task/run/侧表验收、106标签原因表、panel修复方法及registry。服务器预检只读核验现有mask、基因和细胞轴、上级规则及资源。源矩阵、细胞/样本行和完整日志留私有；公开仅审阅汇总和通用代码。

拟验证：固定V7索引hash及当前QC继承；全候选测量分母；独立组件同细胞证据；父类竞争与生物学相关priming分开；清理名单互斥和counts精确子集；真实terminal；所有marker统计独立复算；图标对应身份/候选状态及每群n；工作侧表保留旧字段读回一致。未完成者不写成功，细分与大类接受分开。

失败/纠正：首次读取reference表误用reference_id，实际列为id，已按schema修正；未影响科学计算。后续实际参数、结果、失败、验证、尚未完成项与发布收据持续追加。建议commit：Review remaining V7 unresolved rare lineages with source-aware markers。

预检03:47:59 UTC确认轴4676787×95038、QC-v2继承、可用内存832.7GiB／磁盘3666GiB。A于03:53:23—03:56:32 UTC / 11:53:23—11:56:32 UTC+8真实运行并terminal0：298674细胞、213基因／35模块；旧父类外候选Ery78／Platelet80／HSPC16／ILC403，只是召回入口。

B成熟谱系配对复查03:59:48 UTC结束，915联合候选（913拟操作性排除、2血液肿瘤语境holdout）；已逐对审阅来源、原样本分位和各侧表达份额。T+NK/T+ILC、B-ASC过渡及HSPC priming均不自动删；HSPC对内皮比较须SPINK2/HLF＋CD34/PROM1＋PTPRC，不只凭CD34/MEIS1/GATA2。血液背景regex采用边界约束，避免把Gallbladder误当ALL。原QC decision中的04:01:00为准备记录的名义时刻，精确人工审阅时刻未独立捕获，实际执行时间以producer UTC/SSH退出收据为准。

C于04:06:13 UTC / 12:06:13 UTC+8真实完成913排除、2暂存增量，全球QC-v3为2916排除／238暂存／4673633保留；四输入Ery9098／Platelet2318／HSPC14072／ILC45290。未删旧文件或覆盖V7；独立输入审计另存。当前只称QC/输入完成，不称大类或细分已接受。日志补写一次UTF-8解析错误，改用patch补录，不影响科学输出。

D四模型04:08:43 UTC启动；Ery04:10:42、Platelet04:09:30、HSPC04:11:21真实terminal0。E首轮233 marker显示多类父标签误标，追加肝/腺泡/纤毛/髓系面板作为E2新证据，不覆盖E。独立审计Platelet/HSPC初轮通过；Ery初轮CP10K复算失败，定位审计表达式10000乘int32溢出，生产统计本来使用浮点缩放。修正审计为10000.0，保留失败历史，再独立验证E2。

ILC于04:17:21 UTC / 12:17:21 UTC+8真实terminal0；四模型共70,778细胞、59个res0.4诊断分区，不能当59个生物学亚型。四类E2均260个不同marker，15,340条记录；04:17:49—04:18:44 UTC全部独立复算PASS，最大CP10K绝对误差2.28e-12以内，所有分辨率/来源计数、HVG与DE共同测量资格一致。无正向DE群：Ery1/3、HSPC0/1/3、ILC1/2，保留而不强配细分。R locale及min.diff.pct警告保留。

逐群核读RNA、正向组件、DE、来源集中度、旧CellTypist/scATOMIC辅助表与疾病语境后，写人工cluster_decisions；04:23:00 UTC / 12:23:00 UTC+8实际生成四类及合并工作侧表并读回；04:23:41 UTC独立侧表审计PASS。70,778行与模型细胞轴精确相同且互不重叠，旧V7标签/父类保留；27,074大类工作记录，43,704身份待复查。冻结文件SHA仍e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e。身份侧表不自动改变全球QC，QC保留不等于身份通过。

关键观察：Platelet池1,479个血管壁簇级工作记录、149个Platelet/MK谱系群；ILC池6,363个正向T路由及496个黑色素谱系簇级记录，保留逐细胞竞争复核。HSPC_like簇12,546个，其中8,937血液疾病语境，不能确证正常HSC/MPP/GMP/MEP。Ery12的189个有强非HB/结构＋cycling背景，另605个从待定分区经非HB组件进入红系工作路由；不是全数据只有189个红系。原HB/NRGN/KIT/CD34等共享或环境RNA可能性均不能代替父身份。

G于04:31:30 UTC / 12:31:30 UTC+8追加166个ILC-like候选组件核查与新侧表读回：32来源、67类型2/99类型3正向RNA程序，无两类同时满足者；未接受ILC2/ILC3，未用TCR零值证明Lin阴性。背景226,981个入口尚未做身份重建；未宣称所有V7细分或其他谱系已完成。

本轮文献核查沿用已读Lodi作者代码与registry，并正常访问Björklund原始PubMed摘要、Velten可获取原文信息；Nature全文访问失败、PMC直接页验证限制，不宣称未获取全文已完整核读。具体阈值是项目诊断，不是论文统一分类器。完整论证见方法及结果。

命令概括：经已配置跳板以BatchMode SSH执行固定环境Python/R，A/B/C顺序生产与独立输入审计；四R新目录并行16线程；真实terminal0后E/E2和独立marker审计；人工决策后F及独立侧表审计；G只做候选组件跟进。完整命令和私有启动收据留工作区；公开脚本使用显式PROJECT_ROOT与逻辑RUN_ROOT。运行命令、参数和版本见正式run。没有重跑原全V7 doublet工具或新TME任务。

补充失败/纠正：一次大红系只读汇总SSH超时且只收到部分字节，未采用不完整输出，另一次只读取证成功；没有重跑科学任务。几次本地读取在等待文件生成前尝试导致missing-file，以及默认Python无pandas，改用已验证绘图环境；不影响服务器科学产物。工作侧表草稿在执行前发现字段名不符，另存v2再执行，未运行旧草稿。两次文档patch上下文不匹配，工具原子拒绝无改动，按真实上下文修正。所有失败历史保留私有，不将传输失败当科学失败。

交付修改：结果目录包含59群完整中文原因表、每类260 marker/组件/原分类器/DE/来源及分辨率汇总、工作身份聚合和旧新对照、QC及旁表审计、逻辑视图合同、policy、人工决策、PNG/SVG各16；新增通用生产/审计/绘图代码、方法、run、daily及run/artifact索引，并更新README/STATUS/TODO/ROADMAP/CHANGELOG。没有上传细胞/样本/患者行、矩阵、主机地址或凭据。

验证：四代表性科学PNG实际查看（红系身份、血小板竞争、HSPC身份、ILC竞争），16个SVG逐行标签/计数/图注及全部图SHA由交付检查验证；marker对应身份/状态面板显式。真实数据交付检查PASS：59群/70,778行汇总/15,340 marker/16图，manifest hash、TSV无细胞/患者/样本行key，新增Python AST通过。仓库验证首轮0 errors/1既有路径warning，候选新文件专项地址/凭据检查无命中；最终合并及提交前继续核验。

同步发现远端main已推进到d130faae（并发V5通讯与TME来源决策），将正常合并保留其提交及状态章节，不重写历史或覆盖主工作区并发修改。此前88文件已同步2c19567，本批审阅后按用户授权与AGENTS长期授权普通commit/push并独立比对远端main；最终SHA见交接收据。

未验证/后续：RNA污染机制、蛋白/Lin阴性、HSPC正常性/恶性与功能、来源/参数稳定性、剩余背景及全V7最终细分均未接受；目前仍需进一步来源内子集重建与相容参考。建议commit：Review V7 unresolved rare lineages with actual QC and source-aware working annotations。

发布合并：本轮160文件科学提交f1a21e16已本地保存；普通合并origin/main保留c63e3ef3与d130faae。仅STATUS顶部发生内容冲突，合并保留V7四类、并发V5通讯及旧QC-v2全部章节。首次stage检查发现完整表EOF多空行，因组合命令末段返回0未及时阻止本地commit；在合并提交中规范EOF并刷新artifact hash，保留原提交历史，再单独检查退出码。该修正不改变科学数据。

最终检查纠正：合并后验证1 error/1 warning；error来自交付自检脚本里的密钥检测正则字面量触发仓库扫描，并非存在凭据。也发现该自检模式包含不必要的本机特定身份/主机匹配字符串，已改通用模式。为避免这些草稿字符串进入公共历史，f1a21e16及其合并仅保留本地审阅历史，另从最新远端d130faae新建public_sync_checkout_v2，复制已泛化、审阅通过的本轮160文件并普通提交；远端并发提交作为父历史完整保留。不删除或改写本地草稿提交，不推送含草稿的历史。科学数据与模型不变，重新验证发布副本。

发布副本v2最终验证：validate_repository.sh真实退出0、0 errors/1既有路径warning；无新FAIL。真实数据交付检查再次PASS，160文件的diff --cached --check退出0；14个新增代码文件逐文件SHA与交付清单一致，专项本机路径/身份/主机模式检查无命中。发布将仅包含从d130faae直接形成的新提交，不将本地草稿历史作为祖先；普通push之后独立ls-remote核验并将最终SHA写入私有传输收据及交接消息。主工作区既有tracked/untracked修改保留，仅复制本轮全新唯一结果路径供本机查阅。
