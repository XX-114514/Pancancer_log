# V7 组织上皮面板扩展与原始组合审计

用户要求核对 Claude 整理的组织特异上皮组合、补充原始文献并纳入当前 V7 面板。用户已说明原文件没有逐项或逐基因参考文献，不将其视为作者原始签名。

初始状态：已读取仓库规则、当前项目状态、V7 修复方法/资源及运行日志；112 模块的 v7_panel_repair_v1 已发布，但全 V7 细胞重新注释未完成。主工作区并发修改保留，在干净发布副本上从 main 7656d03db957501cf7768f3a3774d731ec5e09b4 继续。源文件 methods/annotation_gene_panel_v1.tsv 实际 100 行，其中 role=tissue 且 compartment=epithelial 的只有 26 行，不是 27 行。

计划：保留源文件，逐行映射到拆分后的候选模块；逐基因区分直接文献出现证据、此前已整理证据和未核实项；生成独立版本与兼容原引擎的合并 registry；用反例验证身份/器官/状态/良恶性分离；更新方法、运行与必要索引后正常发布。

已检查：源 TSV、旧 registry.json、gene_reference_trace.tsv、evidence_engine.py；肺/肝/前列腺既有参考及肠、胃、胰腺、乳腺、肾、膀胱、甲状腺、输卵管、子宫内膜、皮肤/食管原始研究。正文/图注确认的范围和定位随新资源保存，不声称全文补充表全部逐行核验。

初步发现：胃主/壁细胞、前列腺基底/club、近端小管/TAL、集合管/足细胞需要拆分；原乳腺管腔组合偏激素响应；肠 SMOC2 不能作为人肠干细胞独占证据；器官相容性不等于唯一来源，正常参考不能判定良恶性。

失败与边界：若干 PMC 请求返回 reCAPTCHA，部分 Nature PDF 跳转失败；改用作者机构 PDF、其他可访问原始研究及可核查正文。未确认的原基因保留待核实状态，不靠引用整篇论文冒充逐基因核实。未执行服务器全矩阵重注释，未改变任何细胞标签。

验证、时间、文件与发布结果在完成时补录。建议 commit：Add literature-traced epithelial panel extension and split mixed identities。

## 完成的资源与验证

2026-11-05 07:04:08 UTC / 15:04:08 UTC+8补录取证。新增scripts/v7_epithelial_extension_v1构建器与9项语义反例/字节重建测试；resources/annotation/v7_epithelial_extension_v1包含32项完整Markdown表、原121条基因逐项处置、JSON原始快照/映射/当前基因证据/合并registry与manifest。新增方法、run、daily，并为README/STATUS/TODO和方法/run索引添加窄范围入口。

构建结果：25新增+7复用=32目标，112旧+25新增=137总模块，473不同正向基因；119条目标模块×基因记录中85条本轮定位、34条继承旧trace。42条原记录未导入对应目标，未宣称全部原基因核实。只在合并副本调整旧近端小管模块，源文件及旧registry不改。

实际命令：build.py生成；test_extension.py测试；资源全量字节重建比较通过。首次出现两个未关闭输入句柄的ResourceWarning，已修正with后9/9通过。文献图号复核发现肾图谱基因标签位于Fig.3d p.589，已更正构建器及所有生成记录；原图像URL尝试404，未保存或发布伪造图像。

Spreadsheets技能用于源TSV只读分析；输出为代码资源JSON和仓库Markdown，未新建Excel。Pages write-page技能用于保留仓库既有目的地与文档体例。仍未执行生物学全量注释；42条未导入记录和未覆盖上皮类型需后续核实。下一步按新registry检查真实上皮谱系的测量、同细胞表达与完整DEG。

发布前检查：9项语义/重建测试通过；仓库检查SUMMARY: 0 error(s), 1 warning(s)（历史内部路径警告保留，非新文件）；任务15个新增资源/代码/说明文件的路径与凭据模式扫描0命中，JSON解析通过，20个暂存文件diff空白检查通过。原始人类细胞行、矩阵、凭据均未纳入。

## 并发TME v2兼容补充

2026-11-05 07:09:54 UTC / 15:09:54 UTC+8取证：首次普通push因远端并发提交3400d49而被Git拒绝（non-fast-forward），正常fetch/merge保留双方历史。新TME v2的179模块中仍有混合上皮项目组合，因此另加compose_tme_v2.py和185模块兼容registry：19个项目上皮模块替换为25项，近端小管更新，159个模块内容保持相同，41个Lodi作者代码/参照模块逐字段保留。原v2 JSON不改，README增加新版本入口。两项兼容/字节重放测试加入，总计11/11通过。

兼容registry SHA256为`c23a96b51a2d3394020ebef1c326cc1621e8a3ad23055878bebb0cbba9765cfe`；来源、替换映射与两个输入SHA见资源目录中的tme_v2_overlay_manifest.json。仅统一资源，没有切换正在运行的任务或更新任何细胞标签。

兼容补充后重新验证：11/11测试通过；完整仓库检查仍为0错误/1既有警告；新资源/代码隐私扫描及暂存diff检查通过。

独立发布核验（2026-11-05 07:11:42 UTC / 15:11:42 UTC+8取证）：内容及兼容提交`771024a5538ec9f64dd87a3a7b47cf9555d8f2b3`已普通push到用户指定public仓库。随后独立`git ls-remote origin refs/heads/main`与本地HEAD完整一致。发布包括原扩展提交1dd7b73、保留并发3400d49的正常合并及兼容修订；无force push或历史改写。公开资源完成不改变全V7最终注释仍待执行的状态。
