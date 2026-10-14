# HNSCC固定MCP临床与配对验证

## 要求与起点

持续推进原9CM的生物学/临床及细胞协同意义。上一节点完成用户原Table S1接入与精确GEO患者/取样键核对；二分类标签存在冲突，旧程序完整测量待，原权威身份/9CM及并发atlas/state任务保留。

## 实际执行与审计时间

- 2026-10-14 03:53:50 UTC /11:53:50 UTC+8：基因轴v2稳定ID存储别名核查完成；唯一可测物理列恢复，歧义及碰撞不猜。原权重/95%门槛未变，完整分数仍缺失，可测部分单独登记。
- 03:57:14 UTC /11:57:14 UTC+8：临床定义、原生RNA固定投影及CRC发现/留出校准完成；此时未算临床效应。使用保留原生counts及全来源RNA分母，不以ambient-corrected X代替冻结CRC输入口径。所有原权重、标准化、质量系数未重学；最多50cell/库/端点的固定随机抽样独立于结局。
- 主检验v1在任何效应计算前因治疗字段合并后缀退出；v2逐行确认两个字段一致后恢复，输入分数不重算、门槛/患者/检验合同不变，旧失败保留。
- 03:59:27 UTC /11:59:27 UTC+8：全部主系数/t值、20000残差置换、2000块bootstrap、辅助效应与各6槽BH独立复算PASS。
- 04:01:42 UTC /12:01:42 UTC+8：原预声明6问题×2程序槽配对家族完成，全部可估计效应、置换、区间及12槽BH独立复算PASS。
- 所有组成/捕获/留一患者描述变体及原生RNA质量/年龄描述变体实际完成；不追加P/BH、不替代主分析。完整科学报告保留服务器，真实统计仅获准当前对话审阅。

## 临床口径与范围

主要结局为作者连续Tumor.regression，同患者各取样一致；二分类原标签保持不动。按论文病理退缩>=50%另建明确标注的派生标签，作为辅助；只有作者标签跨访视一致且与阈值一致的患者进入作者标签敏感性，其响应分组覆盖不足时不做检验。未将冲突直接写回权威registry。

基线原生程序与连续病理退缩的患者级秩模型调整D/D+T和HPV；至少20患者、至少10残差df、满秩，不自动删协变量。3端点×原2槽=6主检验；原MCP2缺基因证据，未产生新MCP2，内部P=1参与校正、展示仍缺失。辅助连续/派生二分类/作者一致标签家族分开。不使用治疗后分数作基线预测，不使用RECIST或生存口径。

新atlas端点为多基因CD8与巨噬亚型候选；尚无旧172细身份直接映射，因此不是原9CM完整临床复现。来源内共变、程序临床关联、直接互作/因果机制各自登记，不能相互替代。

## 科学结论边界与后续

当前出现值得进一步追踪的巨噬端临床关联线索；T端/两端合并和配对协同结论并未同时获得支持。全阴性/缺覆盖/不支持槽位保留。组成、原生质量、年龄和留一患者是事后描述性诊断，不是新增显著证据。不能以巨噬一端关联命名原CM为耐药模块，也不能宣称已证实T—巨噬直接通讯。

下一以来源限定的程序临床线索继续：审病例别名/临床协变量与身份内状态贡献；在完整冻结候选家族下与发送端及接收响应联系；寻找满足RNA/患者/端点的独立临床或空间验证。当前空间靶向panel不满足完整程序权重覆盖，不降低门槛。整体目标仍active。

## 可复现输入与产物

服务器环境scanpy Python；运行代码通过PANCANCER_PROJECT_ROOT/PANCANCER_PREFLIGHT_ROOT配置。运行check_HNSCC_frozen_gene_axis_v2.py、prepare_HNSCC_frozen_MCP_clinical_v1.py、test_HNSCC_fixed_MCP_clinical_v2.py、audit_HNSCC_fixed_MCP_clinical_v1.py、test_HNSCC_fixed_MCP_paired_v1.py、audit_HNSCC_fixed_MCP_paired_v1.py、diagnose_HNSCC_clinical_composition_v1.py、diagnose_HNSCC_native_quality_v1.py、report_HNSCC_fixed_MCP_private_v1.py；各生产/审计终态均exit0，失败v1另保留。

逻辑输出：${PREFLIGHT_ROOT}/HNSCC_frozen_MCP_gene_axis_20261014_v2、HNSCC_frozen_MCP_clinical_preparation_20261014_v1、HNSCC_fixed_MCP_clinical_results_20261014_v2、HNSCC_fixed_MCP_paired_results_20261014_v1、HNSCC_clinical_descriptive_composition_20261014_v1、HNSCC_native_quality_descriptive_20261014_v1、HNSCC_fixed_MCP_private_evidence_20261014_v1。每版含输入/代码SHA与冻结计划/终态。主统计产物SHA：6a0fbc7c37e93f8cd9eb3f2cc351867ffc31c9abd8e7942e65efe2de04b587b8；配对统计SHA：44e232069376fb70e811f7d4b34fde625f3028421c2f90870104ef6bcb2d67a4。

无新表达下载、数据删除或历史覆盖；用户原表/患者/细胞行、科学统计、基因/权重不进公共GitHub。代码与方法/执行轨迹即时普通push，云端提交另独立核验。建议commit：analysis: complete fixed HNSCC programme clinical and paired validation with private evidence。
