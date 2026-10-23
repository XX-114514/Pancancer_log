# 固定程序分布读出的既有覆盖审计 P0.1

本轮是[Pelka启发路线](../methods/Pelka-informed-program-distribution-roadmap-20261024.md)的首项真实输入审查。只读取旧验收的材料计数/来源上下文/患者块标志/程序基因轴资格及捕获池元数据；不读原生表达、基因权重、逐细胞分数或临床结局。个体键只在服务器内存用于去重，返回来源级汇总留私有。本轮尚未完整审查所有宽父类、外部状态字典和泛癌来源。

## 冻结输入与两次实际执行

旧机会producer SHA `2a37f7af2f73ac5ceddcb065626df411ca6f016271fce7f0830a1bd658f83cd1`，完整旧机会audit SHA `7f181b10ce21de7b84b0ebb0f34ca88a34946afcaabba38a0c48db2b24157e0e`。旧全捕获native producer SHA `2463d66299d7b1f1be96d99a19554a834735e688440716edfa00cd1aa01575e9`，对应独立audit SHA `f7a6db048cfccac2a98de84023ce5048d4faa595353d2589767718d0e661cfc0`。只加载这两个producer绑定的少量必要TSV，每文件最多1MiB，全部读入前及结束时校验SHA。

- v1：[元数据审计代码](../scripts/audit_program_distribution_coverage_RAM_v1.py)，SHA `708dbaece110a141fd422c50ff844156080248aaa19723e897ec8b23387b29ce`；计划SHA `695784cab2c107d273f3b190b89f199b225a9b020bae345159c4e1aa70cce01f`。实际源执行UTC2026-10-23T14:14:31.143037至14:14:46.636367，exit0；按完整解剖/材料/阶段上下文生成覆盖汇总。
- v2：[两层上下文代码](../scripts/audit_program_distribution_cohort_scope_RAM_v2.py)，SHA `eeea37a60e5956b1949fda11e86f587ccacee018bebe167117bca1f105f5e0b2`。在观察到解剖子部位切分后，先记录原因及合同，再追加同来源×癌种×材料类别×阶段×捕获×ex-vivo口径，省去literal tissue；保留全部原细分结果。UTC14:16:00.151525至14:16:17.192400，exit0。对应UTC+8为22:16，Melbourne为2026-10-24 01:16。没有根据关联P或临床结果选择聚合方案。

源使用现成scanpy环境Python3.10、仅标准库，`-B`禁止字节码。现有SSH、严格主机校验，只读执行，120秒本机观察限时，512KiB汇总返回上限。没有新服务器文件、下载或科学模型；源stderr仅保留字节数/hash。

可复算命令结构：

```bash
python -B scripts/audit_program_distribution_coverage_RAM_v1.py --check-artificial
python -B scripts/audit_program_distribution_cohort_scope_RAM_v2.py --check-artificial
python -B scripts/audit_program_distribution_coverage_RAM_v1.py --root "${PANCANCER_PREFLIGHT_ROOT}"
python -B scripts/audit_program_distribution_cohort_scope_RAM_v2.py --root "${PANCANCER_PREFLIGHT_ROOT}"
```

真实运行通过SSH标准输入传入代码；输入根路径仅私有控制记录。`__file__`为stdin时结果内代码SHA为空，真实代码SHA由事前本机合同、实际传入字节检查和退出收据绑定，不冒称源文件hash已读。

## 审计内容和边界

门槛固定20/50/100/200/500/1000细胞/双方。分别返回材料数、已审来源内病例块、重复材料、未知患者文库、旧MCP发现病例重合，以及完整程序基因轴标志；旧core-T/捕获准入保持。缺某一程序轴不影响计数诊断，但不能据该计数宣称该程序可测。病例数表示同上下文中至少一份材料合格的机会；正式单材料选择或重复材料模型尚未冻结，不把该机会数当已确定的训练队列。

人工检查覆盖重复患者、前后阶段、未知患者、缺基因轴、无合格材料、错误布尔/负计数/重复材料拒绝。v2额外构造不同解剖部位的同一患者及血液材料，确认合并部位不会增加患者数或合并血液。真实每个上下文×门槛使用独立直接Counter循环复算；与既有native捕获池的全细胞数、100细胞上限逐条一致。

本机独立验收所有计数分解/门槛单调性、v2原细分结果与v1完全一致、代码Python3.10 AST。新源输出SHA v1为`a63a43a556b5862861ae89ab46943abef6089203471b8d74ef476589c0f13d28`，v2为`44d02af49401dd8837e36a9dd0e249763198a74689530890f6543103ac354a23`。科学来源计数、排序、解读仅私有/获准对话。

旧native完整测量文件保存全捕获均值及至多100个所选细胞的逐细胞分数；求和/均值不支持恢复全体尾部分布。P1需要现有计数的全细胞流式重算及RAM内分布，不得把旧抽样当完整分布。基因可测、细胞覆盖、功能状态、来源外验证各为独立资格；材料数不当患者数，跨来源重复患者未新审。

## 下一阶段

先按测量覆盖选择P1试点、固定宽父类与RNA捕获的不同问题、既有程序/外部参考及稀释拆半诊断。旧发现病例重合不因新读出消失；相关统计和临床家族仍未启动。其余P0字典/泛癌来源资格继续，不能将这部分审计当完整研究目标达成。主日志与仓库检查见[任务日志](../logs/codex/2026/2026-10-24_Pelka-design-and-program-distributions.md)。公共push此前拒绝保持。
