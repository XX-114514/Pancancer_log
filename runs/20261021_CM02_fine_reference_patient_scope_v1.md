# CM02细参考患者覆盖与计数版本预检

[总主线](../methods/CM02-specific-states-spatial-clinical-line-20261021.md)、[前阶段](20261021_CM02_state_context_v1.md)。对象为2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor），原成员/W/S/H/正式注释保持。

## 冻结问题

全部七CD4×六CD8作者细参考标签，pre/post共84组合；post2只作覆盖。唯一细参考文件父类各至少20、患者/时点唯一、原治疗/HPV完整后核可检验范围。此分母来自处理后作者参考子集，不是全组织/CD4/CD8真值；不足/零子群不反向补标签。

所有当前节点按阶段/药物统计原病例覆盖及唯一参考/交叠/无参考守恒；不能由汇总节点覆盖推断每个患者同样覆盖。原宽Teff对六CD8标签全部63个非空子集，分别在原作者、当前匹配及去交叠三个范围作严格逐文库计数一致检查；其余同名标签另核。按来源/阶段汇总，不输出病例行，不由相似或最小误差选标签对应。

## 执行与验收

使用原scanpy Python3.10、既有R解码器、单BLAS；原项目只读，小输出新任务tmpfs。PANCANCER_PREFLIGHT_ROOT为已审作者表/原计数，PANCANCER_PARENT_OUTPUT为已审参考桥；PANCANCER_RAM_OUTPUT是新任务目录。独立核查另以PANCANCER_PATIENT_SCOPE_OUTPUT指本预检输出。官方细元数据固定SHA，只白名单cell key/Patient/Time/state在RAM；不解析TCR/临床结局值，不下载表达或新病例表。

```text
python -B scripts/qualify_CM02_fine_state_patient_scope_v1.py
python -B scripts/audit_CM02_fine_state_patient_scope_v1.py
```

第一版派发器远端文本SyntaxError发生于解释前，未创建或运行科学任务；原失败保留。派发器v2使用原样多行文本并在发送前语法检查，科学源码未改。科学终态、所有范围和版本机会复核、具体科学解释分别登记，未到终态不称完成。两个原始/匹配计数范围共享固定输入矩阵，独立审核须披露此边界；唯一计数、病例资格/精确分数变异、覆盖分位数和原作者计数使用其他实现复核。

公开Git仅通用代码/方法/技术SHA；来源计数及科学解释保持已授权私有，原拒绝公开push/随访表未重试。整体互作目标active，尚未发现经状态—空间—同终点验证的关系。


## 完整细状态描述阶段

病例/版本预检与复核终态后，新描述于2026-10-21T00:44:55.102325+00:00 UTC启动。七CD4×六CD8×pre/post×来源/命名锚点两范围共168位置，源源码/随机种子/门槛效果前固定。来源范围与锚点范围不互相替代。使用describe_CM02_fine_state_structure_v1.py及audit_CM02_fine_state_structure_v1.py；PANCANCER_PATIENT_SCOPE_OUTPUT/AUDIT指向本预检和接受回执，PANCANCER_ANCHOR_LABELS是私有已核节点名配置，独立审核另设PANCANCER_FINE_STRUCTURE_OUTPUT。具体方法见总合同，本轮新描述与旧56宽计数问题分开。实际数值/科学终态待各自回执，不提前采用。


## 2026-10-21T00:51:28.714957+00:00 UTC /2026-10-21T08:51:28.714957+08:00 UTC+8 细状态全机会与表示差异解释完整接受

完整84细参考病例机会、876版本核对及所有节点/药物/阶段病例覆盖，独立Counter/Fraction/原计数表和排序分位数验收完成。后续完整168描述（来源与原命名锚点范围分开）及所有实际执行的2000随机流/留一/三分母由另一Counter/Pandas秩/SVD重放；原不支持位置、n恰达门槛时留一不足的情形均保留，不降门槛。所有42组合的完整范围/阶段剖面及反例已科学审阅，具体数量/效应只私有。

原宽/细计数版本差异已核实，不能直接把宽名称当功能身份。另在同一患者集合的四种输入/比例/协变量表示、两阶段共8点完成SVD对照，解释原相关对分母定义的依赖，不以此宣布数学错误、直接互作或换分母后的新临床发现。计数描述家族本轮收束；CM02仍优先，具体状态组合尚待独立空间与同终点临床支持。

下一优先核同细胞联合面板、分割/物理坐标和病例/区室层级；GSE300147同癌种Xenium为新结构入口，尚无同ICI阶段/终点资格，官方目录实际HTTP200仅整包tar，未下载大型数据。GSE288199继续作为同治疗场景候选，精确终点/分选/药物和空间资料待。原两版CM/正式注释保持，goal active未完成，拒绝后的公开push/新随访表无重试。


## 技术交付

[病例/版本接受](../inventories/analysis/2026-10-21/CM02_fine_patient_scope_v1/technical_scope_and_version_acceptance.json)、[描述/科学与表示接受](../inventories/analysis/2026-10-21/CM02_fine_patient_scope_v1/technical_fine_description_complete_acceptance.json)绑定实际执行与独立审核源码SHA、输入/结果终态、完整私有科学审阅SHA；公共记录不含科学数量/效应、患者行或学习权重。decompose_CM02_reference_denominators_v1.py以原完整描述/接受回执和私有focal配置运行全部8表示，只作同病例数学解释，原拒绝及未完成事项保留。
