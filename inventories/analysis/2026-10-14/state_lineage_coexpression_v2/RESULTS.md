# 原生同细胞RNA诊断：稀疏信号、状态候选及完整CM输入缺口

## 执行及数值验收

v1真实02:14:32—02:14:34 UTC exit1，categorical默认observed=False产生未出现的空组。v2只改observed=True，新目录/同240533cells和门槛，02:17:12—02:17:42 UTC真实exit0，两个句柄退出；02:20:03独立审计完成。

全部5source×2major×11panel的110记录、55gene coverage、212实际source-node/foreign-panel、990module尾部记录保留。7937589个原生cell/panel计数/检测/fraction值逐gene重新累加复算，实际node key用unique字符串独立检查、全部990尾部Python稳定排序复算通过。v1首source私有cell/panel量与v2完全一致；故本次恢复没有改变表达、细胞、门槛或身份。[审计](independent_audit.json)。

## 实际观察：不能由稳定因子推大范围身份错误

全部55panel实测基因完整。下表是“至少3个panel基因检测，且panel占完整RNA至少1%”口径，**不是生物/双细胞判定门槛**。

|观察|结果|允许解释|
|---|---|---|
|最高source-major强foreign-panel比例|GSE148071 Mac stromal44/6724，0.654%|强的同cell跨谱系RNA混合在此口径下罕见，不支持大范围免疫身份错误 |
|同一源/大类达到5%RNA的stromal事件|4/6724，约0.0595%|少数事件，来源未定，不直接删cell或称doublet |
|T epithelial事件最高source比例|GSE1480716/2333，0.257%|其他来源更低；仍不能排除低水平环境RNA/混合 |
|GSE182109 Mac neural≥3marker|326/78565，约0.415%|其中占RNA≥1%仅1个，存在弱/稀疏RNA而非大量neural身份 |
|T own-core至少2marker检测|各source约61.8%—92.4%|受panel及稀疏度影响，低检测不等于不是T；没有正式重新鉴定 |

多数foreign panel的95/99分位RNA很低；同时表达myeloid和stromal的44cells仍无法由RNA判定双细胞、ambient或摄取RNA。**这是面板口径下的诊断，不是全量细胞身份审阅。**

## 全90模块尾部提供的具体线索

固定各模块normalized usage最高300/6000（5%）只是描述性尾部，不是state阈值。[全990记录](all990_discovery_module_tail_diagnostics.tsv)。

- T K15_05 FOXP3/IL2RA/CTLA4相关programme尾部含30/34个满足Treg_related panel较强RNA口径的发现细胞；支持programme内容，但type与state仍需在原编号内部区分。
- T K15_10 G2M尾部含31/31个该panel较强RNA事件；支持增殖相关候选，不强制每source固定5%cycling。
- Mac K15_05 IFN相关尾部含61/167个该panel较强RNA事件；不证明唯一IFNG刺激来源或胞外通讯。
- T K15_11 stromal尾部只有3/300满足相同strong-panel口径，panel median RNA约0.0673%；Mac K15_11相应15/300，median约0.0609%。稳定胶原因子可能刻画弱或少数RNA，不直接作为“胶原合成免疫细胞状态”。
- 同理neural/keratin/mural候选保留作QC/来源/低水平RNA解释，尚无证据判定其技术或生物来源。cNMF稳定性与cell-level多基因支持是两件不同的事。

## 原9CM的实际输入范围

[逐member映射](all9CM_state_input_member_map_v1.tsv)、[全部9CM范围](all9CM_state_input_scope_summary_v1.tsv)直接匹配冻结graph core的原node标签，不合并编号/相似marker。当前2major native输入仅占各CM约0—26.98%的**原参考core权重**；CM06为0。此权重不是新样本活动、细胞比例或生物重要性。全部state_scored_nodes仍0。

因此后续必须补原单核、树突、Myeloid、T/NK/未决lymphoid、B、Mast/Neutrophil及fibro/endo/mural/epithelial等，不以首批T/M完成9CM目标。现有240533输入及来源限定programme仍有价值，原9CM/172和权威库保持。

## 下一步

固定方法候选字典，明确programme描述、within-method支持、来源/复杂度/identity及QC-sensitive标签；不能把跨K重复57候选称57状态。用原生完整gene轴UCell评分，同时检查每cell/programme实际检测/UMI支持，连续多程序而非强制type替换。扩充缺失CM谱系，再进行身份内状态/样本级LR/接收响应，完成状态后DIALOGUE和新空间验证。原阴性/来源限定阳性保留。

## 网络及发布恢复

SCP和JSON/gzip载荷观察超时；基本TCP/SSH printf仍可用，ssh alias当前无ProxyCommand，不能据此归因Clash或重启科学任务。IPQoS=none及IPv4单file观察仍失败；两个较小聚合先经cat取回，最后将两个缺file压缩为9014字节tar一次取回并全SHA通过。父诊断/独立审计已实际结束，传输失败不代表科学失败。所有私有cell/患者/score量仍在server，公开仅aggregate/genes/无数据代码。
