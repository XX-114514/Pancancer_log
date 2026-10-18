# 优先免疫CM的TCGA与免疫治疗探索

## 用户要求与授权

用户要求结合已有数据、治疗及TCGA继续研究优先模块；数据不足时检索合适来源并申请下载。可复用公开数据和下载公开处理数据；受控申请若需机构/伦理/个人信息则先准备材料，不代填未知信息或发送联系邮件。延续每任务及时普通push，不公开患者行/矩阵/学习权重。

## 起点与范围

取证 2026-10-18T08:50:04.622193+00:00 /Melbourne 2026-10-18T19:50:04.622193+11:00。优先2026-10-14免疫图谱CM（v1，T_whole_TME_primary_tumor）的CM02 CTLA4–GZMB与CM04 GZMK–IL7R–CTSW。原CM/身份不改，10月16日HNSCC30槽阴性家族不重开，旧九CM/旧九月bulk结论不借用。

已读AGENTS、README/STATUS/ROADMAP/TODO、Cloud及RUNTIME_SYNC、最新daily、CM命名方法和原日志。服务器已有bulk_cohorts_20260904，含TCGA表达/临床及5个ICI来源，先盘点和核其原研究与单位/时间/终点。沿用临床转化skill及repo-implementation guardrails，读取旧reader后为本轮新问题写独立版本，不重拟历史模型。

## 分析原则与初始计划

1. 效果盲核实际文件、患者/时点、治疗、应答/时间事件、原表达尺度与研究重叠；发现/同队列和真正独立来源明确。
2. 先冻结原CM的单细胞端点及固定多基因功能读出，在已有全量native统计核表达依据；bulk仅为可测功能代理，不把bulk共表达当精确细群模块复现。
3. TCGA用于癌种内预后/背景，不能默认为ICI治疗队列；临床结局按TCGA-CDR建议，真实ICI来源单独分析。
4. 预先冻结全部主要槽、控制/可识别门槛与校正，所有失败/阴性保留；不按效果更换签名/分组/切点。
5. 条件充分时执行患者级探索与数值复核，交付图表和证据/缺口；公开数据不足则检索下载合适处理数据。真实推断执行尚待预检和合同。

## 初步观察

- 已定位TCGA RNAseq/TPM和临床，ICI包含IMvigor210、Hugo、Riaz、Jung NSCLC、Rose尿路数据；这是库存证明，不是本轮样本验收。旧bulk目录1.2GB，避免重复下载。
- 服务器缺rg，改用目标目录Python文件清单；原数据无变化。文献查询TCGA-CDR、UCSC Xena及GEO官方原队列记录；不依赖二手论文可能错误的GSE/治疗配对。

## 效果前合同冻结

固定两轴多基因代理，7个公开生物学/控制panel，不从bulk结局学习基因。全33癌种×2轴×PFI/OS×基础/分期敏感性共264 TCGA槽；ICI5队列×2轴应答10、OS/PFS20，独有增量10另列完整家族，缺测/失秩保留P1内校正并展示NA。详情frozen_plan_before_effects.json。主要Score为来源/癌种内gene percentile平均，再取两臂min，真实CM投影与该代理明确区别。原HNSCC30家族不重开。

实查Hugo有27pre/1on，Riaz51pre/58on；新readers将排除on、保留SD为客观应答NR（不沿旧reader只保留PD的定义），未知删失保留缺失。IMvigor347均Atezolizumab，RESPONDER49缺失不填。Rose有io.response但时点需证明；T/N/M字段不当AJCC整体分期。

实际data4组配额按文件系统精确匹配，余量约501.9GB（十进制），本轮仅读取现有矩阵，采用100MB输出上限+1GB额外余量、500GB留空；不复写740MB源表达或写大中间矩阵。首次元数据shell包装引号SyntaxError在执行前，改为本地stdin脚本后完整metadata审计成功，无科学数据改变。
