# 免疫治疗专用单细胞验证对象与下游分析

## 用户授权与目标

用户要求：对新完成的3180973细胞预处理Zarr继续下游；在Project_v5查找合格免疫治疗数据，共同组成专用对象，按当前Project_v3方法降维/去批次/聚类、证据注释和逐样本恶性预测，再识别CM。新增输出独立版本，最终原项目卷，备用卷仅临时；沿用逐关键节点普通push授权。未授权删除或修改旧权威对象，新的矩阵公共下载仍不恢复。

## 开始状态与已读证据

已读AGENTS、README、STATUS、ROADMAP、TODO、RUNTIME_SYNC、当前表达注释下游run、V7参考/恶性方法及CM定义合同。服务器读大队列、证据谱系、逐样本inferCNV、R CoVarNet、单样本QC和临床解释skill及其相关执行/储存参考；PanCancer及annotation_rebuild规则、CURRENT_STATE/CURRENT_RUNS。

已有新数据技术合集9来源792文库3180973细胞/49985特征，全行PASS；来源正式AC待。旧免疫四队列的sc/snRNA历史约1447081，与Visium分开。真实矩阵盘点另显示GSE305078297261细胞；ICB作者汇编355941细胞/8原始研究，raw retained counts与native nCount归一化审计PASS，但存在基线重叠与端点语义问题。

## 有界实施顺序

1. 读取研究/患者/样本/时点/治疗/材料、已有QC和计数证据；按原始研究、文库/barcode排重，核发现使用历史；登记独立验证、内部复现、机制参照和隔离角色。富集免疫、完整TME、CAR-T产品、sc/snRNA分背景。
2. 复用已审QC，不全局二次过滤；必要新增输入逐原文库QC/双细胞及污染可辨识性评估。冻结统一counts/来源掩码/稳定细胞ID和资源峰值；最终权威指针需实际验收，不据目录或新命名直接晋升。
3. 原生RNA提取HVG/PCA，保留治疗/癌种/材料，校正技术批次；计算整合前后混合及生物保留指标，全细胞注释与各大类重聚类。遵循多marker/作者/模型证据与边界标签；显示名单marker加数字，node_id单列。
4. 从原始文库QC counts、同文库候选/参考运行已有malignant_inference核心；不足明确skip/uncertain，不把未知当normal；细胞ID+样本键精确回填。肿瘤身份与CNV轴分开。
5. CoVarNet候选发现与既有9CM固定投影分开，rank/重复拟合/来源敏感性；临床验证仅用预先冻结独立患者队列，不以重新学W称旧模块恢复。

## 实际执行、失败及未验证

截至初始记录只读核查；未声称新整合/注释/CNV/CM完成。一次目录检查SSH连接超时，缩小范围重试读取成功；未重启任何科学任务。网络论文/GEO读取遇验证码或Nature跳转限制，优先继续核服务器已缓存官方材料。scanpy/harmony/BBKNN/celltypist/infercnvpy/leiden依赖存在，scanpy环境未发现scvi；内存约1TB，总可用约800GB。有效组配额余量547249078272bytes的当前只读快照，需要按下一阶段实测峰值保留500GB。

## 验证与同步

每阶段记录输入/代码/参数/真实终态/独立检查；不新增无关根审计或重复测试。公共仓库只发布代码、计数、状态和脱敏方法，不发矩阵、患者行或内部绝对路径。关键节点普通commit/push后独立HTTPS核验。

建议初始提交：plan: establish immunotherapy atlas scope and downstream milestones。
