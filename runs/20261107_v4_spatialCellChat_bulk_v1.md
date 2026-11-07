# 当前V4成像空间通讯、TCGA免疫与六ICI

用户要求基于重点CM的空间复现，推进单细胞成像型空间转录组CellChat，并同步做TCGA免疫和ICI队列分析。对象仍为2026-11-06 V4主v2（14因子、45非上皮节点、164样本、13来源），原注释SHA `2a4bc82cb495d0fc2ffb40ecda6b1cbb53802eb2bf87efd2392746371653244c`。原空间与25CMT队列不变。

逻辑运行根`${PROJECT_ROOT}/runs/20261107_v4_spatialCellChat_bulk_v1`；[方法](../methods/V4-spatial-CellChat-TCGA-ICI-20261107.md) · [实际状态与资格](../inventories/analysis/2026-11-07/v4_spatialCellChat_bulk_v1/README.md) · [代码](../scripts/v4_spatialCellChat_bulk_20261107/)。固定原参考训练/留出与原生面板，不安装共享包。实际环境：原scanpy CPU Python3.10、R4.3.3、Seurat5.3.1、CellChat2.2.0；各科学产物另绑定实际版本/代码SHA。

## 真实时间（UTC；UTC+8另加8小时）

- 07:38:37：原Visium8/53完整生产/评分/独立审核通过、另2张拟合；原CellCharter53/53接受。07:38:42原CPDB已有CM01–10 High审核通过，CM11 High在跑，完整Low对比未完成。
- 07:39–07:44：只读定位现有TCGA、六ICI输入机会、5 Xenium/8 CosMx归档，核对父工作区规则、坐标单位、面板及既有处理链；不读取效果来选择模块。
- 07:48:26：作者SU2C公开三成员CRC/SHA接受，全TPM44,557,830字节/152 RNA列，本机保留。原34基因子集不适用于本轮。
- 07:49:48–07:51:15：同560,871细胞的164样本/14,051共同基因pseudo-bulk准备exit0，13个原生空间机会完整登记。
- 07:51:15–07:51:51：全部14CMT签名资格推导exit0，Seurat5.3.1；5个非空因子具完整10基因，CM06不足、CM14仅2High等保持不可评估。07:51:51–07:51:56数据库导出exit0。
- 初次启动SSH客户端超时未收到ACK；后续三份终态与同注册进程证据确认科学准备实际完成，未重复启动。
- 07:59:33：有限两CPU分支已进入首个Xenium身份准备；bulk等待读入精确SU2C本轮投影。原GPU/CPDB不重启。
- SU2C完整上传慢、一次大stdin传输超时，保留部分文件不用。精确88基因派生保留全部152 RNA，源/派生/冻结签名SHA绑定；13小块逐块校验，仅一个跳板超时块重试，08:06:26组装后完整SHA接受。
- 08:03:11：首个乳腺Xenium准备exit0，精确CM03身份仍不可测；原int8组合码导致邻域背景阶段exit1。追加int64修正保留原源码/失败，图/种子/阈值不改。首个CellChat返回不可评估，不把exit0称成功拟合。
- 08:06:32：第二乳腺Xenium准备exit0；08:06:38完整14父类背景统计写出待独立审阅；08:06:54其空间CellChat审阅同样无可评估拟合，原因另核。两份身份训练出现fine分类器迭代上限警告，不能将概率/留出重放当独立身份真值；后续收敛与面板细分资格需审阅。
- 08:06–08:07：bulk正式读取六ICI及TCGA。已精确计数来源内患者为IMvigor347、Riaz51、Hugo26、Jung27、Rose88、SU2C Keep136；模型门槛另适用，小队列不强拟合。08:07:29已进入TCGA癌种逐个拟合，未收取完整终态/效应，不能提前报告临床发现。

## 配置与验收

空间四组合100/250μm×10/20μm、100表达置换、小面板不做DE、蛋白且完整组件实测；几何缩放规则及额外999次FOV内标签置换分别登记。原生细胞/坐标/矩阵与FOV明细全部私有。5 Xenium均缺IGHD及XCL1，不能用大类B/内皮/NK冒充精确CM03；CosMx逐归档另核。

临床2,100槽与免疫4,620槽的完整BH/BY、病例得分和拟合设计绑定、直接Cox/logistic数值及免疫统计另审。临床程序为癌种/队列内完整10基因百分位均值，非原H或作者AddModuleScore原样重放；无阈值优化。未知事件/阶段/应答保留缺失，原空网络和缺测槽保留。

只公开已审阅的配置、代码、来源/标签汇总和脱敏收据；不公开患者/细胞/FOV坐标行、矩阵、设计或模型。生产、数值审计、身份资格、科学结论与公开同步分别记录；失败历史不删除，不自动推送。

## 较晚完整bulk验收与真实空间变体

08:06:34–08:08:31 UTC临床2100槽和免疫4620槽生产exit0；395设计可拟合。08:09:35 UTC／16:09:35 UTC+8全部395设计另算法/直接Efron与logistic重拟合和完整BH/BY通过，最大差2.81e-6。首轮测量审核因CSV患者索引类型推断失败；仅按原字面字符串读取，不做模糊匹配。随后统计读回发现默认浮点解析改变并列秩；round_trip读数在08:20:17完整验收236136得分条目、395设计、1305原RNA统计与签名FC，最大P/统计误差3.33e-16、FC6.22e-15。生产分数和模型未改，旧失败保留；SSH ACK超时后压缩部署只读回已存在PASS，没有重复运行。

另查原始免疫亚型的六个完整描述，而非裸C1–C6，发现初版类别过滤未覆盖原字符串。按六个原标签精确映射并用独立正则复核；08:28:21 UTC／16:28:21 UTC+8当前完整4620免疫家族、20可检验类别槽和1305 RNA槽通过，H/P最大差2.28e-13。原RNA统计/P与临床模型不变，所有BH/BY依完整家族重算，旧v1保留，v2为当前免疫结果。

当前结果：TCGA350可估计/5 BH/2 BY；ICI应答20可估计/0 BH，生存25/0 BH；免疫1325可估计/814 BH/688 BY。CM03各可估计ICI终点未获BH支持；LGG OS基础模型HR1.590、BH0.00679/BY0.0550是单癌种程序关联，缺分级/IDH等背景，不称泛癌或治疗标志物。全14/33癌种/六ICI三图已看过，修正重复行名与colorbar拥挤，底层统计不改。

空间预检0×0问题准确定位于安装版computeRegionDistance读meta$group、初版却只传labels。v2仅补相同分组字段，原CellChat_v1空矩阵/不可评估收据保持；不是数据没有邻近群。首张int64背景于08:15:15恢复exit0，其他成功准备/背景不重跑。08:29:38 UTC首乳腺与卵巢各2/4个真实空间变体已完成，两个R CPU时间约860/914秒持续增长；全部4变体独立审计及后续归档仍待。5个Xenium全部准备，fine均迭代到上限，三IO大类也到上限，审阅版精确CM资格因此不晋级。

CosMx在进入原归档前用小型tar流/Pandas互操作检查确认seekable接口问题，以相同成员字节的BytesIO缓冲修正，未改变CSV数据；cell_ID0按原定非细胞背景先排再做字面键核验。缓存只在同运行、完全相同面板/参考/顺序时复用分类器，每单元仍重算留出诊断；原版和新SHA完整保留。

08:36:27 UTC／16:36:27 UTC+8首乳腺4/4空间v2变体独立通过；08:40:23 UTC／16:40:23 UTC+8卵巢4/4数学审核通过，后者大类未收敛，仍不晋级身份/生物解释。首个CosMx归档实际表头含IGHD/XCL1/PECAM1/VWF，fov/cell_ID与global_px坐标字段齐全；不是全部8归档已通过。后续相同L2/C1目标预算5000、未收敛缓存warm-start续优化，原预测不覆盖；原型1000次和新SHA均保留。CosMx分视野生产/身份及全空间结论仍待。

## 脱敏可复现命令形状

工作区和输入目录先按原冻结配置准备；以下宏由既有用户环境绑定，不含服务器地址/账户。Python阶段通过runpy注入相同PROJECT_ROOT，不直接裸运行缺少该上下文的脚本。

```text
${SCANPY_PYTHON} -B ${RUN_ROOT}/launch/finite_V4_extension_controller.py ${RUN_ROOT} ${SCRNA_R_RSCRIPT}
${SCRNA_R_RSCRIPT} --vanilla ${RUN_ROOT}/scripts/derive_V4_CMT_signatures.R ${RUN_ROOT}
${SCRNA_R_RSCRIPT} --vanilla ${RUN_ROOT}/scripts/run_spatial_CellChat_V4_v2.R ${RUN_ROOT} ${UNIT_CODE}
${SCRNA_R_RSCRIPT} --vanilla ${RUN_ROOT}/scripts/audit_spatial_CellChat_V4_v2.R ${RUN_ROOT} ${UNIT_CODE}
${SCANPY_PYTHON} -B ${RUN_ROOT}/scripts/audit_V4_bulk_models.py --run ${RUN_ROOT}/bulk/analysis_v1 --output ${RUN_ROOT}/audit/bulk_models --expected-slots 2100
${SCANPY_PYTHON} -B ${RUN_ROOT}/launch/finite_V4_interface_recovery.py ${RUN_ROOT} ${SCRNA_R_RSCRIPT}
V4_REVIEW_ROOT=${REVIEW_ROOT} ${PLOT_PYTHON} render_V4_bulk_aggregates.py
```

原v1准备/生产和失败日志按各自时间保留；上述命令不得在已有产物目录盲重跑，控制器只创建新阶段或复用已接受输入。实际generic脚本及部署SHA见公开源码清单；连接与私有输入映射仅用户研究环境保存。

08:53:57 UTC取得首乳腺大类背景聚合表：133517模型细胞、11类、15完整蛋白互作，主参数273名义候选/原软件BH265。08:57:18 UTC／16:57:18 UTC+8全四变体非零Monte Carlo下限诊断通过：100次重采样原P0经(100P+1)/101和全1815假设BH后各0候选，最小q0.0593–0.0678。原数组/P未改、无细胞替换/模型重跑，不将候选或软件q0称稳健通讯。通用诊断R、聚合表和SHA v3保留；此诊断不是空间零模型或完整tie校准，精确CM03仍未验证。

09:00:11 UTC／17:00:11 UTC+8最终取证：I001/I002/I005共3个原生Xenium单元各4/4距离变体及独立数组审核通过；I006首CosMx为10/112个视野×参数变体完成，两有限控制器均有活跃R子进程。3/13计算完成不表示精确CM03身份或通讯验证，112不是患者数。
