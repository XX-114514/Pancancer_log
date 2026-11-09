# ccRCC补充探索：独立身份、蛋白与局部血管邻域

本轮针对**2026-11-06 V4 CM03**的血管—淋巴细胞关系做单供者ccRCC补充。它不是IGHD B—XCL1 NK精确亚型复现，不替代多患者乳腺癌或新TME发现。

## 先纠正两个测量问题

旧EC规则把ACKR1与PECAM1/VWF/EMCN共同计数，再比较ACKR1检出组，存在身份选择依赖。原规则及结果保留，原代码/输入/距离/零模型先在服务器原样复算；新主规则只用**PECAM1、VWF、EMCN中至少2个实测阳性**定义EC，至少2个淋巴管marker阳性者另列，ACKR1不进入身份门槛或RNA聚类特征。主视图排除规则间谱系冲突，但保留原始对象与复查池；这不是新增自动doublet判定。

蛋白矩阵`protein_metric=MEAN_PER_CELL_STAIN`，文件内`protein_scaling_factor=10`。读入后除以存储缩放因子，得到任意单位的平均荧光强度；RNA与蛋白按feature_type独立提取，即使名称相同也不混算。旧报告94、28、588等是**存储缩放值**，除10后的同抗体值才是本轮表格口径，不能当作RNA UMI，也不能跨抗体作绝对强弱比较。[10x官方蛋白输出说明](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/protein-data)

## 输入、恢复与空间坐标

只需要官方56.7GB归档中的六个原成员，加58.6MB的`cell_boundaries.parquet`，共120,420,096个未压缩字节。六个核心成员必须匹配此前完整下载的独立SHA、ZIP CRC与大小；新增轮廓成员校验ZIP库存CRC及大小，并记录新SHA，不冒称已有独立SHA或整包MD5验收。块缓存的接收不等于完整成员验收。有限重试、持久块SHA缓存、每块写锁、失败保留及原子发布沿用既有恢复工具；四连接轮廓预取与三连接核心预取共享缓存，各有期限，不强制解锁、不取消旧作业、不改持久代理。

矩阵barcode与cells.parquet逐一核对；原生x/y_centroid保持微米，不乘图像pixel_size。轮廓vertex_x/y也按原生微米核对到对应EC质心。显示图像原点在上方。[官方输出及坐标说明](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/outputs/xoa-output-understanding-outputs)

## 身份审阅

沿用RNA≥50 UMI、≥20基因的QC。B至少2/3 CD19/MS4A1/CD79A，T至少2/3 CD3D/CD3E/TRAC；NK需KLRD1与NKG7或GNLY，同时无上述T marker并排除B/EC。血管周规则沿用PDGFRB/ACTA2/MYH11。marker图同时标基因和预期身份，缺测用×，不把缺测记为阴性。

将QC合格且有免疫程序、血管/基质程序的细胞分别作RNA-only诊断聚类，包括规则冲突细胞。只用实测RNA，排除ACKR1及所有蛋白、空间坐标、厂商空间语义标签；CP10K/log1p、30PC/30邻居、Leiden0.5，种子0/1/2，报告ARI和全部簇的marker、规则组成与同抗体蛋白分布。聚类编号不是新冻结亚型；先审阅混合/弱程序，不按空间检验结果改注释。

## 完整空间探索家族

先固定全部QC非EC/非淋巴管接收细胞，再按最近EC互斥归属，比较ACKR1≥1 RNA UMI与未检出EC的邻域；未检出不能命名为阴性血管亚型。

- 几何：最近EC质心；最近**分割EC细胞区域**作为敏感性。后者不是完整血管壁重建。轮廓距离非唯一或为0的接收细胞均独立于ACKR1状态排除，并报告数量。
- 半径：25/50/100/200μm；接收者B/T/NK；分母为全部QC非EC或仅B/T/NK。
- 局部零模型：250或500μm网格×EC RNA总量floor(log2)分层，在层内置换ACKR1检出标签，固定坐标、身份及每层目标EC数量；999次，共享置换用于全端点。
- 统计量：加0.5的邻域组成log-odds。单侧增强P=(1+null≥observed)/1000，含并列、不报P=0；2几何×4半径×3身份×2分母×2局部尺度=**完整96槽**，统一BH/BY。预留无法评估项并以P=1参与家族，不据最小P选展示。
- 同时报原始比例、百分点差、条件零模型期望与超出期望的量；原始效应与条件P的解释不能混淆。

这份设计在读取本轮新效应前保存，但此前ACKR1描述性结果已看过，因此仍是探索性随访，不能称独立确认。位置与粗RNA深度分层不等于已控制所有组织学/TLS或连续检测偏倚。

## 验收与边界

原试跑重放检查690,322输入、全部QC/身份数量、所有距离及旧置换值；新分析独立检查RNA/蛋白分离、缩放、完整96端点、P/FDR、距离与聚类特征/ARI。合成反例检查ACKR1泄漏、EC/T冲突、轮廓并列/内部点与多重校正。

单供者/单切片，细胞、血管、网格均不当患者数。477 RNA/27蛋白缺IGHD/XCL1/NCR1/CCL21，不能验证精确亚型或CCL21信号；尚无病理TLS真值。厂商带血管/空间关系语义的标签不作身份真值；蛋白交叉信号需要结合分割审阅，不能直接断言doublet。此次不新增通讯、TCGA或ICI拟合。

[本轮实际结果与状态](../inventories/analysis/2026-11-09/ccrcc_followup_v1/README.md) · [运行记录](../runs/20261109_ccRCC_followup_v1.md)。大型源、缓存、逐细胞对象留科研服务器。
