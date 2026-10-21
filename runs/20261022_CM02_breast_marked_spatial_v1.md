# CM02乳腺来源：预固定连续标记的空间配对

## 前态与目标

前轮3b20e86已完成固定Danenberg v2.1.0全成员原生资格。真实配对与预留子集不足以直接拟合复杂临床增量，但可检验有限的状态/空间问题。保留2026-10-14免疫图谱CM02 v1（T_whole_TME_primary_tumor）优先级，不将当前蛋白代理称原CM02投影、纯Treg、受体通讯或临床验证。

## 在新效应之前固定

- 只用明确主原发肿瘤、非hotAggregate、固定marker/XY有限和正面积的原细胞；临床键连接且作者混合Reg/Tex和CD8各至少20细胞的病例，无事件/ER方向选择。纯CD4比较只用其中也达到20的病例。
- 以每ROI/父类内平均秩定义FOXP3+CD4、ICOS+OX40+GITR两个连续轴，Ki67作一般增殖参照；CD8端为Ki67、PD1。等权预定义、无二分阈值或结果选择marker，混合身份不改。
- 固定H3、DNA1、DNA2、log1p面积测量背景；非Ki67轴另外调整本细胞Ki67，Ki67从不回归自身。完整设计需满秩、condition<=1e8且残差df>=8；ROI两父类各>=10、病例各>=20，不自动删协变量。
- 主图为ROI内25微米中心距离边（作者像素标尺），边>=20、连接细胞双方各>=10。按双方度数平衡边权，计算标记残差在异类边两端的加权Pearson相关。它是中心邻近状态关联，不是膜接触或LR通信。
- 所有QC细胞centroid包围盒固定2x2宏观网格；同ROI/父类/网格内联合置换整个标记残差向量，保持父类数量、图、位置、同cell多marker关系。双方各至少10可交换细胞，2000置换，病例等权汇总；point为观察相关减null均值。病例bootstrap2000给区间，seed20261022。
- 两父类对×两个source轴×两个CD8轴8机会，加四个同病例Reg-minus-CD4差值，共固定12主槽BH；n<20或不可估计保留缺失/原因，不能换检验。Ki67 source四槽仅诊断。作者预留标志/ER只报告覆盖方向，不做临床模型或重新选子集。
- radius50/grid2、radius25/grid4及radius25/grid2+作者interface/perivascular置换块为三项敏感性，仅效应/病例bootstrap CI，不新增P或挑最好设置。完整结果后再规划外部来源，而非把同队列文献再分析称独立发现。

## 执行与隐私

单CSV流读重新核固定SHA 9f693514acc02bbb6b2dabd14c15c806348c5b882659aea071a83ed16e2d5df0、CRC、前后目录/表头。原先冻结传输不改；最多768MiB/2h/请求45s/8GiB RAM，不写任何个体/源payload/矩阵。只保存来源级效应、覆盖、不可估计原因、component一致性和技术收据于公共Git外。原科学家族、W/S/H/权威标签不变。

当前协议先冻结，计算核/测试与真实运行均尚未开始。GSE288199 Table S1待用户补，三个受拒动作不重试、云端push待同步。

观察 2026-10-21T18:16:08.896844+00:00 UTC / 2026-10-22T02:16:08.896844+08:00 UTC+8。


## 实现与合成核查先冻结

计算核与输入worker实际实现，八项合成检查通过：秩ties/常量、背景混杂移除、共线与df不删协变量、已知正负空间关系与独立相关复算、网格联合置换及vectorized-null逐draw重放、未知组织标志独立块、来源病例不足保留不可估计、BH固定含缺失槽分母。源实际正/负边和null首draw重放及真实不可估计df检查内置，须运行后核是否触发。核信号、面积及非Ki67轴的Ki67背景均在ROI/父类内秩标准化；source多marker轴等权，无学习权重。

输入源重读及任何新效应此时尚未开始。科学启动前保存实现SHA与本地提交；外部设计/实现收据独立。运行环境为既有本机Python3.13/NumPy/pandas，OPENBLAS_NUM_THREADS=1，不安装包。命令（逻辑路径）为：

```bash
OPENBLAS_NUM_THREADS=1 PANCANCER_EXISTING_CLASH_PORT="${EXISTING_LOCAL_PROXY_PORT}" python scripts/CM02_Danenberg_marked_spatial_v1.py --qualification "${PREVIOUS_RUN_ROOT}/private_native_capture_scope_v2.json" --frozen "${RUN_ROOT}/private_frozen_marked_spatial_implementation_v1.json" --output "${RUN_ROOT}/private_marked_spatial_results_v1.json"
```

主12槽必须完整报告；优先候选解释需同时考虑BH、病例bootstrap区间、敏感性和同病例CD4差值，不以任一个P或单个marker状态名证明特异机制。残差置换的条件可交换性是模型假设，不是随机干预或因果证明。所有非主配置/子组/一般Ki67参照禁止新增P。
观察 2026-10-21T18:20:31.897074+00:00 UTC / 2026-10-22T02:20:31.897074+08:00 UTC+8。


### 在真实输入前收紧配对ROI合同

同病例Reg-minus-CD4差值额外固定为同一eligible ROI交集内先作差再平均，防止同病例的不同切片位置被当成状态特异差异。新增一个不相交ROI不得配对的合成检查，九项通过；实现冻结v2收据替代v1用于执行，v1保留历史。此前未启动任何真实输入或新效果。命令的--frozen使用private_frozen_marked_spatial_implementation_v2.json，其余协议不变。
观察 2026-10-21T18:21:20.639685+00:00 UTC / 2026-10-22T02:21:20.639685+08:00 UTC+8。
