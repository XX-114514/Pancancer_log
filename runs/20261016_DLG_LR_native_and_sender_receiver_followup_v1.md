# 原CM通讯候选：原生RNA与发送/接收分解

## 输入和当前终态

2026年10月CM（v20261010）的原172身份、14个已拟合上下文、冻结DIALOGUE基因贡献和完整每样本LIANA保持。沿用[原LR/TF run](20261016_DLG_LR_receiver_function_v1.md)，本run负责原生复合物与条件RNA描述。科学值只在授权对话及服务器私有目录审阅，不在本记录复制患者、效应表或学习签名。

当前已完成独立全点复算、全部原生面板读回及完整条件描述验收。不能据此确认直接互作、正式状态、跨癌种独立验证或临床价值。原P/BH和候选家族未改变。

## 实际时间线（UTC /UTC+8）

- 06:11:52 /14:11:52：原全LR生产终态，当前任务于08:54后读取该既有结果，没有重复生产。
- 09:08:52—09:09:37 /17:08:52—17:09:37：独立所有可估计点效应、敏感性、内部留出及BH复算完成。仅少数常量敏感性向量产生预期NaN警告，原主要结果不改；未逐项重放全部原随机化分布。
- 09:11:59 /17:11:59：原生面板v1启动，编号同时为index及column导致merge歧义，原日志和输出保留。
- 09:18:12—09:18:38 /17:18:12—17:18:38：v2只重置右表index解决连接；随后发现一个来源的固定内皮标志不可测，汇总全缺失列求和会误成零。
- 09:22:04 /17:22:04：v3保留缺测集合为NA并显式列可测标志，同时区分完整上下文与原推断观测分母；原始计数/候选/统计不变。
- 09:23:28—09:24:17 /17:23:28—17:24:17：独立全部原生选定细胞/基因/单元和所有病例汇总验收通过，未宣称表面蛋白。
- 09:25:07—09:25:09 /17:25:07—17:25:09：原14条各两个条件RNA描述完成，无新P。
- 09:26:31 /17:26:31：全部28个描述及1,088项留一块独立OLS通过；本轮取证晚于上述实际事件，未虚构早期公开时间。

## 版本和证据

代码SHA固定于各服务器launch_receipt，源码为本提交`scripts/`同名文件。环境继承scanpy Python/NumPy/pandas/SciPy/Statsmodels、只读canonical Zarr；计算线程固定1。新建私有目录，不覆盖或删除旧版本，无新表达下载。

|证据|逻辑位置与SHA|
|---|---|
|完整点/BH复算|`${PREFLIGHT_ROOT}/DLG_LR_point_binding_review_20261016_v1/private_terminal.json`；原效应表SHA `e01cf5519c54af81e7669a271148534cf768849a1fe3120de0ef6760ded2d0b8`|
|原生面板v3|`${PREFLIGHT_ROOT}/DLG_LR_native_receptor_panels_20261016_v3/private_terminal.json`；终态SHA `82ea416611f80ca4a58cc95938b23a7b37895f77d5b9e3ab419c19245684a740`|
|原生独立验收|`${PREFLIGHT_ROOT}/DLG_LR_native_receptor_panel_audit_20261016_v1/private_terminal.json`；SHA `20626974010968afb6c49ab3e4bb43198985fbc2e6051c0a64b3ce2b7c4e478f`|
|完整条件描述|`${PREFLIGHT_ROOT}/DLG_LR_sender_receiver_decomposition_20261016_v1/all28_private_sender_receiver_RNA_diagnostics.tsv`；SHA `d322628a23d51cf5f96fbb596e56d8b96bd0ad62a7dfe33385f6c67e1632c1bc`|
|独立OLS及全部留一|同目录`private_independent_audit.json`；状态PASS_ALL28_STATS_MODEL_OLS_POINTS_ALL_INTERNAL_HOLDOUT_ALL_LEAVE_ONE_BLOCK_OUT_DESCRIPTIONS|

完整步骤见[方法](../methods/DLG-LR-native-receptor-and-conditional-RNA-v1.md)。直接运行各Python源码须提供`PANCANCER_PREFLIGHT_ROOT`及需用时的`PANCANCER_PROJECT_ROOT`，拒绝复用输出目录；正式生产通过已有隔离detached launcher，PID/startticks和私有stdout留原控制目录。原生v1/v2不能用于后续已修正汇总。

## 下一动作

审阅已有TF任务的真实终态和独立验证，不重复派发。结合全部条件描述区分发送端线索与接收端RNA主导的评分关系。原程序功能和空间测量仍需验证；不因为统计关联或RNA共现命名“耐药”“免疫抑制”。[NSCLC原生输入](20261016_NSCLC_frozen_programme_measurement_v1.md)用于下一来源/临床入口，先固定谱系与冻结变换的支持范围，不能按当前病例疗效挑程序或删控制。

后续取证更新：09:51既定TF生产和其全部官方活动/点效应/内部留出/BH独立审计已完成，当前任务只读其终态及本原生家族对应记录，科学判定仍按[TF主run](20261016_DLG_LR_receiver_function_v1.md)。本run前述09:32证据卡是当时快照，不能由“功能待运行”的历史文字覆盖之后的实际终态。各程序和接口的支持/不支持保留原编号，不以其他程序的临床结果背书。
