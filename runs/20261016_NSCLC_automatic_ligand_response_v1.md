# NSCLC完整自动配体—接收程序家族：依赖、校准与正式任务

## 原生及完整输入事件

原生受体生产实际2026-10-16T10:56:36.342578+00:00 /18:56:36.342578 UTC+8完成，全部2010550708原坐标、压缩SHA/CRC及整数/排序门槛通过。独立原生审计10:57:12.041542 UTC /18:57:12.041542 UTC+8通过全单位边际/受体细胞集合/谱系汇总及预固定64原细胞不同解析器复读；不是全细胞独立原始重放。[原生run](20261016_NSCLC_prior_joint_RNA_and_response_preflight_v1.md)。

原生生产终态SHA `0d6e0b7608d52670dde93734016ea2bde01814912fa6fa514be97816a4e6fd50`；原生独立审计SHA `66d45f0e9e1dff7da72e08ab0e2aabbaf48012cec96e371a77b8d131ac021963`；完整接口资格SHA `8602cd810f01f4f409b12824fa63aa7e1bc71bb9ead54237f7d25fa5a481d927`。

完整输入预检10:58:11.131942 UTC /18:58:11.131942 UTC+8结束，45760逻辑槽/44528不重复组合和3339初始合格设计保持。终态SHA `7be8a884664fe9a8f13cfe769810b1de9b363d98ededc9df738229ee047b5a2f`，范围SHA `d6871a3075251c5bdffdaa6984ad5027aefcdd314e80dc11a7a4af195808f73d`。输入独立审计10:51:47.496708 UTC /18:51:47.496708 UTC+8已启动，源码SHA `ca4ffa6399f2449d6076e4b03dca6636c440e9619aded9cefc6b2b9384103990`；最近观察仍R，已核26000组合，完整终态尚待。不可估计原因为接口、基因剔除、覆盖、共线等，不能当阴性。

## 效果前统计校准和实际依赖启动

统计模块10:54:40.619059—10:54:43.690117 UTC /18:54:40.619059—18:54:43.690117 UTC+8完成真实校准。源码SHA `8f1579dc2aa44ac14bea2f80883799abe9029a7f2f7807b80ff287516537b70c`；校准终态SHA `a044c9e48af3dfd8834d1fd98aa11a4c562d6cb3f60c9d2af688472cf6861a96`。有限全部路径、批内准确停止、Statsmodels HC3、预设合成空模型/注入信号检查通过；不声称普遍wild精确有效或科学结果已验证。

正式任务10:56:59.223899 UTC /18:56:59.223899 UTC+8实际启动，源码SHA `cdc0dfc4e1c7c50196298e6fce2ee6732c94c17392d8b058989c40857b2f87fb`，当前等待完整输入独立审计；尚未读取效应。独立结果审计另已建立等待链，实际启动收据另补。两者均通过原PID/startticks和终态依赖，不能因短时观察超时重启。

命令和依赖：

```bash
python scripts/audit_NSCLC_automatic_ligand_response_inputs_v1.py
python scripts/NSCLC_automatic_ligand_response_statistics_v1.py
python scripts/test_NSCLC_automatic_ligand_response_v1.py
python scripts/audit_NSCLC_automatic_ligand_response_results_v1.py
```

方法及固定规则：[完整自动推断合同](../methods/NSCLC-automatic-ligand-response-inference-v1.md)。正式任务只在全部门槛通过后写`plan_before_any_ligand_response_effect.json`并执行；固定h100、999999上限、2000bootstrap、1900区间门槛、全部不重复组合BH及BY敏感性，不按结果改候选/控制或切点。4进程、BLAS1；环境沿用Python3.10.19/NumPy1.26.4/pandas2.3.3/SciPy1.11.4/Statsmodels0.14.6。

逻辑根目录：`${RUN_ROOT}/NSCLC_ligand_response_statistics_calibration_20261016_v1`、`NSCLC_automatic_ligand_response_preflight_20261016_v1`、`NSCLC_automatic_ligand_response_results_20261016_v1`。实际源码按启动SHA固定，公开commit发生在启动之后，不伪造运行前已提交。

## 当前边界与失败

当前正式效果/全部数值验收仍待。数据未复制本机或公共仓库；代码/方法/过程和SHA可公开，科学数值及患者行仅服务器/授权对话。单端临床、两端协调、先验、同细胞RNA、接收程序关联及直接机制分别登记，不相互替代。原9CM/172身份及旧阴性、失败历史保持，goal active。

输入审计遇到真实失秩设计的零奇异值除法warning，只记录无穷条件数及不可估计，不删控制或更改阈值；审计进程仍正常推进。PMC统计论文页面一度浏览器验证、原出版页直链500；已读取作者NUS讲义第12页及可读原发表文献页面，不声称绕过访问。无新表达/先验下载。

## 2026-10-16 10:59:46—11:01 UTC /18:59:46—19:01 UTC+8 门槛通过及真实计算

10:59:46.647331 UTC完整输入独立审计通过，SHA `c408b62e3ff94932d20fefa63395a51886dc198f9a259706a4dafe4c6fff5016`，全部3339输入设计逐文件SHA表为 `9a47d0d0204d4a00267b10f1cb0f598ec0986353ae4ea5d50757b169db7526b4`。原病例、基因逐项和、全部协变量/门槛均核，无效果筛选。

11:00:00.415191 UTC /19:00:00.415191 UTC+8真实写出效果前统计合同，SHA `df4eb54fb631766b8badecb5dc10055005a3014aa40177ee9384301c266d52c6`，随后4进程开始完整家族。11:01:03 UTC检查点已完成228/3339个输入合格模型；主进程原PID/startticks相符，独立结果审计等待原终态。此为执行进度，尚未审阅阶段性效果或发布科学结论。

独立结果审计实际10:59:19.079309 UTC /18:59:19.079309 UTC+8启动，源码SHA `580e758efb72b281d4648954d1ec2813162e35b3f15b9c9f69c771ae8c764063`。计划核全部点/敏感性/BH/BY/留一、首/中/首个达到抽样上限（无则末）记录的完整顺序随机流和bootstrap；不能预称审计已通过。
