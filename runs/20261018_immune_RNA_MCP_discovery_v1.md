# 20261018 免疫RNA MCP发现 v1

## 范围、执行与当前状态

[方法](../methods/immune-RNA-MCP-discovery-design-20261018.md)。真实执行Python为scanpy环境、R为scrna_r及已验证隔离库；源码SHA绑定而非假称未提交代码已在某Git commit运行。标准入口：

```bash
python scripts/freeze_immune_RNA_MCP_discovery_scope_v1.py
python scripts/prepare_immune_RNA_MCP_native_expression_inputs_v1.py
python scripts/audit_immune_RNA_MCP_native_expression_inputs_v1.py
python scripts/run_immune_RNA_MCP_official_family_v1.py
```

运行时设置PANCANCER_PREFLIGHT_ROOT与PANCANCER_DIALOGUE_RSCRIPT为实际环境路径；根目录不公开。对应原生输入、病例行和基因图均留服务器。

## 时间线（事件UTC；括号为UTC+8）

- 11:47:51（19:47:51）：官方TNBC SOFT完整3,841字节读取，SHA e24c261f5c44157f376c7c1259b779bd77dfe5ef3bbe7eb78dff2ee3c3182d49；只下载公开metadata，无新表达矩阵。原32材料与GSM唯一对应，临床患者别名仍待，不升级患者审阅。
- 11:50:37–11:50:40（19:50:37–19:50:40）：冻结无结局病例/来源和基因机会，终态0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196。
- 11:55:03–11:56:52（19:55:03–19:56:52）：19,018细胞、38发现病例/28保留来源样本、三视图各5,000候选及20PC写完；终态2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656。
- 11:59:59–12:01:06（19:59:59–20:01:06）：另实现从原obs/counts和物理基因轴全复读，精确整数1%、全部表达/变换/控制重放，终态58385bca98d751470d7f7026b0e6ff4ca7f1c584dacbfa6b7eb5ac644f7b5d4b。所有视图发现设计秩6；未训练临床模型。
- 12:04:32（20:04:32）：首次driver因R工作脚本缺失断言退出，模型目录尚不存在。先前Python专用部署工具对R语法编译失败且未写文件，不是科学模型阴性。
- 12:06:31.734393（20:06:31.734393）：直接部署同SHA R源码，经真实R parser通过，旧driver进程已不存在且未生成模型后，在新控制句柄launchfix2启动原不变Python driver。PID392166、startticks314029788；首R子进程392334。原错误日志保留，无科学重算或参数改变。
- 12:09:30（20:09:30）：取证driver同PID/startticks存活、R child存在；首K2主模型RUNNING，terminal尚未生成。随后日志确认已进入多层模型步骤，但不把部分输出当稳定结果。

当前节点状态：输入完整技术接受；全部五个预声明官方拟合已派发，由同一driver顺序执行。完整拟合、稳定性、保留来源投射与生物学解释待登记。直接互作总目标未完成。

## 关键源码及私有凭证

源码精确SHA见[技术checkpoint](../inventories/analysis/2026-10-18/immune_RNA_MCP_discovery_v1/checkpoint_v1.json)。

- `${RUN_ROOT}/immune_RNA_MCP_discovery_scope_20261018_v1/private_terminal.json`
- `${RUN_ROOT}/immune_RNA_MCP_native_expression_inputs_20261018_v1/private_terminal.json`
- `${RUN_ROOT}/immune_RNA_MCP_native_expression_input_audit_20261018_v1/private_terminal.json`
- `${RUN_ROOT}/immune_RNA_MCP_official_models_20261018_v1/plan_before_any_official_MCP_fit.json`
- `${RUN_ROOT}/HNSCC_receiver_followup_job_control_20261016_v2/run_immune_RNA_MCP_official_family_v1_launchfix2/launch_receipt.json`

原计划先四项，第五项来源情境敏感性在任何拟合结果前写进官方模型plan。旧病例scope的字符串不等计数把空字符串与unknown计为差异；不当作真实阶段冲突数量。最新TNBC pre只在新旁表纠正旧空阶段，原CM表不覆盖。

公开仅方法/源码/脱敏执行轨迹和技术尺寸/SHA。学习基因、权重、病例/细胞记录、科学效应和图只服务器及获准汇总对话。无新增受控数据申请、身份标签或原闭合统计家族。


## 官方名称接口的预先核查

12:12 UTC（20:12 UTC+8）在读取学习基因/效应之前检查官方源码，确认make.cell.type会删除视图名的下划线。v1 worker的便捷gene-slot计数仍用原带下划线标签，因此两个T视图可能被错误报告为空；不据此作生物阴性或重跑模型。后续独立结果reader按明确原名→官方规范名一一映射，重算全部gene-slot及官方分数一致性，保留v1便捷摘要和原官方模型文件。此问题在官方模型之后的摘要读取层，不改变正在执行的拟合。


## 2026-10-18T12:18:08.557710+00:00 UTC /2026-10-18T20:18:08.557710+08:00 模型终态门槛与后续合同已冻结

首主要K2拟合12:14:13.214109 UTC（20:14:13.214109 UTC+8）exit0，第二种子12:14:13.230541 UTC真实开始；12:17:12取证仍运行，未查看学习基因/评分或依据部分效果选模型。五模型全终态仍待。

后续合同于12:15:03.156283 UTC（20:15:03.156283 UTC+8）实际冻结，plan SHA1172d887cffb7c3b0c28a21583154fbe3718a1f13d79edb85500343ead5e2b09，terminal SHA04a862f1126bae8290c20f6f6f7a79fbbba6f571a6407abb6254e2767879cffa。全部33投射位置、三视图官方规范名一一连接、全部官方分数重建、K2种子全矩阵匹配及K3/来源敏感性均固定。未来只对主要锚点2程序×3视图对的6个保留来源关联统一BH；不可估计保留，不用P挑seed/K。固定19999残差置换/2000块bootstrap、完整原丰度/数量和质量控制、n>=20及满秩/df>=10门槛；已审患者子集无覆盖则如实不可估计。临床测试不在本合同内。

新来源的sourceCC不能当发现来源：原metadata保持NA，投射仅减去已知控制项，保留不可识别的来源常数；同来源相关/秩不受该平移影响，禁止比较跨来源绝对分数。新来源批次交互、阶段/临床病例映射限制不会被这一步消除。稳定性0.70相关/0.50权重cosine为预声明工程诊断，不当临床或统计验证阈值。

投射Rreader已真实parser通过，Python控制器及另实现全分数/病例均值审核器已AST/服务器编译并逐SHA核对；12:17:37确认投射目录尚不存在。它们是已准备代码，尚未真实执行，不声称数值通过。等五拟合终态再用同一冻结输入，先核官方全33位置、所有supported分数/权重和病例聚合，独立接受后才做稳定性/来源效果。

本轮早期内容39eefe077de35b8b7ec10fc63b70c9b69bf9dfaa已12:13:16.182487 UTC普通push，12:13:22.732918 UTC独立ls-remote/fetch核云，bundle fec9e572df4efce58fa535c4ea4ecbed6e9de9ee9f8916c7e8f5eb457d50b789。完整goal active，科学结果/逐病例行/学习权重只server。
