# DIALOGUE完整方法入口已修复，真实MCP发现尚待

## 验证内容

官方固定commit `9c146ccf28d7706aaa60d00947a9126b4e75fd69`，未修改DIALOGUE算法。私有R4.3.3依赖和源码SHA固定；共享环境未变。完整合成流程UTC2026-10-14 00:54:56—00:55:09退出0，00:55:16核实驱动/子进程消失。

两类各720模拟细胞／24模拟样本，两个预设共同因子：PMD→HLM→NNLS最终评分都执行。另一次独立读取保存产物，检查2880个分数有限、两个非空MCP、两类25/24行基因证据以及与已知因子的最小绝对相关.99820/.99846。结果只能证明运行和合成信号恢复，不是临床或生物学验证。

[验收](runtime_acceptance_v1.json)、[独立检查](synthetic_independent_output_checks_v2.tsv)、[依赖](runtime_dependencies_v7.tsv)。真实细胞使用0；未给项目赋新状态/MCP或临床p。

## 已解决故障及保留历史

|尝试|具体结果／修正|
|---|---|
|v1|PMA编译器不在PATH；其他私有依赖安装成功|
|v2|补PATH后官方包安装完成；合成X缺PC列名，官方筛选索引失败|
|v3|列名正确，进入HLM；共享lme41.1-34与Matrix1.6-5接口不匹配，模型错误被tryCatch吞成NA|
|v4|重编译同版本失败，旧lme4源码的Cholmod/Eigen头文件不兼容|
|v5|私有lme41.1-35.5编译成功；路径检查未正确处理新目录及/home→实际路径规范化，误拒绝私有库|
|v6|复用成功安装库及正确规范化检查，PMD/HLM通过；最终NNLS依赖未在官方DESCRIPTION声明，缺包停止|
|v7|加入私有nnls1.6；完整三步通过；独立产物复核通过|

首次独立产物检查把`length(MCPs)`当实际模块数；官方`sig2MCP`默认保留MCP1—5的列表槽位，其中3—5为NULL。修正为检查非空槽位＋官方k=2，数值不改变，失败检查保留。原模拟例完善为两个共同因子及具名PC以覆盖官方接口；没有修改真实数据、临床标签或选择通讯阳性。

上游对兼容问题已有[Matrix接口讨论](https://github.com/lme4/lme4/issues/763)和[相同旧源码编译错误](https://github.com/lme4/lme4/issues/756)。本机通过实测定位，按私有兼容安装解决，不以搜索结果当本机成功证据。

## 复用入口与真实分析边界

[bootstrap](../../../../scripts/DIALOGUE_private_runtime_bootstrap.R)通过`PANCANCER_DIALOGUE_RUNTIME_ROOT`加载私有目录，顺序为library_v7、library_v5、library_v2、library，再使用原R库。每次运行校验私有lme4版本和真实规范路径，避免静默回落至旧ABI。真实根目录留服务器/调用环境，不写入公共报告。

现有5源240533原生单细胞输入v1在末尾把212文库与424 sample-major单位混用而退出；矩阵已经写出，独立全量读回验收正在按[另任务合同](../state_singlecell_readback_recovery_plan_v2.json)恢复。未验收前不把它直接升级为训练输入，也不重新装配同一数据。

下一真实MCP研究保持来源内、患者整块、冻结发现/留出、原编号内变化与类型内PC；原CM01/03的CRC细胞接口是来源限定探索入口。现有数据均参与过项目发现，留出只相对新MCP训练成立，不能冒称全项目独立验证。随后保留其他谱系与全9CM映射、可审临床队列及空间验证；本方法准备不能替代完整目标。
