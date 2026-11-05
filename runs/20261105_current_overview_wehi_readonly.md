# 当前研究总览：WEHI实时只读补核

本记录补充[当前研究总览](../docs/homepage/07_current_work.md)的运行新鲜度，不替代[00:49 UTC系统审计](../docs/reports/20261105_server_analysis_audit/README.md)。用户在Edge School的Open OnDemand终端手动登录后，本次只读查看调度、指定运行的小型日志/配置/完成报告和两个已知Git副本；未提交、取消或重跑科学作业，未修改服务器文件。

## 时间、范围与证据

取证窗口为**2026-11-05 01:50:36—01:52:26 UTC / 09:50:36—09:52:26 UTC+8 / 12:50:36—12:52:26 Melbourne**。不同读取不是原子快照，尤其正在增加的完成报告数不能混作同一时刻。脱敏的手工转录摘要见[结构化收据](../inventories/runtime/2026-11-05/wehi_overview_followup.json)；逻辑路径相对于`${PROJECT_ROOT}`，真实路径、账户与完整作业日志不公开。

使用`squeue`、`scontrol show job`、`sacct -X`；Python标准库只读小JSON、配置、文件列表、日志末尾和SHA256；Git仅`status`、`rev-parse`、`branch`、`rev-list`及`log`。没有加载大矩阵、重算指标或全盘搜索。

## TME run_v2：Leiden已完成，最终对象仍未出现

01:50:36 UTC队列和随后`scontrol`确认作业30660430为RUNNING。作业名称沿用`integrate-tme-v1`，但SubmitLine明确指向`tme_v1/run_v2`、`TME_HARMONY=source:2`，不能按作业名称误判成run_v1。30660431为PENDING，依赖`afterok:30660430`；后处理脚本将执行恶性决定和范围补丁，此时尚未启动。

生产日志保留00:31:34 UTC的`harmony_done`，新增**01:27:57 UTC / 09:27:57 UTC+8的`leiden_done`，78簇**。这只是生产阶段记录，不是最终图的科学验收。01:52:12 UTC指定位置的`tme_v1/run_v2/tme_v1.zarr`和`report/tme_summary.json`均不存在；report目录只有HVG及按来源恶性汇总。尚不能声称新版UMAP、最终对象或后处理完成。运行中的`ExitCode=0:0`也不是最终退出成功。

下一步等待真实终态，再核排除范围、ID、层语义、mask与同图UMAP/marker。此次没有因为缺summary而重复提交。

## V7：K4计算收据已出现，来源敏感性分析正在运行

`Project_v5/runs/20261103_V7_CM_blisa_v1/outputs/covarnet_primary_strict_K4/CoVarNet_execution_report_v1.json`在01:52:26 UTC被实际读取。它声明NMF及作者网络计算完成：K4、47节点、1,207分析观察、53来源、unit权重、2个final seed模型；报告明确不声称独立来源/患者验证。文件mtime为01:19:02 UTC / 09:19:02 UTC+8，这只能定位文件修改时刻，不能替代调度器结束时间。本次按`cm_K4`名称查询accounting未返回记录，未据此推断任务失败或确认其退出码。

这纠正早先“尚无K4完成收据”的当前状态，但不证明K4优于K9。仍需读取成员、活动、拆分/合并、网络与尺度不变匹配；K9历史受控分支保留。

来源等权作业30661687（`primary_strict_K9_source_equal`）和逐来源留出数组30661686正在执行。配置的分析范围仍为历史V7 strict，不是扩展72来源TME。来源等权固定K9，seed为20261103/04，final nrun30、max_iter2000；留出清单53来源，数组0–52、并发上限27，抽查GSE115978配置固定K9、seed20261103、final nrun30、max_iter2000、unit权重。它们不会重新选择一个全局最优K。

01:51:41 UTC读取到25份LOSO完成报告，accounting为25 COMPLETED（退出0:0）、27 RUNNING、1 PENDING（数组并发限额）；01:52:12 UTC再次枚举时已变为26份。26份都声明NMF/作者网络计算完成、K9和52来源；25份为47节点，1份为43节点，均明确不声称独立验证。**节点集合会随留出改变**，后续必须报告覆盖损失及共同节点对应，不能只比较同形W。来源等权此时有中间文件、尚无最终报告；尚未观察到全53来源汇总或科学稳定性验收。

## Git与存储边界

WEHI主代码副本为`main`，HEAD `e6cf388ff5ca53c4e00cd4fcaa7559d72028a71e`，有三个未跟踪入口：`decisions/ADR-0004-v7-covarnet-rank-k9-comparison.md`、`scripts/nmf_benchmark/`、`scripts/v7_cm_blisa/`；本次未改动。独立发布副本HEAD `7617443f90cb42bf4e772d21774731cc9fc020c6`，工作树干净。两者相对于**各自缓存**的origin/main均0/0且无ahead提交；本次未在服务器fetch，不能把0/0解释成与GitHub实时同步。本机联网fetch确认GitHub为`b268012ef7c4678d49ce1117867a4e00a29a2aba`。这些未跟踪入口也不能直接称“从未公开”，须与已归档脚本逐项比较。

本次还确认配置PC_ROOT指向持久项目区，而当前作业SubmitLine仍显式使用scratch树。配置改变不证明运行/环境均已迁移；不改变此前两个数据层的迁移验收，也不替代全资源持久性检查。

## 结论与未检查项

本次改变的是运行状态与计算收据：TME推进到78簇；V7 K4计算有收据；来源等权/LOSO真实运行。CT75主空间验证、V7空间资格、NMF尺度问题、注释身份及临床推广结论均不因这些变化而自动升级。未验收完整W/H、全53来源、注释全量、新TME对象、全部环境和个体级数据。旧审计文件及hash保持原样。
