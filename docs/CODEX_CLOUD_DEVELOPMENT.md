# Codex Cloud开发指南

## 本项目的分工

| 环境 | 承担工作 | 输入/输出 |
|---|---|---|
| Codex Cloud | 代码、配置、合成测试、来源级公开信息整理、方法和任务日志、PR审阅 | 这个GitHub仓库；轻量文件 |
| 数据服务器 | 下载、计数矩阵、逐样本QC/doublet、去重、Zarr、全量CoVarNet/LIANA | 外部数据；运行合同、逐来源分片、审计摘要 |
| 本仓库 | 决策、当前快照、版本与任务交接 | 不包含患者行、大型对象或直接服务器凭据 |

先看[当前报告](../inventories/runtime/2026-10-11/REPORT.md)，再看[任务板](../TODO.md)。Cloud是代码开发入口；这里没有把数据服务器搬入Cloud，也没有声称已配置云端账号环境。

## 首次配置（用户在Codex界面执行）

在 **Work in → Cloud → Select environment → Create environment** 选择 `XX-114514/Pancancer_log`。若列表缺少仓库，先为GitHub连接授权该仓库。Get started中要求环境检查依赖、执行下述检查，然后 **Publish** 环境。当前配置使用Install script/Start skill；不按旧版setup/maintenance教程猜测字段。[官方Cloud环境说明](https://learn.chatgpt.com/docs/environments/cloud-environments)

本仓库工具仅依赖Python标准库和bash，运行要求Python ≥3.10、git、rg及sha256sum。让环境测试执行：

```bash
bash scripts/cloud_check.sh
```

若机器的python3较旧但python满足版本，使用 `PYTHON=python bash scripts/cloud_check.sh`。记录实际解释器版本。所有轻量检查离线执行，不需要OpenAI API凭据、患者表、SSH私钥或GEO大文件。

本机个人skills不会自动同步到Cloud；应将有必要的方法合同写入仓库。Cloud环境有资源限制且不能假定连接私有服务器。环境保存的VM状态不能代替Git提交。[官方Cloud环境说明](https://learn.chatgpt.com/docs/environments/cloud-environments)

## 每个任务如何交接

1. 使用[首个Cloud任务提示](../prompts/cloud_next_task.md)，每次只选择任务板的一项。
2. Cloud读取AGENTS、当前状态和输入合同，修改可移植代码，运行合成测试与仓库检查，记录真实数据执行缺口。
3. 用PR交付代码与运行说明；审阅输入/输出兼容性和科学版本边界后合入main。
4. 服务器操作者检出明确commit，按运行合同先只读预检，检查空间、环境与现有作业，再执行已授权阶段。跨研究整合仍需用户决定。
5. 服务器独立读回验证，导出脱敏汇总；按照[同步说明](RUNTIME_SYNC.md)回写新快照、日志和任务状态。

不要让Cloud直接修改旧冻结版本、通过猜测补患者疗效、重新拟合验证集或把QC分片当新增验收。公共GitHub资料可检索；受控/患者级数据必须保留在服务器，测试使用合成fixtures。

## 建议的下一项

**MAP-001：来源与临床映射验证器**。服务于现有CM来源验证和新临床队列验收，分别处理study/library/sample/patient、组织、时间点、治疗、疗效、生存时间与事件。先建立验证器与缺失/冲突报告；真实患者表服务器执行，不复制到公共PR。完成映射后再推进K稳定性、临床验证与新数据整合。
