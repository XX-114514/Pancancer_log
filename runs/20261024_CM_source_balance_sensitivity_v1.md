# 来源权重敏感性 v1

版本日期使用Australia/Melbourne。真实源执行为2026-10-23T13:22:12.142673+00:00至13:23:49.539449+00:00 UTC，即UTC+8的21:22:12至21:23:49、Melbourne 2026-10-24 00:22:12至00:23:49。取证与日志在此后补录，不将补录时间当冻结时刻。

## 输入和预设

原免疫图谱固定T范围的组成/metadata、原nsNMF和成员文件仅源内存读取。完整路径、输入SHA、固定维度/rank/参数及来源统计都保留在私有合同；原CM标签、图及权威对象不改。来源权重敏感性和外部空间任务是独立有限工作。

合同在真实拟合前冻结，SHA256：`ee0a4e3819c534772480539aa6f0020c121c080a253e5bf91f269e34ca491618`。通用R代码SHA256：`cccce702ea8095f5950a8996e5d7a4fa1d60e5343135f8e41f541819d7d12a3c`。

方法见[说明](../methods/CM-source-balanced-sensitivity-20261024.md)。源R4.3.3、NMF0.28，现有环境，单线程，R地址空间上限8GiB、CPU/墙钟上限1200秒。两权重方案、四个固定种子、每次30次初始化、原2000迭代上限。所有原生/新模型数组仅RAM，允许的源级短汇总单独守卫后留私有；不新增下载、模型/矩阵缓存、临床或P值家族。

可复算命令结构（真实控制合同和结果位置仅私有）：

```bash
Rscript --vanilla scripts/measure_CM_source_balance_RAM_v1.R --check-artificial
Rscript --vanilla scripts/measure_CM_source_balance_RAM_v1.R "${PRIVATE_FROZEN_PLAN}"
```

## 实际执行和验收

源人工KL/匹配/实际nsNMF接口检查actual R exit0；真实有限计算actual R exit0、SSH exit0。所有预设位置完整，原种子模型数值重放、归一化/WSH与标量余弦检查通过；输入末次SHA一致。独立汇总验收核两个总目标与逐来源损失、完整方案/种子位置及迭代记录。所选拟合未触达迭代上限，不能据此宣称全局最优。

私有源实际退出收据SHA256：`e4b25301d28fe19f4a89af0167ac2519468e9df8612255cddcfbae70fdf5b553`。科学相似度、权重相关统计、排序及解读仅获准私有/对话，不写入公共仓库。

同期既有半小时上传heartbeat实际view确认启用；只读stat记录保留私有，没有读取仍在写入的用户RDS。真实外部空间输入/模型仍待完整稳定及全hash，沿用[冻结输入run](20261024_CM02_uploaded_fixed_spatial_preparation_v1.md)。

## 失败、未验证和后续

本新真实计算没有非零退出或重试。既有论文正文认证/CAPTCHA限制继续停止，不以搜索摘要当完整方法；正常原始摘要/图注与作者工具文档各标作用域。不存在本次新增患者/空间/临床或干预证明。

原结果和并行atlas工作保留。技术工作归档和全仓库检查随后补记；公共push先前自动审批拒绝，不重试，云端待同步。goal active。


最终仓库全量检查actual exit0、0 errors/1既有warning，输出SHA256：`ae19d66d087e8ad4b547094098602beff0234bedeb894e4cc23dd6715762e67e`。显式新文件的内部路径/来源ID/凭据模式和科学数值排查、人工diff及git diff --check通过；冻结R SHA保持。只归档本轮文件；公共push旧拒绝保持，云端待同步。
