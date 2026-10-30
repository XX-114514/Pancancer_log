# 原九来源UMAP受控比较后台启动

用户要求在已保存历史基线后后台启动UMAP重绘。范围沿用[原九来源登记](../../../methods/CM02-original9-integration-baseline-20261029.md)，包括未校正PCA、Harmony、BBKNN、scVI；不恢复临床或LR任务。

开始前检查AGENTS、README、STATUS、ROADMAP、TODO、相关基线方法和日志及Git状态；已有其他TME流程的修改全部保留。读取原九来源Zarr必要轴、计数层属性与少量计数，未回传细胞矩阵或患者资料。源工作区禁止删除，采用全新运行目录并保留失败环境。此前公共push拒绝保持，本轮不重试。

已观察：九来源共同可测、身份无歧义、原冻结HVG有3000个；归档主分析标记筛出2,667,542个候选细胞。逐来源小样本stored counts为非负整数，但完整计数审核仍是后台第一阶段，候选数不是最终共同入模数。该资格不表示原发灶或临床入组。scRNA/snRNA混合及部分assay字段缺失如实保留，不推断临床阶段或协议。

方法：新四臂共用raw counts和特征／细胞轴；PCA用存储全基因轴raw总量归一化后log1p CP10K，全部合格细胞拟合50维稀疏PCA。相对历史ambient-adjusted chosen-count、landmark PCA的差别单独记录，历史对照不能作为纯算法差别。Harmony theta=1、100中心、最多20轮；BBKNN每来源3邻居、9来源、trim=27、PyNNDescent；普通图k=27，另报实际度数。scVI使用NB、gene-batch dispersion、50潜变量、两层128、batch2048、30轮、10轮KL预热、无early stopping。训练前GPU合成检查，保存后重载核对潜变量。UMAP消费各自全量图，同一PCA二维初始化、100轮、min_dist=0.3、固定种子；新Leiden仅数字标签，旧12群作参照。

执行准备：核心环境真实导入scanpy、Harmony和BBKNN通过；一个旧scVI环境因SparseDataset兼容失败，未修改；另一个scVI 1.4.1／Torch 2.5.1环境导入通过。实际GPU训练尚待调度后的合成门禁。源端Slurm可用；CPU任务每个8核／64GiB，scVI训练另申请1 GPU，训练完成后独立CPU作业构图。六阶段依赖可防止输入失败后下游继续。阶段事件、代码SHA、输入轴SHA、真实退出与错误各留新文件；历史输出不覆盖。

验证与限制：本地worker AST通过。运行脚本包含归一化／PCA／Harmony／BBKNN／UMAP合成门禁；真实全量计数、模型、图和UMAP均须各阶段验收。图混合、连通性与旧群保持仅为描述，不能替代原生marker、生物身份和状态保留审阅。100轮UMAP及30轮scVI是本轮固定初始预算，不预先声称收敛或优于历史结果。

首次提交的SSH在banner exchange超时，退出255，没有收到调度凭据。先沿原跳板只读核实目录／调度记录，不能把网络失败记成科学阶段失败，也不盲目重交。

待收尾：记录实际提交凭据、调度状态和观察时间；交付本机启动页。完整结果与新注释仍未完成。建议commit：Launch controlled original-nine UMAP comparison in background。

实际恢复与启动：本机Tailscale状态为Stopped。初次up因要求显式保留非默认参数而拒绝；随后按CLI报告的现有accept-dns=false、accept-routes原值启用，exit0，未重置配置。跳板恢复后只读确认新目录不存在，才重新提交。2026-10-29 10:13:06 UTC／18:13:06 UTC+8收到六个Slurm作业凭据：26412输入与全PCA、26413 PCA图／UMAP、26414 Harmony、26415 BBKNN、26416 scVI GPU训练、26417 scVI图／UMAP。初次队列快照：输入作业PENDING/Priority，其余PENDING/Dependency；已后台提交不表示已开始训练。参数合同SHA256：fdb4f956066dd8acc54b7cde4c9cb82e62a5c30dcf0788da643ff82ea25e0e9f。源端全新目录与阶段依赖均已真实建立。

后续只读取证：2026-10-29 10:13:26 UTC／18:13:26 UTC+8，26412实际RUNNING，其他五项仍PENDING/Dependency，无阶段失败凭据；尚无新UMAP或模型结果。私有启动页、参数、提交状态与文件SHA已生成，本地链接校验通过。启动页为静态观察，不代表自动监控通知。

仓库首轮轻量检查发现本任务日志的相对链接多一级，已按实际目录修正；期间一个本地修正命令因工作目录前缀重复未执行成功，随后使用绝对路径修正。科学参数和后台作业没有改变。

10:14:30 UTC／18:14:30 UTC+8只读复核：worker已有STAGE_STARTED及COUNTS_INPUT_STREAM事件，已顺序审核到253,952／2,670,057归档行，无失败收据。这证明核心合成门禁已实际通过并进入全量输入读取；全量计数接受、PCA完成仍待。新增进度快照单独保留，不改初始启动状态记录。

收尾检查：仓库轻量检查修正后actual exit0，0 error／1 warning；本任务日志和daily的本地相对链接、启动页链接、worker AST、Slurm凭据与RUNNING事件核对通过。git diff --check通过。仅提交本任务两份技术日志；公开push不重试，其他未提交修改保留。
