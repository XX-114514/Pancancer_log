# 历史V5/K9的CoVarNet原文式通讯补分析

事件与取证：2026-11-06，UTC与UTC+8精确时间见各冻结/完成收据；用户显示日期按Australia/Melbourne。范围为20260825 V5/K9全部九CM，不重启注释或CM发现。用户要求实际执行A全局CPDB、B CMT High/Low、C患者稳健性和D Immune Dictionary响应，并公开经审查代码和汇总。

## 初始状态与检查

主本地Pancancer_log有大量并发未提交修改，未改动；在隔离检出fetch最新main并fast-forward，起点2c19567。已读取AGENTS、README、STATUS、ROADMAP、TODO、首页通讯分页、0906证据链、V7 CM/blisa报告、TME临床计划、同步合同和相关run/日志。旧V5身份仍是历史分支，V7及肿瘤限定分支缺本任务可接受的完整冻结扩展资格，未创建混合分支。

## 已执行

- 经原Tailscale跳板读取服务器规则、空闲资源、同用户进程和原始run；无同任务已运行CPDB。没有改代理、停止他人或重提注释。
- 取回开放原文XML、Zenodo v2 Fig3g/3h代码、Immune Dictionary原始S3、MGI同源表、官方CPDB数据库v4.1.0；冻结方法/运行版本与SHA。
- 真数据W/H、网络节点、原生zarr逐细胞索引、全部九CM CMT/成员及样本资格核对；确认H为9×1303、W为87×9，无并列、无非有限；nsNMF theta=0.5且S为9×9。
- 生成71329细胞的固定抽样真实全局输入；每身份最多1000；全部读取块计数域校验通过。
- 20次置换真实CPDB试运行完成。输出SHA、LR-ID对齐后的显著均值与p值一致性通过。71329真实输入结构及充分统计Welch数值核验通过。
- 正式全局1000置换和全样本原生LR/响应统计实际启动；可恢复监督器覆盖全部九CM的High/Low、响应和汇总。最终计算及科学验收状态以run报告后续收据为准。

## 失败与修正

初次受sandbox阻止的SSH经授权执行后恢复；公开资源tar传输曾超时，保留部分文件并按所需文件补传、核验。元数据首版将占位字符串当已知患者，v2加入patient_id_source资格，旧表保留不发布。初次CPDB审计直接按行对齐，而官方输出有不同排序；v2按id_cp_interaction对齐，并按官方p≤0.05规则核对。初次compileall因系统缓存目录权限失败，改用私有tmp缓存通过。未改变科学分组、成员、阈值或选择显著轴。

## 结果边界与后续

患者层在预冻结同来源/癌种/材料/已知时点且两组各10患者的规则下无候选层；不以细胞级P替代。现有交集基因对CPDB2923个互作仅271个满足当前完整亚基映射覆盖，保留不可测背景。旧LIANA只用于历史跨样本复现对照。

待：正式A/B/D终态及独立读回、全图审阅、完整LR和阴性家族解释、源驱动分解及结果发布。C患者级与pseudobulk不可估计原因需保留。代码与汇总公开，患者行、矩阵、完整原始日志留研究环境。

建议commit: Add frozen V5 CoVarNet communication protocol and real-data execution audit.

验证：全仓库validate_repository.sh退出0（既有警告保留）；新增Python AST、JSON和TSV列数检查通过，新增公开文件无私有服务器路径或患者行，git diff --check通过。此次为正式运行中的检查点，不声称B/D结果完成。
