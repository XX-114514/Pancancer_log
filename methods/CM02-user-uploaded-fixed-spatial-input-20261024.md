# 用户上传单文件作为固定空间分析输入

本阶段只修输入与执行生命周期，沿用[固定状态/空间分析](CM02-external-TDLN-fixed-state-spatial-20261023.md)的全部R、参数、单位/父群准入、缺测和15槽家族。不是增加科学重试机会或改变状态。

## 输入

用户本地下载后自行SFTP上传的单个公开RDS保留原位置。程序只读、核owner/regular file、完整长度、同一描述符/路径inode与读前后元信息。首遍完整MD5/SHA通过后，再从同一描述符作第二遍整文件hash，防止已复制前缀变化仅靠时间戳未被区分。第二遍只用块缓冲，不建立第二个整对象或新原生文件。实际SFTP未完成时不启动科学任务。

[读取器](../scripts/read_uploaded_RDS_RAM_v1.py)不下载、写入、重命名或改变输入权限。[工作器](../scripts/read_external_TDLN_uploaded_workflow_RAM_v1.py)将已核压缩RAM通过原[管道](../scripts/pipe_native_RDS_RAM_v1.py)交给原R；展开/对象/个体状态与空间数组仅RAM，输入文件本身由用户保留。不再把“无新增原生缓存”表述成服务器不存在用户文件。原先已接受的开放目录身份与相同整文件hash绑定，运行时不重复公网请求。

## 执行与输出

[启动器](../scripts/launch_external_TDLN_detached_RAM_v1.py)仅部署小型通用代码/冻结合同，同一控制名重试只返回已有启动标识，部分启动要求审查而不重复起任务。[监督器](../scripts/supervise_external_TDLN_RAM_v1.py)脱离SSH会话，实际wait工作器；输入10分钟、R2小时及120秒监督余量，保留既有父/R内存边界。源进程组在正常/失败/超时后清理，先观察工作器退出再清理仍占stdout的后代，避免误等管道EOF直到总超时。

所有来源帧须通过白名单、大小及个体内容守卫后才持久保存汇总；真实退出、冻结代码/完整输入、R和全15槽数值由[消费者](../scripts/accept_external_TDLN_uploaded_fixed_v1.py)核对。消费者复用原数值/NA/家族校正函数，不放宽统计；观察失败不重启工作器。实际源工作器wait与监督进程本身的OS退出作用域分开。

## 人工验证与边界

[输入反例](../scripts/check_uploaded_RDS_RAM_v1.py)覆盖短文件、hash错、增长、原位修改、路径替换、symlink/目录、预算及原R管道。[汇总反例](../scripts/check_external_TDLN_uploaded_acceptance_v1.py)包括正/负/结构NA、坏hash/代码/家族校正、未完成/坏退出及意外个体字段。[源生命周期](../scripts/check_detached_source_supervisor_v1.py)用独立控制连接核延迟成功、exit7、个体帧拦截、超时与孤立后代；所有夹具为人工，无真实原生或临床输入。

人工测试失败和修复保留，源R管道及修复后五任务均实际接受。原完整统计合成此前已通过，本轮未改变或重跑以挑结果。真实上传仍不完整，人工验证不代表真实单位/状态/空间效果。参见[准备run](../runs/20261024_CM02_uploaded_fixed_spatial_preparation_v1.md)。
