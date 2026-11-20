# WEHI R4 全24模型空间敏感性

用户恢复目标并要求继续，沿用公开结果归档授权。开始前旧记录19/24；服务器检查原PID已不存在，24组作者/密度回执齐全，无缺失。未重启模型。精确拟合完成时刻未在本轮逐项提取；本轮取证记录时间为 2026-11-20T01:44:45.937020+00:00 / 2026-11-20T09:44:45.937020+08:00。

读取AGENTS、README/STATUS/ROADMAP/TODO、RUNTIME_SYNC、既有Visium方法、CoVarNet解释及上次归档run/daily/Codex日志。原工作区629项改动保留，使用隔离clone（起点19b2f48）。

执行compare_melanoma_N_sensitivity.py --check后正式运行；scanpy环境，8片×N5/30/60固定设计、means/q05、原始/组成校正/父类比例相关。24组原审计PASS、当前输入SHA和同片坐标全匹配。1440作者分数、7104边、15168差值；本地标准库review_melanoma_all24.py独立重读、逐条核对差值与符号定义通过。plot_melanoma_all24.py出图并目视检查，无遮挡。

修改新增汇总表、图、审计、代码模板、说明及导航。无科学运行失败，无新P值。既有缺失值符号修正在全量读回通过。完整数字与结论边界见[结果](../../../inventories/analysis/2026-11-20/WEHI_R4_melanoma_all24/README.md)。患者桥接、病理图像、全临床外部验证仍未完成，总体目标active。

验证：全部新增文件SHA清单、无原始矩阵/个体列、真实路径脱敏；仓库验证及diff检查后普通commit/push，再独立ls-remote核验。发布回执保存在本地/服务器，未核验前不宣称上传成功。

建议提交：Complete all24 WEHI R4 melanoma N sensitivity review。

发布检查补充：首次diff检查将标准库csv输出CRLF误报为尾随空白；通过目录专属cr-at-eol属性保留已校验原字节，没有改动统计值。新增14内容文件的SHA/大小及真实路径扫描通过。

最终发布检查：validate_repository.sh完成0错误/1项既有路径警告；新增目录路径/凭据扫描无命中，git diff --cached --check通过。首次索引追加在提升权限环境中因默认Python编码失败，改在本地Python执行后同步成功；未影响科学计算。
