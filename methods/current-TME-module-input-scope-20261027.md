# 当前TME模块输入范围与旧白名单的对照

来源资格表是有版本的准备清单。旧版曾通过的单位可能被后续原始门控证据排除，也可能因未闭合的关键字段暂缓；新补证可恢复此前未决分支。历史表保留，不能把其文件名中的qualified解释为永久准入。

[对照实现](../scripts/reconcile_TME_scope_manifest_v1.py)把旧单位稳定键连接当前资格工作单位及明确的父级关系，逐个核对来源、完整档案分母和去向；不以显示行号、GSE整体或近似名称作对应。拆分的作者样本保留独立资格单位，不把父级存储单位白名单直接套到全部细胞。当前逐细胞来源覆盖及全部继承覆盖是实际读取入口，物理库、作者样本、原生存储单位和独立患者继续区分。

交付包含当前通过、暂缓、排除本用途、旧到新crosswalk、新补证通过单位及来源汇总。主状态列明确使用当前决定，历史资格另列，避免较旧source_qualification列覆盖新增暂缓。字段以稳定键对齐，表格中不同可选证据列的空缺不改变记录身份。

[独立审核](../scripts/audit_TME_scope_manifest_v1.py)不导入生产帮助函数，重新连接原输入、精确比较所有集合与分母，并按每个旧父级验证分配守恒。独立核查服务器现行清单SHA和原生集合绑定；同步失败或观察超时不冒充远端更新成功，也不触发科学计算重跑。

运行方式：

```bash
python scripts/reconcile_TME_scope_manifest_v1.py \
  --old "${OLD_SOURCE_TABLE}" --review "${CURRENT_SOURCE_REVIEW}" \
  --live "${LIVE_AUTHORITY_RECEIPT}" --output "${NEW_SCOPE_DIRECTORY}"
python scripts/audit_TME_scope_manifest_v1.py \
  --old "${OLD_SOURCE_TABLE}" --review "${CURRENT_SOURCE_REVIEW}" \
  --live "${LIVE_AUTHORITY_RECEIPT}" --delivery "${NEW_SCOPE_DIRECTORY}" \
  --output "${INDEPENDENT_SCOPE_AUDIT}"
```

通过来源资格意味着可继续身份/QC和实际样本准备。模块输入合同固定为`PREPARATION_ONLY_NOT_READY_FOR_MODULE_FIT`，直到最终共同身份/QC、实际组织与技术重复、来源别名及组成覆盖分别冻结。未决身份留在明确分母中，状态与身份分开；患者未知不能当独立患者。旧带日期CM保持原输入，新范围不自动继承旧模块数或编号。

主问题针对未按谱系富集的肿瘤组织scRNA；原发、转移、受累淋巴结及未知部位保留背景层。正常/癌旁、血液、产品、富集/分选混合和snRNA分别研究，不用捕获细胞数恢复人为分选前比例。来源/癌种依赖可以是值得解释的背景，模块稳定性同时评估来源权重与来源内支持，不要求所有癌种产生相同生态。

私有名单、个体标识、来源科学统计和原始附件不进入本公开方法。公开同步限本代码、方法及技术日志。
