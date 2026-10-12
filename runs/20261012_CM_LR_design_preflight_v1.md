# Run：原CM全候选完整设计预检

UTC 2026-10-12T05:34:57.723231+00:00 / UTC+8 2026-10-12T13:34:57.723231+08:00，PREDECLARED_NOT_STARTED。父LR v2已独立完整验收；5377数据门槛记录/13497原CM映射，双活动口径与患者块均值诊断，共40491记录。

[合同与全部输入/代码SHA](../inventories/analysis/2026-10-12/CM_LR_design_preflight_plan_v1.json)。主仍6列完整rank/condition/leverage，不自动去协变量；结构常数、必要控制共线与CM曝光自身无独立变化分开。完整nuisance span下的曝光唯一性只做数学诊断，不拟合或接受新主检验。患者块平均只作诊断，不改原主模型单位或临床映射。

producer/auditor依次执行（提供PANCANCER_PREFLIGHT_ROOT、PANCANCER_CM_DESIGN_PLAN），新独立目录。公开诊断标量、不含private samples；没有beta/CI/p或验证轴。真实首尾/exit另记。

## 2026-10-12T13:49:04.609306+08:00 真实完成

两阶段exit0、40491全行独立QR/块重建通过；[完整结果](../inventories/analysis/2026-10-12/CM_LR_design_preflight_v1/RESULTS.md)。主3269严格通过/8491控制恒定/1737高影响点，CM目标系数唯一不等于非零或机制。后续另预声明补充参数化；本轮不估beta/p。
