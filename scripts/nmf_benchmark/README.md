# V7 CM NMF比较：运行代码与收据

首次记录：2026-11-04；2026-11-05审计归档。配套[系统报告第6节](../../docs/reports/20261105_server_analysis_audit/README.md#6-nmf比较与本次代码审计新增发现)、[原聚合结果](../../inventories/server_analysis/2026-11-05/tables/nmf/)、[45任务运行信息](../../inventories/server_analysis/2026-11-05/nmf_tasks.tsv)、[源码commit收据](../../inventories/server_analysis/2026-11-05/nmf_source_commits.tsv)。

`build_nmf_envs.sbatch`与`install_nmf_r.R`保存当时环境搭建方法；其中有浮动版本要求，**重执行不会保证复现历史依赖**。复现须先使用已有环境、源commit和原环境收据，不能把新安装HEAD当作原版本。完整环境导出留服务器。

`smoke_nmf.sbatch`运行包级合成输入。新归档版要求`NMF_SMOKE_OUT`指向尚不存在的新目录，保留以往输出。`PC_ROOT`、`CODE`和micromamba module由正常HPC环境配置；MAMBA_ROOT_PREFIX使用数据根下的独立位置，不写硬编码用户路径。

`cm/prep_inputs.py`接受CoVarNet输入目录与新run目录，构造同一X_minmax和counts；`cm/cm_bench.sbatch`按无表头的method/seed/K任务表运行单任务。参数与输出：

```bash
# 从正常HPC终端运行；COV_INPUT和RUN必须明确绑定到相应V7分支/新输出。
python -I scripts/nmf_benchmark/cm/prep_inputs.py "$COV_INPUT" "$RUN"
# tasks.tsv每行method<TAB>seed<TAB>K，本次历史运行45行、K9、5seed。
sbatch --array=0-44 --export=ALL,CODE="$PWD/scripts/nmf_benchmark/cm",PC_ROOT,RUN,TASKS \
  scripts/nmf_benchmark/cm/cm_bench.sbatch
```

R方法worker使用原NMF、RcppML、CoGAPS、RcppPlanc或log1pNMF；Python方法worker使用sklearn或cNMF。`W.tsv`为节点×因子，`H.tsv`含样本维度，只留服务器。归档版拒绝已存在的method/seed目录，防止覆盖；原数值算法未改变。`compare_methods.py`保留原匹配定义，写比较汇总前必须使用新的run目录。

**解释限制：** 原比较按W行份额归一再匹配，不具备因子尺度不变性；见[尺度审计](../server_audit/check_nmf_scaling.py)。K9是受控比较，计数臂及cNMF自身预处理改变输入几何，LIGER按来源分组且有失败起点。原参考nsNMF不是准确性金标准，速度、来源eta²和匹配相似度不能独立证明生物学有效。
