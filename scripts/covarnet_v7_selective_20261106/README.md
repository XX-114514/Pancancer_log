# 执行代码合同

这是本次已执行的科学脚本和配置，字节SHA见`executed_code_manifest.json`；本地聚合绘图单独列入同一清单。输入为既有私有V4注释、来源资格与技术上下文，不随公开仓库提供。请先读[方法](../../methods/CoVarNet-V7-selective-author-CM-20261106.md)与[真实运行](../../runs/20261106_v7_selective_CoVarNet_author_v1.md)。

- `prepare_paper_input.py`/`audit_paper_input.py`：保留的v1来源范围诊断。
- `prepare_paper_input_v2.py`/`audit_paper_input_v2.py`：主v2，加入GSE179994保留矩阵范围暂存，V4标签及其他方法不变。
- `run_author_rank.R`、`assemble_author_modules.R`：固定作者函数及本次配置的秩调查、选秩、CM网络与置换。
- `audit_author_models.R`、`audit_author_network.py`：独立模型/CCC/WSH/置换，以及990对相关/作者校正/特异性和成员边的复算。
- `render_CM_results.py`：仅对已审聚合输出作图/汇总，不修改科学筛选。

Python输入脚本通过`runpy.run_path`注入`PROJECT_ROOT: pathlib.Path`，网络审计注入`RUN_ROOT: pathlib.Path`。R脚本从命令行接收`RUN_ROOT`（逐K脚本再接收K）。配置在每个run保存为`config/paper_config.json`。脚本使用的相对输入版本与本次绑定一致，不能把换输入后运行当成本次复现。

作者源文件取自固定commit`e59b314941e73c8a6eb008ee537cfd461c7e10bd`下`R/CoVarNet.R`，以LF统一换行后放入`${RUN_ROOT}/author/CoVarNet.R`；`author_source_binding.json`声明这一换行规范及实际字节SHA，不冒充Git原始blob字节。这里只抽取需要的函数执行，不需加载作者完整包的全部依赖。兼容性改动及作者P值处理详见方法。

这些脚本会新建输出目录并拒绝覆盖；完整的SSH连接、PID、物理服务器路径与bootstrap不公开。生产流程是有期限的阶段执行，不是长期定时任务。Cloud缺输入时不可宣称真实生产运行成功。
