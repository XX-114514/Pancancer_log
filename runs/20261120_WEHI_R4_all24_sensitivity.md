# 全24组参数比较

[完整报告](../inventories/analysis/2026-11-20/WEHI_R4_melanoma_all24/README.md)；[执行日志](../logs/codex/2026/2026-11-20_WEHI-R4-all24-sensitivity.md)。

生产比较与本地独立重放均exit0。生产环境scanpy，绘图matplotlib；复核Python标准库。命令：`python compare_melanoma_N_sensitivity.py`（配置路径后），`python review_melanoma_all24.py ${RESULTS_ROOT}`，`python plot_melanoma_all24.py ${RESULTS_ROOT}`。原始模型未重跑。本轮环境具体版本未另采集；原模型审计与上次归档保留。
