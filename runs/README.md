# 正式运行记录

本目录保存正式运行的轻量级可复现记录，不保存分析对象或完整程序日志。

创建方式：

```bash
bash scripts/new_run_record.sh "phase07 annotation"
```

记录文件按年份组织，命名为 `YYYYMMDD_HHMMSS_run-name.md`。创建后应：

1. 在运行前填写输入、代码 commit、环境、参数、命令和完成标准；
2. 在运行后填写结束时间、退出状态、输出验证、错误和解释；
3. 在 `inventories/runs.tsv` 增加一行；
4. 重要产物加入 `inventories/artifacts.tsv`；
5. 不复制 `${RUN_ROOT}` 中的完整日志或大型对象。

已有的外部历史运行可先在 inventory 建索引；只有在需要深入审计时才补建对应 Markdown 记录。
