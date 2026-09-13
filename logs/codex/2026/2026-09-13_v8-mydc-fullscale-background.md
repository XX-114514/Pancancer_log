# 2026-09-13 V8 Myeloid/DC 全量后台启动

## 完成事项

- 以统一 5M marker 对象替代 54 个异构 raw adapter，验证其与 V7 cell index 完全一致。
- 从全 uncertainty union 中移除不涉及 Myeloid/DC 的冲突，避免 T/B/上皮伪候选。
- 加入 margin gate、identity/state 分离、结构性缺失防护、样本 checkpoint 与原子输出。
- 三轮真实 smoke 中保留前两轮修正证据，最终 smoke v3 2/2 success。
- 通过低资源 preflight 后以 setsid/flock 启动 54-GSE full run。

## 验证与限制

4 Python tests 和 launcher contract 通过。任务为单 worker review-only，不运行 scANVI、Numbat
或重复 scATOMIC 新队列结果。当前 Git 只保存轻量运行记录和 SHA-256，不复制细胞级输出。


## Git 状态同步：2026-09-13 09:03:20 +0800

- 后台 run 仍为 running；93 success、41 explicit skipped、
  1 running、0 failed、981 pending、
  205 pending-too-few-scope。
- PID 559/988 均存活；worker RSS 2998676 KiB。未把运行进度解释为 accuracy。
- 写入 `runs/snapshots/20260913_084700_mydc_full_v1_20260913_090320.json`（SHA-256 `b50ead5b8dcfd8c405c8b1da2631382adb2fab88ffab68f017b4f9f2c65a2859`），随后执行仓库与 release validators。
